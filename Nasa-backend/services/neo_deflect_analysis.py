"""
neo_deflect_analysis.py

Requirements: numpy, scipy, matplotlib, requests, astropy, poliastro, pyproj, shapely, tqdm, pandas
Install: pip install numpy scipy matplotlib requests astropy poliastro pyproj shapely tqdm pandas

What it does:
- Fetch NEO orbital elements from NASA NEO API (or accept input).
- Build a poliastro Orbit object.
- Run Monte Carlo ensemble with small perturbations to initial state.
- Propagate each sample and detect impacts with Earth (R_earth + h_atm).
- Return impact lat/lon list, compute footprint ranges and impact probability.
- Optionally apply a Δv impulse at t_deflect and re-run to compute the effect.

Notes:
- This uses two-body + optional simple perturbations. For precise long-term orbits include planetary perturbations (poliastro can do some).
- If NASA returns no covariance, we assume an observational uncertainty (tune `sigma_pos`, `sigma_vel`).
"""

import argparse
import requests
import numpy as np
from scipy.optimize import brentq
from scipy.stats import multivariate_normal
from tqdm import tqdm
from datetime import datetime, timedelta
import math
import warnings
warnings.filterwarnings("ignore")

# Astro libraries
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import CartesianRepresentation, CartesianDifferential
from astropy.coordinates import ICRS, GCRS, EarthLocation
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.constants import G, M_earth, R_earth
from poliastro.bodies import Sun, Earth
from poliastro.twobody.orbit import Orbit
from poliastro.twobody.propagation import propagate
from poliastro.twobody.propagation import func_twobody
from poliastro.twobody.propagation import mean_motion
from poliastro.twobody import Ephem
from poliastro.twobody.propagation import propagate
from poliastro.twobody.events import NodeEvent
from poliastro.core.iod import lambert
from poliastro import iod
from poliastro.util import time_range
from poliastro.constants import J2000
import matplotlib.pyplot as plt
import pandas as pd

# For coordinate transforms to geodetic
from pyproj import Transformer
from shapely.geometry import Point, box
from collections import Counter

# ---------------------------
# Utility functions
# ---------------------------
def fetch_neo_by_id(api_key: str, neo_id: str):
    """Fetch NEO lookup from NASA NEO feed by id."""
    url = f"https://api.nasa.gov/neo/rest/v1/neo/{neo_id}"
    params = {"api_key": api_key}
    r = requests.get(url, params=params, timeout=30)
    r.raise_for_status()
    return r.json()

def neo_to_orbit_from_elements(elem):
    """
    Build poliastro Orbit from classical orbital elements dictionary.
    Expect keys: a (AU), e, i (deg), asc_node (deg), arg_peri (deg), mean_anomaly (deg), epoch (ISO)
    The NASA API returns 'orbital_data' with these keys but names differ — map accordingly.
    """
    a = float(elem["semi_major_axis"]) * u.AU
    e = float(elem["eccentricity"]) * u.one
    inc = float(elem["inclination"]) * u.deg
    raan = float(elem["ascending_node_longitude"]) * u.deg
    argp = float(elem["perihelion_argument"]) * u.deg
    M = float(elem["mean_anomaly"]) * u.deg
    epoch_str = elem.get("epoch_osculation") or elem.get("epoch") or elem.get("epoch_jd")
    # Convert epoch to astropy Time
    epoch = Time(epoch_str, scale="utc")
    # Create Orbit using classical elements (from mean anomaly M)
    return Orbit.from_classical(Sun, a, e, inc, raan, argp, M, epoch)

def state_vector_from_orbit(orbit: Orbit):
    """Return (r_vec, v_vec) as astropy Quantity in ICRS frame at orbit.epoch"""
    r = orbit.r.to(u.km)
    v = orbit.v.to(u.km / u.s)
    return r, v, orbit.epoch

def propagate_orbit_rv(r0, v0, t0, times):
    """
    Propagate using poliastro's Orbit.from_vectors and orbit.propagate for each time.
    times: list of astropy.Time objects
    Returns list of r vectors (astropy Quantity km) in ICRS.
    """
    orb = Orbit.from_vectors(Sun, r0, v0, epoch=t0)
    res = []
    for t in times:
        tof = (t - t0).to(u.s)
        r_new, v_new = orb.propagate(tof).rv()
        res.append((r_new.to(u.km).value, v_new.to(u.km / u.s).value, t))
    return res

def earth_state_at_times(times):
    """
    Simplified: get Earth heliocentric state using poliastro Ephem (built-in)
    Return list of (r_earth_km, v_earth_km_s, time)
    """
    from poliastro.ephem import Ephem
    from poliastro.bodies import Earth
    from poliastro.ephem import build_ephem_interpolant
    # Use poliastro's Earth ephemeris via 'Earth' object & JPL approximations
    res = []
    for t in times:
        # get Earth position in Sun-centered frame via Orbit.from_body_ephem
        from poliastro.ephem import get_body_ephem
        r_e, v_e = get_body_ephem("earth", t)
        res.append((r_e.to(u.km).value, v_e.to(u.km / u.s).value, t))
    return res

def distance_vec(a, b):
    a = np.array(a); b = np.array(b)
    return np.linalg.norm(a - b)

def find_impact_time_between(sample_r_times, earth_r_times, Rlimit_km):
    """
    sample_r_times: list of tuples (r_vec_km, v_vec_km_s, time)
    earth_r_times: same times for Earth
    linear scan for crossing; if bracket found refine with interpolation + root finding
    returns (t_impact, r_impact_icrs) or (None, None)
    """
    n = len(sample_r_times)
    for i in range(n-1):
        r1 = np.array(sample_r_times[i][0])
        r2 = np.array(sample_r_times[i+1][0])
        re1 = np.array(earth_r_times[i][0])
        re2 = np.array(earth_r_times[i+1][0])
        d1 = np.linalg.norm(r1 - re1)
        d2 = np.linalg.norm(r2 - re2)
        if d1 > Rlimit_km and d2 <= Rlimit_km:
            # bracket found between t1 and t2; refine by scalar root on distance-Rlimit
            t1 = sample_r_times[i][2].jd
            t2 = sample_r_times[i+1][2].jd
            def f(t_jd):
                # linear interpolation of vectors (sufficiently small steps)
                alpha = (t_jd - t1) / (t2 - t1)
                r_interp = r1*(1-alpha) + r2*alpha
                re_interp = re1*(1-alpha) + re2*alpha
                return np.linalg.norm(r_interp - re_interp) - Rlimit_km
            try:
                t_root = brentq(f, t1, t2)
                alpha = (t_root - t1) / (t2 - t1)
                r_root = r1*(1-alpha) + r2*alpha
                return Time(t_root, format='jd', scale='utc'), r_root
            except Exception as e:
                # fallback coarse location
                t_est = sample_r_times[i+1][2]
                return t_est, sample_r_times[i+1][0]
    return None, None

def inertial_to_geodetic(latlon_transformer, r_icrs_km, t_astropy_time):
    """
    Convert ICRS (heliocentric?) position of asteroid at impact to ECEF geodetic lat/lon.
    Here we assume r_icrs_km is asteroid position wrt Earth center in ITRS already.
    BUT our r vectors are heliocentric; the impact vector we compute is asteroid - earth -> vector in Sun-centered coords.
    After computing r_impact_rel = r_ast - r_earth (sun-centered), it's already Earth-centered inertial (ECI).
    We need to account for Earth's rotation at T to go ECI -> ECEF (ITRS).
    We'll do a simple rotation by GMST (approx) to align frames (sufficient for research demo).
    """
    # Compute Greenwich sidereal angle (approx) in degrees
    # astropy can perform transforms more accurately; use astropy.gcrs/itrs
    from astropy.coordinates import ITRS, CartesianRepresentation
    # Construct a GCRS/ICRS coordinate and transform to ITRS (ECEF)
    # r_icrs_km relative to Earth center => create ICRS coordinate then transform
    cart = CartesianRepresentation(r_icrs_km * u.km)
    # Use GCRS with obstime and then transform to ITRS
    gcrs = GCRS(cart, obstime=t_astropy_time)
    itrs = gcrs.transform_to(ITRS(obstime=t_astropy_time))
    # Extract x,y,z meters
    x = itrs.cartesian.x.to(u.m).value
    y = itrs.cartesian.y.to(u.m).value
    z = itrs.cartesian.z.to(u.m).value
    # Use pyproj to convert ECEF (meters) to lat/lon
    lon, lat, alt = latlon_transformer.transform(x, y, z, radians=False)
    return lat, lon, alt

# ---------------------------
# Main analysis functions
# ---------------------------
def run_ensemble_and_get_impacts(orbit_nominal: Orbit,
                                 t_window_start: Time,
                                 t_window_end: Time,
                                 n_samples=2000,
                                 sigma_pos_km=10.0,
                                 sigma_vel_kms=0.001,
                                 dt_minutes=60,
                                 h_atm_km=100.0,
                                 apply_dv=None):
    """
    Run Monte Carlo ensemble:
    - orbit_nominal: poliastro Orbit object (heliocentric) at epoch
    - t_window_start, t_window_end: astropy Time
    - n_samples: ensemble size
    - sigma_pos_km, sigma_vel_kms: assumed 1-sigma uncertainties to sample around nominal state
    - dt_minutes: time step for coarse propagation
    - h_atm_km: atmosphere edge radius offset
    - apply_dv: None or dict {'t_deflect': Time, 'delta_v': [vx,vy,vz] km/s} in inertial frame
    Returns:
      hits: list of dicts {'time':Time, 'lat':float, 'lon':float, 'r_ecef':..., 'impact_energy_J':...}
      stats: dict with fraction_hit etc.
    """
    # get nominal state
    r0_q, v0_q = orbit_nominal.r, orbit_nominal.v
    t0 = orbit_nominal.epoch

    # convert to numpy arrays in km and km/s
    r0 = r0_q.to(u.km).value
    v0 = v0_q.to(u.km / u.s).value

    # build time grid
    total_hours = (t_window_end - t_window_start).to(u.hour).value
    steps = max(2, int(total_hours * 60 / dt_minutes))
    times = time_range(t_window_start, end=t_window_end, periods=steps)
    # Get Earth states at times (heliocentric)
    # Note: poliastro provides body vectors via get_body_ephem; use that
    from poliastro.ephem import get_body_ephem
    earth_states = [get_body_ephem("earth", t) for t in times]
    earth_r_times = [ (st[0].to(u.km).value, st[1].to(u.km/u.s).value, t) for st,t in zip(earth_states,times) ]

    Rlimit_km = R_earth.to(u.km).value + h_atm_km

    # Transformer for ECEF -> latlon: pyproj expects (x,y,z) in meters and returns lon,lat,alt
    transformer = Transformer.from_crs("epsg:4978", "epsg:4326", always_xy=True)

    hits = []
    hits_count = 0

    # Pre-generate gaussian samples in state space
    # Simple block-diagonal covariance in (x,y,z,vx,vy,vz)
    cov = np.diag([sigma_pos_km**2]*3 + [sigma_vel_kms**2]*3)
    mvn = multivariate_normal(mean=np.concatenate([r0, v0]), cov=cov)

    # sample
    # To speed up: we can sample in chunks and propagate; here sequential for clarity
    for s in tqdm(range(n_samples), desc="MC samples"):
        sample = mvn.rvs()
        r_s = sample[:3]
        v_s = sample[3:]

        # If applying dv at some time before window, modify velocity accordingly
        if apply_dv:
            # If deflection at t_deflect = orbit epoch or a given time, we approximate delta-v applied directly to sample initial vel if t_deflect == t0.
            # For general deflection time, a more accurate approach: propagate to t_deflect, apply dv, repropagate.
            t_def = apply_dv.get("t_deflect")
            dv = np.array(apply_dv.get("delta_v"))  # km/s
            if t_def is None or t_def == t0:
                v_s = v_s + dv

        # build orbit object from vectors
        try:
            orb_s = Orbit.from_vectors(Sun, r_s * u.km, v_s * (u.km / u.s), epoch=t0)
        except Exception as e:
            # skip invalid sample
            continue

        # propagate sample along times, collect r vectors
        sample_r_times = []
        for t in times:
            tof = (t - t0).to(u.s)
            try:
                rv = orb_s.propagate(tof)
                rvec = rv.r.to(u.km).value
                vvec = rv.v.to(u.km / u.s).value
            except Exception as e:
                # fallback: linear approx
                dt = (t - t0).to(u.s).value
                rvec = r_s + v_s * dt
                vvec = v_s
            sample_r_times.append((rvec, vvec, t))

        # Search for impact
        t_imp, r_imp = find_impact_time_between(sample_r_times, earth_r_times, Rlimit_km)
        if t_imp is not None:
            # compute Earth position at t_imp to get relative vector
            # linear interp of earth states
            # find nearest index
            # compute r_ast - r_earth
            # careful: get Earth position at t_imp using get_body_ephem
            earth_r, earth_v = get_body_ephem("earth", t_imp)
            earth_r_km = earth_r.to(u.km).value
            rel = np.array(r_imp) - earth_r_km  # vector in km in GCRS/ICRS approx
            # convert to geodetic lat/lon
            lat, lon, alt = inertial_to_geodetic(transformer, rel, t_imp)
            hits.append({'time': t_imp.iso, 'lat': lat, 'lon': lon, 'alt_m': alt, 'r_rel_km': rel})
            hits_count += 1

    fraction = hits_count / n_samples
    stats = {'n_samples': n_samples, 'hits': hits_count, 'fraction': fraction}
    return hits, stats

# ---------------------------
# CLI and orchestration
# ---------------------------
def main_cli():
    parser = argparse.ArgumentParser(description="NEO deflection impact analysis")
    parser.add_argument("--api_key", required=True, help="NASA API key")
    parser.add_argument("--neo_id", required=True, help="NEO id from NASA")
    parser.add_argument("--n_samples", type=int, default=2000)
    parser.add_argument("--window_days", type=float, default=365.0, help="look ahead window in days")
    parser.add_argument("--sigma_pos_km", type=float, default=20.0)
    parser.add_argument("--sigma_vel_kms", type=float, default=0.001)
    parser.add_argument("--dt_minutes", type=float, default=60.0)
    parser.add_argument("--h_atm_km", type=float, default=100.0)
    parser.add_argument("--apply_dv_magnitude_ms", type=float, default=0.0, help="m/s delta-v magnitude to apply along velocity at epoch")
    args = parser.parse_args()

    data = fetch_neo_by_id(args.api_key, args.neo_id)
    elem = data.get('orbital_data')
    if elem is None:
        print("No orbital_data found for NEO. Exiting.")
        return
    print("NEO name:", data.get("name"))
    orbit_nom = neo_to_orbit_from_elements(elem)
    print("Nominal epoch:", orbit_nom.epoch.iso)
    # set window
    t0 = orbit_nom.epoch
    t_window_start = t0
    t_window_end = t0 + (args.window_days * u.day)

    apply_dv = None
    if args.apply_dv_magnitude_ms > 0:
        # build dv along current velocity direction (unit vector)
        v_dir = (orbit_nom.v / orbit_nom.v.norm()).to(u.dimensionless_unscaled).value
        dv_kms = (args.apply_dv_magnitude_ms / 1000.0) * np.array(v_dir)
        apply_dv = {'t_deflect': t0, 'delta_v': dv_kms}

    hits, stats = run_ensemble_and_get_impacts(orbit_nom,
                                               t_window_start,
                                               t_window_end,
                                               n_samples=args.n_samples,
                                               sigma_pos_km=args.sigma_pos_km,
                                               sigma_vel_kms=args.sigma_vel_kms,
                                               dt_minutes=args.dt_minutes,
                                               h_atm_km=args.h_atm_km,
                                               apply_dv=apply_dv)
    print("Results:", stats)
    if stats['hits'] == 0:
        print("No impacts in ensemble. Impact probability ~ 0 for assumed uncertainties.")
    else:
        # create dataframe and simple summary ranges
        df = pd.DataFrame(hits)
        df['lat'] = df['lat'].astype(float)
        df['lon'] = df['lon'].astype(float)
        print("Impact points (sample):")
        print(df[['time','lat','lon']].head())
        # compute bounding box
        latmin, latmax = df['lat'].min(), df['lat'].max()
        lonmin, lonmax = df['lon'].min(), df['lon'].max()
        print(f"Impact footprint lat range: {latmin:.3f} to {latmax:.3f}")
        print(f"Impact footprint lon range: {lonmin:.3f} to {lonmax:.3f}")
        print(f"Estimated impact probability: {stats['fraction']*100:.4f} % based on ensemble and assumed uncertainties.")

        # cluster coarse subregions by 1-degree grid
        df['lat_r'] = df['lat'].round(1)  # 0.1 deg bins
        df['lon_r'] = df['lon'].round(1)
        counts = df.groupby(['lat_r','lon_r']).size().reset_index(name='cnt')
        counts = counts.sort_values('cnt', ascending=False).head(10)
        print("Top candidate subregions (lat,lon bins) with relative chance in hits:")
        total_hits = stats['hits']
        for _, row in counts.iterrows():
            print(f"  bin {row['lat_r']:.1f}, {row['lon_r']:.1f} -> {row['cnt']} hits -> {(row['cnt']/total_hits)*100:.2f}% of impacts")

if __name__ == "__main__":
    main_cli()
