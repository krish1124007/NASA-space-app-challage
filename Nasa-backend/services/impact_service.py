from neo_deflect_analysis import neo_to_orbit_from_elements, run_ensemble_and_get_impacts, fetch_neo_by_id
from astropy.time import Time
from astropy import units as u

def simulate_impact_neo(api_key, neo_id, n_samples=500, window_days=365, dv_ms=0.0):
    """
    Run NEO deflection/impact Monte Carlo simulation using your neo_deflect_analysis.py logic
    """
    # 1️⃣ Fetch NEO orbital elements
    data = fetch_neo_by_id(api_key, neo_id)
    elem = data.get("orbital_data")
    if not elem:
        return {"error": "No orbital_data found for NEO"}

    # 2️⃣ Convert to Orbit
    orbit_nom = neo_to_orbit_from_elements(elem)
    t0 = orbit_nom.epoch
    t_window_start = t0
    t_window_end = t0 + (window_days * u.day)

    # 3️⃣ Delta-v setup
    apply_dv = None
    if dv_ms > 0:
        v_dir = (orbit_nom.v / orbit_nom.v.norm()).to(u.dimensionless_unscaled).value
        dv_kms = (dv_ms / 1000.0) * v_dir
        apply_dv = {"t_deflect": t0, "delta_v": dv_kms}

    # 4️⃣ Run simulation
    hits, stats = run_ensemble_and_get_impacts(
        orbit_nom,
        t_window_start,
        t_window_end,
        n_samples=n_samples,
        apply_dv=apply_dv
    )

    # 5️⃣ Prepare response
    response = {
        "hits_count": stats["hits"],
        "n_samples": stats["n_samples"],
        "impact_probability": stats["fraction"],
        "impact_points": hits[:10],  # sample top 10 points
    }
    return response
