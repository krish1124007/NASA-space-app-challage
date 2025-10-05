import math
import requests
import folium

# -------------------------------
# Calculate asteroid impact parameters
# -------------------------------
def asteroid_impact(diameter_m, velocity_kms, density=3000, burst_altitude_m=0):
    radius = diameter_m / 2
    mass = density * (4 / 3) * math.pi * radius ** 3
    velocity_ms = velocity_kms * 1000
    kinetic_energy_j = 0.5 * mass * velocity_ms ** 2
    tnt_kg = kinetic_energy_j / 4.184e6
    tnt_kt = tnt_kg / 1e3
    thermal_fraction = 0.35
    thermal_energy_j = kinetic_energy_j * thermal_fraction
    thermal_energy_tnt_kt = thermal_energy_j / 4.184e9

    # Shockwave radius for different overpressure
    overpressure_Z = {
        "0.15_psi_Loud_Boom": 55,
        "1_psi_Window_Break": 40,
        "3_psi_Moderate_Damage": 20,
        "5_psi_Heavy_Damage": 15,
        "10_psi_Severe_Structural": 10
    }

    shockwave_radii_km = {}
    for key, Z in overpressure_Z.items():
        R = Z * (tnt_kg ** (1 / 3))
        if burst_altitude_m > 0:
            R *= 1 / (1 + burst_altitude_m / 10000)
        shockwave_radii_km[key] = R / 1000

    return {
        "mass_kg": mass,
        "kinetic_energy_j": kinetic_energy_j,
        "tnt_equivalent_kt": tnt_kt,
        "thermal_energy_j": thermal_energy_j,
        "thermal_energy_tnt_kt": thermal_energy_tnt_kt,
        "shockwave_radius_km": shockwave_radii_km
    }

# -------------------------------
# Convert impact to earthquake magnitude
# -------------------------------
def impact_to_magnitude(kinetic_energy_j, seismic_fraction=0.01):
    seismic_energy_j = kinetic_energy_j * seismic_fraction
    magnitude = (2 / 3) * math.log10(seismic_energy_j) - 10.7
    return magnitude

# -------------------------------
# Fetch similar earthquakes from USGS
# -------------------------------
def get_similar_earthquakes(magnitude):
    url = f"https://earthquake.usgs.gov/fdsnws/event/1/query?format=geojson&minmagnitude={magnitude}&maxmagnitude={magnitude+0.2}&limit=10"
    response = requests.get(url)
    data = response.json()
    earthquakes = []
    for feature in data.get('features', []):
        mag = feature['properties']['mag']
        place = feature['properties']['place']
        coords = feature['geometry']['coordinates']
        earthquakes.append({
            'magnitude': mag,
            'location': place,
            'longitude': coords[0],
            'latitude': coords[1]
        })
    return earthquakes

# -------------------------------
# Plot map with all attractive visuals
# -------------------------------
def plot_impact(lat, lon, impact_result, earthquakes):
    m = folium.Map(location=[lat, lon], zoom_start=4, tiles='cartodb dark_matter')

    # --- Impact point (Red) ---
    folium.CircleMarker(
        location=[lat, lon],
        radius=10,
        color='red',
        fill=True,
        fill_color='red',
        popup="🚀 Asteroid Impact Point"
    ).add_to(m)

    # --- Thermal energy radius (Yellow translucent) ---
    thermal_radius_km = (impact_result['thermal_energy_tnt_kt'] ** (1 / 3)) * 5
    folium.Circle(
        location=[lat, lon],
        radius=thermal_radius_km * 1000,
        color='yellow',
        fill=True,
        fill_opacity=0.25,
        popup=f"🔥 Thermal Zone (~{thermal_radius_km:.2f} km)"
    ).add_to(m)

    # --- Shockwave zones (Orange to Red gradients) ---
    shockwave_colors = {
        "0.15_psi_Loud_Boom": "#ffe699",
        "1_psi_Window_Break": "#ffcc66",
        "3_psi_Moderate_Damage": "#ff9933",
        "5_psi_Heavy_Damage": "#ff6600",
        "10_psi_Severe_Structural": "#ff3300"
    }

    for key, radius_km in impact_result["shockwave_radius_km"].items():
        folium.Circle(
            location=[lat, lon],
            radius=radius_km * 1000,
            color=shockwave_colors[key],
            fill=False,
            popup=f"💥 {key.replace('_', ' ')} - {radius_km:.1f} km"
        ).add_to(m)

    # --- Earthquake zones (Blue shades) ---
    eq_colors = ['#99ccff', '#6699ff', '#3366ff', '#0033cc']
    for idx, eq in enumerate(earthquakes):
        folium.CircleMarker(
            location=[eq['latitude'], eq['longitude']],
            radius=7,
            color=eq_colors[idx % len(eq_colors)],
            fill=True,
            fill_opacity=0.6,
            popup=f"🌍 Earthquake<br>Magnitude: {eq['magnitude']}<br>Region: {eq['location']}"
        ).add_to(m)

    # --- Legend ---
    legend_html = """
    <div style="
        position: fixed;
        bottom: 30px;
        left: 30px;
        width: 270px;
        background: rgba(0,0,0,0.75);
        color: #fff;
        padding: 15px;
        border-radius: 10px;
        z-index:9999;
        font-size: 14px;
        box-shadow: 0 0 10px #000;
    ">
        <h4 style='text-align:center; color:#FFD700;'>🌍 Impact Visualization Guide</h4>
        <b style='color:red;'>●</b> Impact Point<br>
        <b style='color:yellow;'>●</b> Thermal Heat Zone<br>
        <b style='color:#ff9933;'>●</b> Shockwave Range (varied shades)<br>
        <b style='color:#3366ff;'>●</b> Earthquake Epicenters<br>
        <hr style="border:1px solid #fff;">
        <i>Hover or click zones for details.<br>Each color represents intensity.</i>
    </div>
    """
    m.get_root().html.add_child(folium.Element(legend_html))

    # --- Floating Title Banner ---
    title_html = """
    <div style="
        position: fixed;
        top: 20px;
        left: 50%;
        transform: translateX(-50%);
        background: rgba(255,255,255,0.1);
        color: #FFD700;
        font-size: 24px;
        font-weight: bold;
        padding: 10px 25px;
        border-radius: 15px;
        text-align: center;
        box-shadow: 0px 0px 15px rgba(255,255,255,0.3);
        z-index: 9999;
    ">
        ☄️ Asteroid Impact Simulation Map ☄️
    </div>
    """
    m.get_root().html.add_child(folium.Element(title_html))

    # --- Layer Control ---
    folium.LayerControl().add_to(m)
    m.save("asteroid_impact_attractive_map.html")
    print("✅ Map saved as 'asteroid_impact_attractive_map.html'")

# -------------------------------
# Main script
# -------------------------------
if __name__ == "__main__":
    print("🪐 Asteroid Impact Simulator - Powered by Folium + USGS\n")

    diameter = float(input("Enter Asteroid Diameter (meters): "))
    velocity = float(input("Enter Impact Velocity (km/s): "))
    density = float(input("Enter Density (kg/m³) [default=3000]: ") or 3000)
    burst_altitude = float(input("Enter Burst Altitude (m) [default=0]: ") or 0)
    lat = float(input("Impact Latitude: "))
    lon = float(input("Impact Longitude: "))

    result = asteroid_impact(diameter, velocity, density, burst_altitude)
    magnitude = impact_to_magnitude(result["kinetic_energy_j"])
    quakes = get_similar_earthquakes(magnitude)

    print("\n--- 🌋 Impact Summary ---")
    print(f"Mass: {result['mass_kg']:.2e} kg")
    print(f"Kinetic Energy: {result['kinetic_energy_j']:.2e} J")
    print(f"TNT Equivalent: {result['tnt_equivalent_kt']:.2f} kt")
    print(f"Thermal Energy: {result['thermal_energy_tnt_kt']:.2f} kt TNT")
    print(f"Predicted Earthquake Magnitude: {magnitude:.2f}")

    print("\nShockwave Radii (km):")
    for k, v in result["shockwave_radius_km"].items():
        print(f" - {k}: {v:.2f} km")

    print("\nFetching Similar Earthquakes from USGS...")
    for eq in quakes:
        print(f"  🌎 {eq['location']} - M{eq['magnitude']}")

    plot_impact(lat, lon, result, quakes)
