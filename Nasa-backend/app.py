# backend/app.py
from flask import Flask, request, jsonify
from services.nasa_service import get_neo_feed, get_asteroid_details
from services.usgs_service import get_recent_earthquakes, get_tsunami_alerts
from services.simulation import simulate_impact
from services.orbit import kepler_to_cartesian
from flask_cors import CORS
import numpy as np
import requests
from services.impact_service import simulate_impact_neo


app = Flask(__name__)
CORS(app)

AU_TO_KM = 149597870.7
EARTH_RADIUS_KM = 6371.0


NASA_API_KEY = "assxgDZway7KraOMkma6Ex3BSFIdDhOIvFoWr4a2"
NEO_WS_URL = "https://api.nasa.gov/neo/rest/v1/neo/{}"


@app.route("/")
def home():
    return jsonify({"message": "Asteroid Defense Backend Running 🚀"})

# ===== NASA NEO API =====
@app.route("/api/neo/feed", methods=["GET"])
def neo_feed():
    start_date = request.args.get("start_date", "2025-10-01")
    end_date = request.args.get("end_date", "2025-10-07")
    data = get_neo_feed(start_date, end_date)
    return jsonify(data)



@app.route("/api/neo/<asteroid_id>", methods=["GET"])
def neo_details(asteroid_id):
    data = get_asteroid_details(asteroid_id)
    return jsonify(data)

# ===== USGS APIs =====
@app.route("/api/earthquakes", methods=["GET"])
def earthquakes():
    data = get_recent_earthquakes()
    return jsonify(data)

@app.route("/api/tsunami", methods=["GET"])
def tsunami():
    data = get_tsunami_alerts()
    return jsonify(data)

@app.route("/api/simulate_impact", methods=["POST"])
def simulate_impact_endpoint():
    """
    Expects JSON:
    {
        "neo_id": "3542519",
        "n_samples": 500,
        "window_days": 365,
        "apply_dv_magnitude_ms": 5
    }
    """
    payload = request.json
    neo_id = payload.get("neo_id")
    n_samples = payload.get("n_samples", 500)
    window_days = payload.get("window_days", 365)
    dv_ms = payload.get("apply_dv_magnitude_ms", 0.0)

    if not neo_id:
        return jsonify({"error": "neo_id is required"}), 400

    results = simulate_impact_neo(
        api_key=NASA_API_KEY,
        neo_id=neo_id,
        n_samples=n_samples,
        window_days=window_days,
        dv_ms=dv_ms
    )

    return jsonify(results)
# ===== Simulation Engine =====
@app.route("/api/simulate", methods=["POST"])
def simulate():
    payload = request.json
    results = simulate_impact(
        size=payload["size"],
        velocity=payload["velocity"],
        density=payload.get("density", 3000),
        impact_angle=payload.get("impact_angle", 45),
        lat=payload.get("lat", 0),
        lon=payload.get("lon", 0)
    )
    return jsonify(results)

@app.route("/api/orbit", methods=["POST"])
def orbit():
    data = request.json
    trajectory = kepler_to_cartesian(
        a=data["a"], e=data["e"], i=data["i"],
        Omega=data["Omega"], w=data["w"], M=data["M"],
        steps=data.get("steps", 500)
    )
    return jsonify({"trajectory": trajectory})

# @app.route("/asteroid/<asteroid_id>", methods=["GET"])
# def get_asteroid_data(asteroid_id):
#     try:
#         # Fetch asteroid details from NASA NEO API
#         url = f"https://api.nasa.gov/neo/rest/v1/neo/{asteroid_id}?api_key={NASA_API_KEY}"
#         response = requests.get(url)
#         data = response.json()

#         if "orbital_data" not in data:
#             return jsonify({"error": "No orbital data found"}), 404

#         orbital_data = data["orbital_data"]

#         # Extract only the required orbital elements
#         result = {
#             "id": data.get("id"),
#             "name": data.get("name"),
#             "semi_major_axis": orbital_data.get("semi_major_axis"),
#             "eccentricity": orbital_data.get("eccentricity"),
#             "inclination": orbital_data.get("inclination"),
#             "ascending_node_longitude": orbital_data.get("ascending_node_longitude"),
#             "perihelion_argument": orbital_data.get("perihelion_argument"),
#             "mean_anomaly": orbital_data.get("mean_anomaly"),
#             "epoch": orbital_data.get("epoch_osculation"),
#         }

#         return jsonify(result)

#     except Exception as e:
#         return jsonify({"error": str(e)}), 500
    

@app.route('/asteroid/<asteroid_id>', methods=['GET'])
def get_asteroid_data(asteroid_id):
    """
    Fetch asteroid orbital + close approach data from NASA NeoWs
    and return cleaned JSON for frontend animation.
    """

    # Step 1: Fetch asteroid data from NASA
    print("api is called")
    response = requests.get(NEO_WS_URL.format(asteroid_id), params={"api_key": NASA_API_KEY})
    print(response.json)
    
    if response.status_code != 200:
        return jsonify({"error": "Failed to fetch data from NASA API"}), 500

    data = response.json()

    # Step 2: Extract orbital elements (for propagation if needed)
    orbital_data = data.get("orbital_data", {})
    elements = {
        "epoch": orbital_data.get("epoch_osculation"),
        "a_au": orbital_data.get("semi_major_axis"),
        "e": orbital_data.get("eccentricity"),
        "i_deg": orbital_data.get("inclination"),
        "raan_deg": orbital_data.get("ascending_node_longitude"),
        "arg_peri_deg": orbital_data.get("perihelion_argument"),
        "mean_anomaly_deg": orbital_data.get("mean_anomaly"),
        "orbit_class": orbital_data.get("orbit_class", {}).get("orbit_class_description")
    }

    # Step 3: Extract close approach info (next Earth encounter)
    close_approaches = data.get("close_approach_data", [])
    close_approach = None
    if close_approaches:
        ca = close_approaches[0]
        close_approach = {
            "date": ca.get("close_approach_date_full"),
            "relative_velocity_km_s": ca.get("relative_velocity", {}).get("kilometers_per_second"),
            "miss_distance_km": ca.get("miss_distance", {}).get("kilometers"),
            "orbiting_body": ca.get("orbiting_body")
        }

    # Step 4: Extract physical properties
    est_diameter = data.get("estimated_diameter", {}).get("kilometers", {})
    physical = {
        "absolute_magnitude_H": data.get("absolute_magnitude_h"),
        "diameter_min_km": est_diameter.get("estimated_diameter_min"),
        "diameter_max_km": est_diameter.get("estimated_diameter_max"),
        "hazardous": data.get("is_potentially_hazardous_asteroid")
    }

    # Step 5: Final JSON for frontend
    asteroid_json = {
        "id": data.get("id"),
        "name": data.get("name"),
        "orbital_elements": elements,
        "close_approach": close_approach,
        "physical": physical
    }

    return jsonify(asteroid_json)


def propagate_orbit(a, e, i, Omega, omega, M, v, n_samples=1000):
    """
    Simple propagation with Monte Carlo sampling to estimate impact lat/lon range.
    """
    # Convert angles to radians
    i = np.radians(i)
    Omega = np.radians(Omega)
    omega = np.radians(omega)
    
    lat_list = []
    lon_list = []

    for _ in range(n_samples):
        # Slight random variation for Monte Carlo
        a_var = a + np.random.normal(0, 0.001)
        e_var = e + np.random.normal(0, 0.001)
        i_var = i + np.random.normal(0, 0.001)

        # Simple Keplerian approximation for distance r
        theta = np.random.uniform(0, 2*np.pi)  # true anomaly
        r = a_var * (1 - e_var**2) / (1 + e_var * np.cos(theta))
        
        # 3D position in Earth-centered frame (simplified)
        x = r * (np.cos(Omega)*np.cos(theta + omega) - np.sin(Omega)*np.sin(theta + omega)*np.cos(i_var))
        y = r * (np.sin(Omega)*np.cos(theta + omega) + np.cos(Omega)*np.sin(theta + omega)*np.cos(i_var))
        z = r * (np.sin(theta + omega) * np.sin(i_var))
        
        # Convert to lat/lon on Earth (assuming distance = 1 AU ~ impact)
        lat = np.degrees(np.arcsin(z / r))
        lon = np.degrees(np.arctan2(y, x))
        
        # Randomly shift for Earth rotation (simplified)
        lon += np.random.uniform(-10, 10)
        
        lat_list.append(lat)
        lon_list.append(lon)

    lat_min, lat_max = min(lat_list), max(lat_list)
    lon_min, lon_max = min(lon_list), max(lon_list)

    return {
        "lat_range": [lat_min, lat_max],
        "lon_range": [lon_min, lon_max],
        "sample_count": n_samples
    }

@app.route('/calculate_impact', methods=['POST'])
def calculate_impact():
    data = request.get_json()
    a = float(data.get('a'))
    e = float(data.get('e'))
    i = float(data.get('i'))
    Omega = float(data.get('Omega'))
    omega = float(data.get('omega'))
    M = float(data.get('M'))
    v = float(data.get('v'))

    result = propagate_orbit(a, e, i, Omega, omega, M, v, n_samples=1000)
    return jsonify(result)





if __name__ == "__main__":
    app.run(debug=True)
