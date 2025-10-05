from flask import Flask, jsonify, request
import requests
import os

app = Flask(__name__)

# NASA API Key (set in environment variable, else use DEMO_KEY)
NASA_API_KEY = os.getenv("NASA_API_KEY", "assxgDZway7KraOMkma6Ex3BSFIdDhOIvFoWr4a2")
NEO_WS_URL = "https://api.nasa.gov/neo/rest/v1/neo/{}"

@app.route('/asteroid/<asteroid_id>', methods=['GET'])
def get_asteroid_data(asteroid_id):
    """
    Fetch asteroid orbital + close approach data from NASA NeoWs
    and return cleaned JSON for frontend animation.
    """

    # Step 1: Fetch asteroid data from NASA
    print("api is cooled")
    response = requests.get(NEO_WS_URL.format(asteroid_id), params={"api_key": NASA_API_KEY})
    
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


if __name__ == "__main__":
    app.run(debug=True)
