# backend/app.py
from flask import Flask, request, jsonify
from services.nasa_service import get_neo_feed, get_asteroid_details
from services.usgs_service import get_recent_earthquakes, get_tsunami_alerts
from services.simulation import simulate_impact
from services.orbit import kepler_to_cartesian

app = Flask(__name__)

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


if __name__ == "__main__":
    app.run(debug=True)
