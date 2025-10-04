from flask import Flask, jsonify, request
from flask_cors import CORS  # Add this import
import requests

app = Flask(__name__)
CORS(app)  # Add this line to enable CORS for all routes

@app.route("/")
def home():
    start_date = request.args.get('start_date', '2025-10-01')
    end_date = request.args.get('end_date', '2025-10-05')
    
    url = f"https://api.nasa.gov/neo/rest/v1/feed?start_date={start_date}&end_date={end_date}&api_key=NrrTMw0Od7hKSsd2JeDaZ42pIvwF4xMM2kCShfSX"
    
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raises an HTTPError for bad responses
        data = response.json()
        return jsonify(data)
    except requests.exceptions.RequestException as e:
        return jsonify({"error": str(e)}), 500

@app.route("/usingbyid")
def getobject():
    query = request.args.get("query")
    object_id = request.args.get("id")
    try:
        url = f"https://api.nasa.gov/neo/rest/v1/neo/{object_id}?api_key=NrrTMw0Od7hKSsd2JeDaZ42pIvwF4xMM2kCShfSX"
        fetch_data = requests.get(url)
        fetch_data.raise_for_status()
        data = fetch_data.json()
        return jsonify(data)
    except Exception as e:
        return jsonify({"error": str(e)}), 500

if __name__ == "__main__":
    app.run(debug=True)