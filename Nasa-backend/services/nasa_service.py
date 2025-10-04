# backend/services/nasa_service.py
import requests
import os
from dotenv import load_dotenv

load_dotenv()

NASA_API_KEY = os.getenv("NASA_API_KEY", "DEMO_KEY")
BASE_URL = "https://api.nasa.gov/neo/rest/v1"

def get_neo_feed(start_date, end_date):
    print("nasa is",NASA_API_KEY)
    url = f"{BASE_URL}/feed?start_date={start_date}&end_date={end_date}&api_key={NASA_API_KEY}"
    r = requests.get(url)
    return r.json()

def get_asteroid_details(asteroid_id):
    url = f"{BASE_URL}/neo/{asteroid_id}?api_key={NASA_API_KEY}"
    r = requests.get(url)
    return r.json()
