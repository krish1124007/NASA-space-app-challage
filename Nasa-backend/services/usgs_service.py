# backend/services/usgs_service.py
import requests

# USGS Earthquake API
USGS_QUAKE_URL = "https://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/all_day.geojson"
# NOAA Tsunami API
NOAA_TSUNAMI_URL = "https://www.weather.gov/wwamap/json/tsunami"

def get_recent_earthquakes():
    r = requests.get(USGS_QUAKE_URL)
    return r.json()

def get_tsunami_alerts():
    r = requests.get(NOAA_TSUNAMI_URL)
    return r.json()
