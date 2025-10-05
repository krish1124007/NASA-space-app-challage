import requests

def get_similar_earthquakes(magnitude):
    url = f"https://earthquake.usgs.gov/fdsnws/event/1/query?format=geojson&minmagnitude={magnitude}&maxmagnitude={magnitude+0.1}"
    response = requests.get(url)
    data = response.json()
    
    # Extract relevant information
    earthquakes = []
    for feature in data['features']:
        mag = feature['properties']['mag']
        place = feature['properties']['place']
        time = feature['properties']['time']
        earthquakes.append({'magnitude': mag, 'location': place, 'time': time})
    
    return earthquakes

# Example usage
magnitude = 7.5  # Replace with your calculated magnitude
earthquakes = get_similar_earthquakes(magnitude)
for eq in earthquakes:
    print(f"Mag: {eq['magnitude']}, Location: {eq['location']}, Time: {eq['time']}")
