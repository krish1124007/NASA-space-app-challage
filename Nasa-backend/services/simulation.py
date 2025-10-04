# backend/services/simulation.py
import math

def impact_energy(diameter, velocity, density=3000):
    radius = diameter / 2
    volume = (4/3) * math.pi * (radius**3)
    mass = density * volume
    return 0.5 * mass * velocity**2

def crater_size(energy):
    return 1.161 * (energy ** 0.294) / 1000  # km

def seismic_magnitude(energy):
    return (2/3) * (math.log10(energy) - 4.8)

def tsunami_height(energy, distance=500):
    return min(100, (energy / 1e15) / (distance**2))

def simulate_impact(size, velocity, density=3000, impact_angle=45, lat=0, lon=0):
    energy = impact_energy(size, velocity, density)
    crater = crater_size(energy)
    magnitude = seismic_magnitude(energy)
    tsunami = tsunami_height(energy, 500)

    return {
        "impact_energy_joules": energy,
        "crater_diameter_km": crater,
        "seismic_magnitude": magnitude,
        "tsunami_height_m": tsunami,
        "impact_location": {"lat": lat, "lon": lon}
    }
