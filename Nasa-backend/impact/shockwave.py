import math

# -------------------------------
# Function to calculate asteroid impact
# -------------------------------
def asteroid_impact(diameter_m, velocity_kms, density=3000, burst_altitude_m=0):
    """
    Parameters:
    diameter_m       : Diameter of asteroid in meters
    velocity_kms     : Velocity in km/s
    density          : Density in kg/m^3 (default 3000)
    burst_altitude_m : Burst altitude in meters (default 0 for ground impact)
    
    Returns:
    Dictionary with mass, kinetic energy, TNT equivalent, thermal energy, shockwave radii
    """
    
    # --- Step 1: Mass of asteroid ---
    radius = diameter_m / 2
    mass = density * (4/3) * math.pi * radius**3  # kg
    
    # --- Step 2: Kinetic energy ---
    velocity_ms = velocity_kms * 1000  # convert km/s to m/s
    kinetic_energy_j = 0.5 * mass * velocity_ms**2  # Joules
    
    # TNT equivalent
    tnt_kg = kinetic_energy_j / 4.184e6
    tnt_kt = tnt_kg / 1e3  # kilotons
    
    # --- Step 3: Thermal energy ---
    thermal_fraction = 0.35
    thermal_energy_j = kinetic_energy_j * thermal_fraction
    thermal_energy_tnt_kg = thermal_energy_j / 4.184e6
    thermal_energy_tnt_kt = thermal_energy_tnt_kg / 1e3
    
    # --- Step 4: Shockwave radius calculation ---
    # Using Hopkinson–Cranz scaled distance: R = Z * W^(1/3)
    # W = TNT mass in kg
    # Z depends on overpressure (typical values)
    overpressure_Z = {
        "0.15_psi_loud_boom": 55,   # approx m/kg^(1/3)
        "1_psi_window_break": 40,
        "3_psi_moderate_damage": 20,
        "5_psi_heavy_damage": 15,
        "10_psi_severe_structural": 10
    }
    
    shockwave_radii_m = {}
    for key, Z in overpressure_Z.items():
        R = Z * (tnt_kg ** (1/3))
        # Apply burst altitude scaling
        if burst_altitude_m > 0:
            R *= 1 / (1 + burst_altitude_m / 10000)  # rough correction factor
        shockwave_radii_m[key] = R  # meters
    
    # Convert to km
    shockwave_radii_km = {k: v/1000 for k, v in shockwave_radii_m.items()}
    
    # --- Step 5: Return results ---
    return {
        "mass_kg": mass,
        "kinetic_energy_j": kinetic_energy_j,
        "tnt_equivalent_kt": tnt_kt,
        "thermal_energy_j": thermal_energy_j,
        "thermal_energy_tnt_kt": thermal_energy_tnt_kt,
        "shockwave_radius_km": shockwave_radii_km
    }

# -------------------------------
# Example usage
# -------------------------------
if __name__ == "__main__":
    # Example asteroid parameters
    diameter = 20        # meters
    velocity = 20        # km/s
    density = 3000       # kg/m³
    burst_altitude = 0   # meters (ground impact)
    
    impact_result = asteroid_impact(diameter, velocity, density, burst_altitude)
    
    print("\n--- Asteroid Impact Results ---")
    print(f"Mass: {impact_result['mass_kg']:.2e} kg")
    print(f"Kinetic Energy: {impact_result['kinetic_energy_j']:.2e} J")
    print(f"TNT Equivalent: {impact_result['tnt_equivalent_kt']:.2f} kt")
    print(f"Thermal Energy: {impact_result['thermal_energy_j']:.2e} J ({impact_result['thermal_energy_tnt_kt']:.2f} kt TNT)")
    
    print("\nShockwave Radii (approximate) for different overpressures:")
    for key, radius_km in impact_result['shockwave_radius_km'].items():
        print(f"{key}: {radius_km:.2f} km")
