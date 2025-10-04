# backend/services/orbit.py
import math
import numpy as np

AU = 1.496e+11  # meters
MU = 1.32712440018e20  # m^3/s^2, Sun’s GM

def kepler_to_cartesian(a, e, i, Omega, w, M, steps=500):
    """
    Convert Keplerian elements to Cartesian orbit trajectory
    Returns list of positions [x, y, z] in meters
    """
    a = a * AU  # convert AU to meters
    i = math.radians(i)
    Omega = math.radians(Omega)
    w = math.radians(w)
    M = math.radians(M)

    positions = []

    for k in range(steps):
        # Mean anomaly update (simplified)
        Mk = M + 2*math.pi * (k/steps)
        
        # Solve Kepler’s Equation for E (eccentric anomaly)
        E = Mk
        for _ in range(10):  # Newton-Raphson iteration
            E = E - (E - e*math.sin(E) - Mk) / (1 - e*math.cos(E))

        # True anomaly
        theta = 2 * math.atan2(math.sqrt(1+e)*math.sin(E/2),
                               math.sqrt(1-e)*math.cos(E/2))

        # Distance
        r = a * (1 - e*math.cos(E))

        # Position in orbital plane
        x_orb = r * math.cos(theta)
        y_orb = r * math.sin(theta)
        z_orb = 0

        # Rotate to 3D space
        R11 = math.cos(Omega)*math.cos(w) - math.sin(Omega)*math.sin(w)*math.cos(i)
        R12 = -math.cos(Omega)*math.sin(w) - math.sin(Omega)*math.cos(w)*math.cos(i)
        R13 = math.sin(Omega)*math.sin(i)
        R21 = math.sin(Omega)*math.cos(w) + math.cos(Omega)*math.sin(w)*math.cos(i)
        R22 = -math.sin(Omega)*math.sin(w) + math.cos(Omega)*math.cos(w)*math.cos(i)
        R23 = -math.cos(Omega)*math.sin(i)
        R31 = math.sin(w)*math.sin(i)
        R32 = math.cos(w)*math.sin(i)
        R33 = math.cos(i)

        x = R11*x_orb + R12*y_orb + R13*z_orb
        y = R21*x_orb + R22*y_orb + R23*z_orb
        z = R31*x_orb + R32*y_orb + R33*z_orb

        positions.append([x, y, z])

    return positions
