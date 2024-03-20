import math

G = 6.6743e-11 #The universal gravitational constant (m^3 kg^-1 s^-2).
SOLAR_MASS = 1.9885e30 #The mass of the sun (kg).
SOLAR_RADIUS = 6.957e8 #The volumetric mean radius of the sun (m).
EARTH_RADIUS = 6.378e6 #Radius of the earth, metres
AU = 1.496e11


def stellarMass(stellarRadius, surfaceGravity):
    """
    Arguments:
        stellar radius (float) -- Units: `Solar radii`.

        surface gravity (float) -- Units: `ms^-2`.

    ------------------------
    Returns:
        stellar mass (float) -- Units: `Solar masses`.
    """
    return (surfaceGravity * (stellarRadius*SOLAR_RADIUS)**2)/(G*SOLAR_MASS)

# Calculates the transit impact parameter using the Star Radius, Planet Radius, Orbital Period and Transit Duration
def transitImpactParameter(stellarRadius, planetaryRadius, orbitalPeriod, transitDuration):
    return (((((stellarRadius*SOLAR_RADIUS - planetaryRadius*EARTH_RADIUS)**2)-((semiMajorAxis(orbitalPeriod, transitDuration)*AU*math.sin((transitDuration * math.pi)/(orbitalPeriod*86400)))**2))**0.5)/(stellarRadius*SOLAR_RADIUS))

# Calculates the semi major axis using the Star Mass, Orbital Period and Orbital Radius
def semiMajorAxis(MASS, PERIOD):
    #returned in metres
    return ((G*SOLAR_MASS*MASS*(86400*PERIOD)**2/(4*math.pi**2))**(1/3))/AU

# Calculates the orbital inclination using the Star Radius, Planet Radius, Orbital Period and Transit Duration
def planetOrbitalInclination(starRadius, planetRadius, MASS, PERIOD, transitDuration):
    return math.acos((transitImpactParameter(starRadius,planetRadius, PERIOD, transitDuration)*starRadius*SOLAR_RADIUS)/(semiMajorAxis(MASS, PERIOD)*AU))*(180/math.pi)

def planetaryRadius(solar_radius, flux):
    return (solar_radius*SOLAR_RADIUS*(abs(flux)**(1/2)))/EARTH_RADIUS
