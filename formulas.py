import math

G = 6.6743e-11 #The universal gravitational constant (m^3 kg^-1 s^-2).
SOLAR_MASS = 1.9885e30 #The volumetric mean radius of the sun (m).
SOLAR_RADIUS = 6.957e8 #The mass of the sun (kg).
EARTH_RADIUS = 6.378e6 #The CONVERSION of metres to earth radii, to divide the radius in metres by.

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
    return (((((stellarRadius - planetaryRadius)**2)-((semiMajorAxis2(orbitalPeriod, transitDuration)*math.sin((transitDuration * math.pi)/(orbitalPeriod)))**2))**0.5)/(stellarRadius))

# Calculates the semi major axis using the Star Mass, Orbital Period and Orbital Radius
def semiMajorAxis1(stellarMass, orbitalPeriod, orbitalRadius):
    #returned in metres
    return round((6.67*(10**-11)*stellarMass*(orbitalPeriod**2))/((2*orbitalRadius*math.pi)**2),2)

# Calculates the semi major axis using the Orbital Period and Transit Duration
def semiMajorAxis2(orbitalPeriod, transitDuration):
    return round((2*math.pi*transitDuration)/orbitalPeriod, 2)

# Calculates the orbital inclination using the Star Radius, Planet Radius, Orbital Period and Transit Duration
def planetOrbitalInclination(starRadius, planetRadius, orbitalPeriod, transitDuration):
    return math.acos((transitImpactParameter(starRadius,planetRadius, orbitalPeriod, transitDuration)*starRadius)/(semiMajorAxis2(orbitalPeriod, transitDuration)))

def PlanetaryRadius(solar_radius, flux):
    return (solar_radius*SOLAR_RADIUS*(abs(flux)**(1/2)))/EARTH_RADIUS
