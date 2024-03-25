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

def transitImpactParameter(stellarRadius, planetaryRadius, orbitalPeriod, transitDuration):
    """
    Arguments:
        stellar radius (float) -- Units: `Solar radii`.

        planetary radius (float) -- Units: `ms^-2`.

        orbital period (float) -- Units: `Days`.

        transit duration (float) -- Units: `Days`.
    
    ------------------------
    Returns:
        transit impact parameter (float) -- Units: `Stellar radius ratio`.
    """
    return (((((stellarRadius*SOLAR_RADIUS - planetaryRadius*EARTH_RADIUS)**2)-((semiMajorAxis(orbitalPeriod, transitDuration)*AU*math.sin((transitDuration * math.pi)/(orbitalPeriod)))**2))**0.5)/(stellarRadius*SOLAR_RADIUS))

def semiMajorAxis(stellarMass, orbitalPeriod):
    """
    Arguments:
        stellar mass (float) -- Units: `Solar radii`.

        orbital period (float) -- Units: `Days`.

    ------------------------
    Returns:
        semi major axis (float) -- Units: `AU`.
    """
    return ((G*SOLAR_MASS*stellarMass*(86400*orbitalPeriod)**2/(4*math.pi**2))**(1/3))/AU

def planetOrbitalInclination(starRadius, planetRadius, stellarMass, orbitalPeriod, transitDuration):
    """
    Arguments:
        stellar radius (float) -- Units: `Solar radii`.

        planetary radius (float) -- Units: `Earth radii`.

        orbital period (float) -- Units: `Days`.

        transit duration (float) -- Units: `Days`.

    ------------------------
    Returns:
        planet orbital inclination (float) -- Units: `Degrees`.
    """
    return math.acos((transitImpactParameter(starRadius,planetRadius, orbitalPeriod, transitDuration)*starRadius*SOLAR_RADIUS)/(semiMajorAxis(stellarMass, orbitalPeriod)*AU))*(180/math.pi)

def planetaryRadius(solarRadius, flux):
    """
    Arguments:
        planetary radius (float) -- Units: `Earth radii`.

        flux (float) -- Units: `Normalised stellar flux`.

    ------------------------
    Returns:
        planetary radius (float) -- Units: `Earth radii`.
    """
    return (solarRadius*SOLAR_RADIUS*(abs(flux)**(1/2)))/EARTH_RADIUS
