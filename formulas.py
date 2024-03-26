import math

G = 6.6743e-11 #The universal gravitational constant (m^3 kg^-1 s^-2).
SOLAR_MASS = 1.9885e30 #The mass of the sun (kg).
SOLAR_RADIUS = 6.957e8 #The volumetric mean radius of the sun (m).
EARTH_RADIUS = 6.378e6 #Radius of the earth, metres
AU = 1.496e11

def stellarMass(stellarRadius, surfaceGravity):
    """
    Arguments:
        stellar radius (float) - The mean volumetric radius of the star `Solar Radii` -- Units: `Solar radii`.

        surface gravity (float) - The acceleration due to gravity on the surface of the star -- Units: `ms^-2`.

    ------------------------
    Returns:
        stellar mass (float) -- Units: `Solar masses`.
    """
    return (surfaceGravity * (stellarRadius*SOLAR_RADIUS)**2)/(G*SOLAR_MASS)

def transitImpactParameter(stellarRadius, stellarMass, planetaryRadius, orbitalPeriod, transitDuration):
    """
    Arguments:
        stellar radius (float) - The mean volumetric radius of the star -- Units: `Solar radii`.

        planetary radius (float) - The radius of the planet -- Units: `ms^-2`.

        semi major axis (float) - The farthest distance between the planet and the star -- Units: `AU`.

        orbital period (float) - The time interval between transits -- Units: `Days`.

        transit duration (float) - The time the planet spends partially blocking the star's light -- Units: `Days`.
    
    ------------------------
    Returns:
        transit impact parameter (float) - The perpendicular distance between the orbit the star's centre, expressed as a ratio of 
        the star's radius -- Units: `Stellar radius ratio`.
    """
    return min((max(((stellarRadius*SOLAR_RADIUS + planetaryRadius*EARTH_RADIUS)**2)-((semiMajorAxis(stellarMass, orbitalPeriod)*AU*math.sin(transitDuration * math.pi/orbitalPeriod))**2), 0)**0.5)/(stellarRadius*SOLAR_RADIUS), 1 + (planetaryRadius*EARTH_RADIUS/stellarRadius*SOLAR_RADIUS))

def semiMajorAxis(stellarMass, orbitalPeriod):
    """
    Arguments:
        stellar mass (float) - The mass of the star -- Units: `Solar radii`.

        orbital period (float) - The time interval between transits -- Units: `Days`.

    ------------------------
    Returns:
        semi major axis (float) - The farthest distance between the planet and the star -- Units: `AU`.
    """
    return ((G*SOLAR_MASS*stellarMass*(86400*orbitalPeriod)**2/(4*math.pi**2))**(1/3))/AU

def planetOrbitalInclination(stellarRadius, semiMajorAxis, transitImpactParameter):
    """
    Arguments:
        stellar radius (float) - The mean volumetric radius of the star -- Units: `Solar Radii`.

        semi major axis (float) - The farthest distance between the planet and the star -- Units: `AU`.

        transit impact parameter (float) - The perpendicular distance between the orbit the star's centre, expressed as a ratio of 
        the star's radius -- Units: `Stellar radius ratio`.

    ------------------------
    Returns:
        planet orbital inclination (float) - The angle between the star's and the planet's orbiting plane -- Units: `degrees`.
    """
    return math.acos((transitImpactParameter*stellarRadius*SOLAR_RADIUS)/(semiMajorAxis*AU))*(180/math.pi)

def planetaryRadius(stellarRadius, peakFlux):
    """
    Arguments:
        stellar radius (float) - The mean volumetric radius of the star -- Units: `Solar Radii`.

        peak flux (float) - The greatest change in flux during a transit -- Units: `Normalised stellar flux`.

    ------------------------
    Returns:
        planetary radius (float) - The radius of the planet -- Units: `ms^-2`.
    """
    return (stellarRadius*SOLAR_RADIUS*(abs(peakFlux)**(1/2)))/EARTH_RADIUS