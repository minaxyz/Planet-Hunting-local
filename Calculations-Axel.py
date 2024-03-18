# Calculates Transit Impact Parameter, semi major axis through two methods
# and the orbital inclination of a planet using the Star Radius, Planet Radius, Orbital Period and Transit Duration

import math
def planetaryRadius(stellarRadius, peakRelativeFlux):
    return stellarRadius * math.sqrt(peakRelativeFlux)
    
# Calculates the transit impact parameter using the Star Radius, Planet Radius, Orbital Period and Transit Duration
def transitImpactParameter(stellarRadius,planetaryRadius, orbitalPeriod, transitDuration):
    return (((((stellarRadius-planetaryRadius)**2)-((SemiMajorAxis2(orbitalPeriod, transitDuration)*math.sin((transitDuration * math.pi)/(orbitalPeriod)))**2))**0.5)/(stellarRadius))

# Calculates the semi major axis using the Star Mass, Orbital Period and Orbital Radius
def semiMajorAxis1(StarMass, orbitalPeriod, OrbitalRadius):
    #returned in metres
    return round((6.67*(10**-11)*StarMass*(orbitalPeriod**2))/((2*OrbitalRadius*math.pi)**2),2)

# Calculates the semi major axis using the Orbital Period and Transit Duration
def semiMajorAxis2(orbitalPeriod, transitDuration):
    return round(((6.67*(10**-11)*StarMass*(orbitalPeriod**2))/(4*math.pi))**1/3,2)

# Calculates the orbital inclination using the Star Radius, Planet Radius, Orbital Period and Transit Duration
def planetOrbitalInclination(stellarRadius,planetaryRadius, orbitalPeriod, transitDuration):
    return math.acos((transitImpactParameter(stellarRadius,planetaryRadius, orbitalPeriod, transitDuration)*stellarRadius)/(SemiMajorAxis2(orbitalPeriod, transitDuration)))
