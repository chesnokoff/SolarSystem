import math
from vpython import *

# Data in units according to the International System of Units
G = 6.67 * math.pow(10, -11)

# Mass of the Earth
ME = 5.973 * math.pow(10, 24)
# Mass of the Moon
MM = 7.347 * math.pow(10, 22)
# Mass of the Sun
MS = 1.989 * math.pow(10, 30)
# Radius Earth-Moon
REM = 384400000
# Radius Sun-Earth
RSE = 149600000000
# Force Earth-Moon
FEM = G * (ME * MM) / math.pow(REM, 2)
# Force Earth-Sun
FES = G * (MS * ME) / math.pow(RSE, 2)

# Angular velocity of the Moon with respect to the Earth (rad/s)
wM = math.sqrt(FEM / (MM * REM))
# Velocity v of the Moon (m/s)
vM = wM * REM
print("Angular velocity of the Moon with respect to the Earth: ", wM, " rad/s")
print("Velocity v of the Moon: ", vM / 1000, " km/s")

# Angular velocity of the Earth with respect to the Sun(rad/s)
wE = math.sqrt(FES / (ME * RSE))
# Velocity v of the Earth (m/s)
vE = wE * RSE
print("Angular velocity of the Earth with respect to the Sun: ", wE, " rad/s")
print("Velocity v of the Earth: ", vE / 1000, " km/s")

# Initial angular position
theta0 = 0


# Position at each time
def positionMoon(t):
    theta = theta0 + wM * t
    return theta


def positionEarth(t):
    theta = theta0 + wE * t
    return theta


def fromDaysToS(d):
    s = d * 24 * 60 * 60
    return s


def fromStoDays(s):
    d = s / 60 / 60 / 24
    return d


def fromDaysToh(d):
    h = d * 24
    return h


# Graphical parameters
print("\nSimulation Earth-Moon-Sun motion\n")
days = 365
seconds = fromDaysToS(days)
print("Days: ", days)
print("Seconds: ", seconds)

v = vector(0.5, 0, 0)
E = sphere(pos=vector(3, 0, 0), color=color.cyan, radius=.3, make_trail=True)
M = sphere(pos=E.pos + v, color=color.white, radius=0.1, make_trail=True)
S = sphere(pos=vector(0, 0, 0), color=color.yellow, radius=1)

t = 0
thetaTerra1 = 0
dt = 5000
dthetaE = positionEarth(t + dt) - positionEarth(t)
dthetaM = positionMoon(t + dt) - positionMoon(t)
print("delta t:", dt, "seconds. Days:", fromStoDays(dt), "hours:",
      fromDaysToh(fromStoDays(dt)), sep=" ")
print("Variation angular position of the Earth:", dthetaE,
      "rad/s that's to say", degrees(dthetaE), "degrees", sep=" ")
print("Variation angular position of the Moon:", dthetaM,
      "rad/s that's to say", degrees(dthetaM), "degrees", sep=" ")

while t < seconds:
    rate(50)
    thetaEarth = positionEarth(t + dt) - positionEarth(t)
    thetaMoon = positionMoon(t + dt) - positionMoon(t)

    # Rotation only around z axis (0,0,1)
    E.pos = rotate(E.pos, angle=thetaEarth, axis=(0, 0, 1))
    v = rotate(v, angle=thetaMoon, axis=(0, 0, 1))
    M.pos = E.pos + v
    t += dt