from FF import cableParameters as cp

# Configuration Definitions
cm = 1.0E-2
mils = 2.54E-5
freq = 1.0E+9
epsr = 1.0
tand = 0.0
sigma = 58E+6

d = 50*mils
h = 1*cm
rw = 7.5*mils

config = cp.Config()
config.freq = freq
config.Ns = 200
config.plane.sigma0 = sigma

# Tube 1
# config.N = 4
# wire1 = cp.Wire(1, -(3/2)*d, h, rw, 1.2*rw, epsr, tand, sigma)
# wire2 = cp.Wire(2, -(1/2)*d, h, rw, 1.2*rw, epsr, tand, sigma)
# wire3 = cp.Wire(3, +(1/2)*d, h, rw, 1.2*rw, epsr, tand, sigma)
# wire4 = cp.Wire(4, +(3/2)*d, h, rw, 1.2*rw, epsr, tand, sigma)
# config.wires = [wire1, wire2, wire3, wire4]

# Tube 2 & 3
# config.N = 2
# wire1 = cp.Wire(1, -(1/2)*d, h, rw, 1.2*rw, epsr, tand, sigma)
# wire2 = cp.Wire(2, +(1/2)*d, h, rw, 1.2*rw, epsr, tand, sigma)
# config.wires = [wire1, wire2]

config.check()
print(config) # Optional

# Compute Parameters
cp.computeParameters(config)

