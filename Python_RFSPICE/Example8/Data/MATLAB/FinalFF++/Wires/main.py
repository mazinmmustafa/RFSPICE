from FF import cableParameters as cp

# Configuration Definitions
mm = 1.0E-3
freq = 1.0E+9
epsr = 4.0
tand = 1.0E-3
sigma = 58E+6

config = cp.Config()
config.freq = freq
config.N = 3
config.Ns = 50

wire0 = cp.Wire(0, +0.0*mm, +0.0*mm, 1.0*mm, 2.0*mm, epsr, tand, sigma)
wire1 = cp.Wire(1, -5.0*mm, +10.0*mm, 1.0*mm, 2.0*mm, epsr, tand, sigma)
wire2 = cp.Wire(2, +5.0*mm, +10.0*mm, 1.0*mm, 2.0*mm, epsr, tand, sigma)
config.wires = [wire0, wire1, wire2]
config.check()
print(config) # Optional

# Compute Parameters
cp.computeParameters(config)

