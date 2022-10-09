import rfspice as rf
import numpy as np

# Declare the overall S-matrix
SM1 = np.matrix([[0.11, 0.12], [0.21, 0.22]], dtype=np.complex128)

# Define devices: Index, Number of ports
d1 = rf.Device(1, 2)

# Set devies S-matrices
d1.setSMatrix(SM1)

# Create a network of devices
n1 = rf.Network(1, [d1])

# Define excitation ports: Device, Local Port, Global Port
n1.excitation(d1, 1, 1)
n1.excitation(d1, 2, 2)

# Solve the system
n1.solve(logFlag=True)

# Print results
print(n1.SMatrix)



