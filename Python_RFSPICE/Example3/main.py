import rfspice as rf
import numpy as np

# Declare the overall S-matrix
SM1 = np.matrix([[0.11, 0.12], [0.21, 0.22]], dtype=np.complex128)
SM2 = np.matrix([[-1.0]], dtype=np.complex128)

# Define devices: Index, Number of ports
d1 = rf.Device(1, 2)
d2 = rf.Device(2, 1)

# Set devies S-matrices
d1.setSMatrix(SM1)
d2.setSMatrix(SM2)

# Create a network of devices
n1 = rf.Network(1, [d1, d2])

# Add connections between ports
n1.connect(d1, 2, d2, 1)

# Define excitation ports: Device, Local Port, Global Port
n1.excitation(d1, 1, 1)

# Solve the system
n1.solve(logFlag=True)

# Print results
print(n1.SMatrix)



