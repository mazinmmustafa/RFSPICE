import rfspice as rf
import numpy as np

# Declare the overall S-matrix
SM1 = np.matrix([[0.111, 0.112], [0.121, 0.122]], dtype=np.complex128)
SM2 = np.matrix([[0.211, 0.212], [0.221, 0.222]], dtype=np.complex128)

# Define devices: Index, Number of ports
d1 = rf.Device(1, 2)
d2 = rf.Device(2, 2)

# Set devies S-matrices
d1.setSMatrix(SM1)
d2.setSMatrix(SM2)

# Create a network of devices
n1 = rf.Network(1, [d1, d2])

# Add connections between ports
n1.connect(d1, 2, d2, 1)

# Define excitation ports: Device, Local Port, Global Port
n1.excitation(d1, 1, 1)
n1.excitation(d2, 2, 2)

# Solve the system
n1.solve(logFlag=True)

# Print results
print(n1.SMatrix)



