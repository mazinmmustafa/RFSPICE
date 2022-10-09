import rfspice as rf
import numpy as np

# Declare the overall S-matrix
SM1 = np.matrix([[0.111, 0.112], [0.121, 0.122]], dtype=np.complex128)
SM2 = np.matrix([[0.211, 0.212], [0.221, 0.222]], dtype=np.complex128)
SM3 = np.matrix([[0.311, 0.312], [0.321, 0.322]], dtype=np.complex128)
SM4 = np.matrix([[-1/3, +2/3, +2/3], [+2/3, -1/3, +2/3], [+2/3, +2/3, -1/3]], dtype=np.complex128)

# Define devices: Index, Number of ports
d1 = rf.Device(1, 2)
d2 = rf.Device(2, 2)
d3 = rf.Device(3, 2)
d4 = rf.Device(4, 3)

# Set devies S-matrices
d1.setSMatrix(SM1)
d2.setSMatrix(SM2)
d3.setSMatrix(SM3)
d4.setSMatrix(SM4)

# Create a network of devices
n1 = rf.Network(1, [d1, d2, d3, d4])

# Add connections between ports
n1.connect(d1, 2, d4, 1)
n1.connect(d2, 1, d4, 2)
n1.connect(d3, 1, d4, 3)

# Define excitation ports: Device, Local Port, Global Port
n1.excitation(d1, 1, 1)
n1.excitation(d2, 2, 2)
n1.excitation(d3, 2, 3)

# Solve the system
n1.solve(logFlag=True)

# Print results
print(n1.SMatrix)



