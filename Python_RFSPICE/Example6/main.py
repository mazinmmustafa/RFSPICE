import rfspice as rf
import numpy as np

# Define frequency vector
Z0 = 50.0
Ns = 401
freq_min = 1.0E+9
freq_max = 2.0E+9
freq = np.linspace(freq_min, freq_max, Ns)
R = 60.0
L = 5.0E-9
C = 8.0E-12
G = 40.0E-3

# Define devices
TL1 = rf.TransmissionLine(75.0, 1.0, 3.0E8)
TL2 = rf.TransmissionLine(35.0, 1.5, 3.0E8)
SPL1 = rf.Splice(3)
SC1 = rf.ShortCircuit()

# Declare the overall S-matrix
SMatrix = np.zeros([Ns, 2, 2], dtype=np.complex128)

# Set timer
T = rf.Timer()
T.tic()

# Loop over the iterations
for i in range(0, Ns):

    # Define devices: Index, Number of ports
    d1 = rf.Device(1, 3)
    d2 = rf.Device(2, 2)
    d3 = rf.Device(3, 2)
    d4 = rf.Device(4, 1)

    # Set devies S-matrices
    d1.setSMatrix(SPL1.getSMatrix())
    d2.setSMatrix(TL1.getSMatrix(Z0, freq[i]))
    d3.setSMatrix(TL2.getSMatrix(Z0, freq[i]))
    d4.setSMatrix(SC1.getSMatrix())

    # Create a network of devices
    n1 = rf.Network(1, [d1, d2, d3, d4])

    # Add connections between ports
    n1.connect(d2, 2, d1, 1)
    n1.connect(d1, 2, d3, 1)
    n1.connect(d3, 2, d4, 1)

    # Define excitation ports: Device, Local Port, Global Port
    n1.excitation(d2, 1, 1)
    n1.excitation(d1, 3, 2)

    # Solve the system
    n1.solve()

    # Update the results S-matrix
    SMatrix[i, :, :] = n1.SMatrix

    continue

# Unser timer
T.toc()

# Save results into data files: final S-matrix, filename
rf.saveSMatrix(SMatrix, freq, "Data/ResultsData.csv")

