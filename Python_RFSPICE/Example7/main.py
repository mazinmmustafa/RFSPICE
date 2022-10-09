import rfspice as rf
import numpy as np

# Read S-parameters files for devices
Device1 = rf.readSParemetersFiles("Data/Device1.csv")
Device2 = rf.readSParemetersFiles("Data/Device2.csv")
Device3 = rf.readSParemetersFiles("Data/Device3.csv")

# Obtain the frequnecy vector and number of iterations
freq, Ns = rf.checkDimensions([Device1, Device2, Device3])

# Declare the overall S-matrix
SMatrix = np.zeros([Ns, 4, 4], dtype=np.complex128)

# Set timer
T = rf.Timer()
T.tic()

# Loop over the iterations
for i in range(0, Ns):

    # Get S-matrix for each iteration
    SM1 = rf.getSMatrix(Device1[i, :])
    SM2 = rf.getSMatrix(Device2[i, :])
    SM3 = rf.getSMatrix(Device3[i, :])

    # Define devices: Index, Number of ports
    d1 = rf.Device(1, 4)
    d2 = rf.Device(2, 2)
    d3 = rf.Device(3, 2)

    # Set devies S-matrices
    d1.setSMatrix(SM1)
    d2.setSMatrix(SM2)
    d3.setSMatrix(SM3)

    # Create a network of devices
    n1 = rf.Network(1, [d1, d2, d3])

    # Add connections between ports
    n1.connect(d1, 3, d2, 1)
    n1.connect(d1, 4, d3, 1)

    # Define excitation ports: Device, Local Port, Global Port
    n1.excitation(d1, 1, 1)
    n1.excitation(d1, 2, 2)
    n1.excitation(d2, 2, 3)
    n1.excitation(d3, 2, 4)

    # Solve the system
    n1.solve()

    # Update the results S-matrix
    SMatrix[i, :, :] = n1.SMatrix

    continue

# Unser timer
T.toc()

# Save results into data files: final S-matrix, filename
rf.saveSMatrix(SMatrix, freq, "Data/ResultsData.csv")

