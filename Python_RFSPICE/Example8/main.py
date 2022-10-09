import rfspice as rf
import numpy as np

# Read S-parameters files for devices
Device1 = rf.readSParemetersFiles("Data/Device1.csv")
Device2 = rf.readSParemetersFiles("Data/Device2.csv")
Device3 = rf.readSParemetersFiles("Data/Device3.csv")

# Obtain the frequnecy vector and number of iterations
freq, Ns = rf.checkDimensions([Device1, Device2, Device3])

# Declare the overall S-matrix
SMatrix = np.zeros([Ns, 2, 2], dtype=np.complex128)

# Define ideal devices
T1 = rf.TerminationImpedance()
T2 = rf.TerminationImpedance()
T3 = rf.TerminationImpedance()
T4 = rf.TerminationImpedance()
T5 = rf.TerminationImpedance()
T6 = rf.TerminationImpedance()
SC1 = rf.ShortCircuit()
OC1 = rf.OpenCircuit()
SR1 = rf.SeriesImpedance()
SH1 = rf.ShuntImpedance()

# Set timer
T = rf.Timer()
T.tic()

Z0 = 50

# Loop over the iterations
for i in range(0, Ns):

    # Get S-matrix for each iteration
    SM1 = rf.getSMatrix(Device1[i, :])
    SM2 = rf.getSMatrix(Device2[i, :])
    SM3 = rf.getSMatrix(Device3[i, :])

    # Define devices: Index, Number of ports
    d1 = rf.Device(1, 8)
    d2 = rf.Device(2, 4)
    d3 = rf.Device(3, 4)

    d4 = rf.Device(4, 1)
    d5 = rf.Device(5, 1)
    d6 = rf.Device(6, 1)
    d7 = rf.Device(7, 1)
    d8 = rf.Device(8, 1)
    d9 = rf.Device(9, 1)
    d10 = rf.Device(10, 1)
    d11 = rf.Device(11, 1)
    d12 = rf.Device(12, 2)
    d13 = rf.Device(13, 2)


    # Set devies S-matrices
    list = []
    d1.setSMatrix(SM1)
    list.append(d1)
    d2.setSMatrix(SM2)
    list.append(d2)
    d3.setSMatrix(SM3)
    list.append(d3)
    d4.setSMatrix(T1.getSMatrix(Z0, 1000.0))
    list.append(d4)
    d5.setSMatrix(T2.getSMatrix(Z0, 500.0))
    list.append(d5)
    d6.setSMatrix(T3.getSMatrix(Z0, 10.0))
    list.append(d6)
    d7.setSMatrix(T4.getSMatrix(Z0, 5000.0))
    list.append(d7)
    d8.setSMatrix(T5.getSMatrix(Z0, 50.0))
    list.append(d8)
    d9.setSMatrix(T6.getSMatrix(Z0, 10.0))
    list.append(d9)
    d10.setSMatrix(SC1.getSMatrix())
    list.append(d10)
    d11.setSMatrix(OC1.getSMatrix())
    list.append(d11)
    d12.setSMatrix(SR1.getSMatrix(Z0, 100.0))
    list.append(d12)
    d13.setSMatrix(SH1.getSMatrix(Z0, 1.0/50.0))
    list.append(d13)

    # Create a network of devices
    n1 = rf.Network(1, list)

    # Add connections between ports
    n1.connect(d1, 2, d4, 1)
    n1.connect(d1, 3, d5, 1)
    n1.connect(d1, 4, d6, 1)
    n1.connect(d1, 5, d12, 1)
    n1.connect(d1, 6, d10, 1)
    n1.connect(d1, 7, d13, 1)
    n1.connect(d1, 8, d3, 1)

    n1.connect(d2, 1, d12, 2)
    n1.connect(d2, 2, d11, 1)
    n1.connect(d2, 3, d7, 1)
    n1.connect(d2, 4, d8, 1)

    n1.connect(d3, 2, d13, 2)
    n1.connect(d3, 4, d9, 1)

    # Define excitation ports: Device, Local Port, Global Port
    n1.excitation(d1, 1, 1)
    n1.excitation(d3, 3, 2)

    # Solve the system
    n1.solve()

    # Update the results S-matrix
    SMatrix[i, :, :] = n1.SMatrix

    continue

# Unser timer
T.toc()

# Save results into data files: final S-matrix, filename
rf.saveSMatrix(SMatrix, freq, "Data/ResultsData.csv")











