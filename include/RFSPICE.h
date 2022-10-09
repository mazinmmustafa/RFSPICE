#ifndef RFSPICE_H
#define RFSPICE_H

// Definitions
#include "myLib.h"
#include "Utilities.h"
#include "Matrix.h"

// Connect 
typedef struct Connect Connect;
struct Connect{
	char type;
	int global_index;
	complex double weight;
};

// Node
typedef struct Node Node;
struct Node{
	char type; 
	int port_local_index, port_global_index;
	int device_index, network_index;
	int n_connects;
	Connect *inputs;
};

// Port
typedef struct Port Port;
struct Port{
	int local_index, global_index, global_excitation_port_index;
	int device_index, network_index;
	int isExcitation, isConnected;
	Node *nodes;
};

// Device
typedef struct Device Device;
struct Device{
	int n_ports;
	int device_index, network_index;
	Matrix S_Matrix;
	Port *ports;
	int isChecked;
};

// Network
typedef struct Network Network;
struct Network{
	int network_index;
	int n_devices;
	Device *devices;
	Matrix S_Matrix;
	int *isDeviceSet;
	int isChecked;
	int n_excitations;
};

// Functions
Network createNetwork(int network_index, int n_devices);
void logNetwork(Network network);
Device createDevice(int network_index, int device_index, 
	int n_ports, Matrix S_Matrix);
void addDevice(Network *network, Device device);
void connectPorts(Network *network, int port_s_index, 
	int device_s_index, int port_d_index, 
	int device_d_index);
void addExcitation(Network *network, int device_index, 
	int port_index, int port_global_index);
void solveNetwork(Network *network);
void deleteNetwork(Network network);
//
void getSMatrixTL(double Z0, complex double Zc, 
	double L, double vp, double alpha, double freq, 
	Matrix *S_Matrix);
void getSMatrixSeriesImp(double Z0, complex double Z, 
	Matrix *S_Matrix);
void getSMatrixShuntImp(double Z0, complex double Y, 
	Matrix *S_Matrix);
void getSMatrixTerminationImp(double Z0, 
	complex double Z, Matrix *S_Matrix);
void getSMatrixOpenCircuit(Matrix *S_Matrix);
void getSMatrixShortCircuit(Matrix *S_Matrix);
void getSMatrixSplice(int N, Matrix *S_Matrix);	

#endif