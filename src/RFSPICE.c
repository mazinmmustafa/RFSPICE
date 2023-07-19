//
#include "../include/RFSPICE.h"

// Ideal Devices Library
void getSMatrixTL(double Z0, complex double Zc, double L, double vp, 
	double alpha, double freq, Matrix *S_Matrix){
	assert((*S_Matrix).rows=2);
	assert((*S_Matrix).cols=2);
	double omega=2.0*pi*freq;
	double beta=omega/vp;
	complex double j=csqrt(-1.0);
	complex double k=beta-j*alpha;
	complex double Gamma_s=(Zc-Z0)/(Zc+Z0);
	complex double Gamma_l=-Gamma_s;
	complex double Gamma_0=Gamma_l*cexp(-j*2.0*k*L);
	complex double S11=(Gamma_s+Gamma_0)/(1.0+Gamma_s*Gamma_0);
	complex double S21=(1.0+Gamma_l)/(1.0+Gamma_s*Gamma_0);
	S21*=(2.0*Zc/(Z0+Zc))*cexp(-j*k*L);
	(*S_Matrix).data[0][0] = S11;
	(*S_Matrix).data[0][1] = S21;
	(*S_Matrix).data[1][0] = S21;
	(*S_Matrix).data[1][1] = S11;
}

void getSMatrixSeriesImp(double Z0, complex double Z, Matrix *S_Matrix){
	assert((*S_Matrix).rows=2);
	assert((*S_Matrix).cols=2);
	complex double S11=Z/(Z+2.0*Z0);
	complex double S21=2.0*Z0/(Z+2.0*Z0);
	(*S_Matrix).data[0][0] = S11;
	(*S_Matrix).data[0][1] = S21;
	(*S_Matrix).data[1][0] = S21;
	(*S_Matrix).data[1][1] = S11;
}

void getSMatrixShuntImp(double Z0, complex double Y, Matrix *S_Matrix){
	assert((*S_Matrix).rows=2);
	assert((*S_Matrix).cols=2);
	complex double S11=-Z0*Y/(2.0+Z0*Y);
	complex double S21=2.0/(2.0+Z0*Y);
	(*S_Matrix).data[0][0] = S11;
	(*S_Matrix).data[0][1] = S21;
	(*S_Matrix).data[1][0] = S21;
	(*S_Matrix).data[1][1] = S11;
}

void getSMatrixTerminationImp(double Z0, complex double Z, Matrix *S_Matrix){
	assert((*S_Matrix).rows=1);
	assert((*S_Matrix).cols=1);
	complex double S11=(Z-Z0)/(Z+Z0);
	(*S_Matrix).data[0][0] = S11;
}

void getSMatrixOpenCircuit(Matrix *S_Matrix){
	assert((*S_Matrix).rows=1);
	assert((*S_Matrix).cols=1);
	(*S_Matrix).data[0][0] = +1.0;
}

void getSMatrixShortCircuit(Matrix *S_Matrix){
	assert((*S_Matrix).rows=1);
	assert((*S_Matrix).cols=1);
	(*S_Matrix).data[0][0] = -1.0;
}

void getSMatrixSplice(int N, Matrix *S_Matrix){
	assert(N>2);
	assert((*S_Matrix).rows=N);
	assert((*S_Matrix).cols=N);
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){
			if (i==j){
				(*S_Matrix).data[i][j] = (2.0/N)-1.0;
			}else{
				(*S_Matrix).data[i][j] = (2.0/N);
			}
		}
	}
}

//

int global_index_counter=1;
void logNetwork(Network network){
	FILE *file=fopen("log.txt", "w");
	assert(file!=NULL);
	fprintf(file, "=================================================\n");
	fprintf(file, "Network %d:\n\n", network.network_index);
	for (int i=0; i<network.n_devices; i++){
		fprintf(file, "------------------------\n");
		fprintf(file, "Device %d:\n", network.devices[i].device_index);
		fprintf(file, "\n");
		for (int j=0; j<network.devices[i].n_ports; j++){
			fprintf(file, "Port %d: ", network.devices[i].ports[j].local_index);
			fprintf(file, "Global Index %d", network.devices[i].ports[j].global_index);
			if (network.devices[i].ports[j].isExcitation==1){
				fprintf(file, " (Excitation Port %d)\n", network.devices[i].ports[j].global_excitation_port_index);
			}else{
				fprintf(file, "\n");
			}
			// Print Node a
			fprintf(file, "Node %c%d: ", network.devices[i].ports[j].nodes[0].type, 
				network.devices[i].ports[j].nodes[0].port_local_index);
			fprintf(file, "Global Index %c%d, inputs from %d nodes\n", network.devices[i].ports[j].nodes[0].type,  
				network.devices[i].ports[j].nodes[0].port_global_index, network.devices[i].ports[j].nodes[0].n_connects);
			for (int k=0; k<network.devices[i].ports[j].nodes[0].n_connects; k++){
				fprintf(file, "Connected to %c%d with weight (%21.14E, %21.14E)\n", 
					network.devices[i].ports[j].nodes[0].inputs[k].type,
					network.devices[i].ports[j].nodes[0].inputs[k].global_index,
					creal(network.devices[i].ports[j].nodes[0].inputs[k].weight),
					cimag(network.devices[i].ports[j].nodes[0].inputs[k].weight));
			}
			// Print Node b
			fprintf(file, "Node %c%d: ", network.devices[i].ports[j].nodes[1].type, 
				network.devices[i].ports[j].nodes[1].port_local_index);
			fprintf(file, "Global Index %c%d, inputs from %d nodes\n", network.devices[i].ports[j].nodes[1].type,  
				network.devices[i].ports[j].nodes[1].port_global_index, network.devices[i].ports[j].nodes[1].n_connects);
			for (int k=0; k<network.devices[i].ports[j].nodes[1].n_connects; k++){
				fprintf(file, "Connected to %c%d with weight (%21.14E, %21.14E)\n", 
					network.devices[i].ports[j].nodes[1].inputs[k].type,
					network.devices[i].ports[j].nodes[1].inputs[k].global_index,
					creal(network.devices[i].ports[j].nodes[1].inputs[k].weight),
					cimag(network.devices[i].ports[j].nodes[1].inputs[k].weight));
			}
		}
		fprintf(file, "\n");
	}
	fprintf(file, "=================================================\n");
	fclose(file);
}

int global_network_index_counter=0;
Network createNetwork(int network_index, int n_devices){
	assert(n_devices>0);
	Network newNetwork;
	newNetwork.network_index = network_index;
	newNetwork.n_devices = n_devices;
	newNetwork.devices = (Device*)malloc(n_devices*sizeof(Device));
	newNetwork.isDeviceSet = (int*)malloc(n_devices*sizeof(int));
	for (int i=0; i<n_devices; i++){
		newNetwork.isDeviceSet[i] = 0;
	}
	newNetwork.isChecked = 0;
	global_network_index_counter++;
	return newNetwork;
}

void deleteNetwork(Network *network){
	for (int i=0; i<network->n_devices; i++){
		for (int j=0; j<network->devices[i].n_ports; j++){
			free(network->devices[i].ports[j].nodes[0].inputs);
			free(network->devices[i].ports[j].nodes[1].inputs);
		}
		for (int j=0; j<network->devices[i].n_ports; j++){
			free(network->devices[i].ports[j].nodes);
			network->devices[i].ports[j].nodes = NULL;
		}
		free(network->devices[i].ports);
		network->devices[i].ports = NULL;
	}
	free(network->devices);
	network->devices = NULL;
	free(network->isDeviceSet);
	deallocateMatrix(&network->S_Matrix);
}

Device createDevice(int network_index, int device_index, 
	int n_ports, Matrix S_Matrix){
	assert(device_index>=0);
	if (network_index<1||network_index>global_network_index_counter){
		printf("ERROR: Wrong Network Index!\n");
		exit(1);
	}
	assert(S_Matrix.rows==S_Matrix.cols);
	assert(n_ports==S_Matrix.rows);
	Device newDevice;
	newDevice.n_ports = n_ports;
	newDevice.network_index = network_index;
	newDevice.device_index = device_index;
	newDevice.S_Matrix = S_Matrix;
	newDevice.ports = (Port*)malloc(n_ports*sizeof(Port));
	for (int i=0; i<n_ports; i++){
		newDevice.ports[i].network_index = network_index;
		newDevice.ports[i].device_index = device_index;
		newDevice.ports[i].local_index = i+1;
		newDevice.ports[i].global_index = global_index_counter;
		newDevice.ports[i].isExcitation = 0;
		newDevice.ports[i].isConnected = 0;
		newDevice.ports[i].global_excitation_port_index = 0;
		newDevice.ports[i].nodes = (Node*)malloc(2*sizeof(Node));
		newDevice.ports[i].nodes[0].type = 'a';
		newDevice.ports[i].nodes[0].network_index = network_index;
		newDevice.ports[i].nodes[0].device_index = device_index;
		newDevice.ports[i].nodes[0].port_local_index = i+1;
		newDevice.ports[i].nodes[0].n_connects = 0;
		newDevice.ports[i].nodes[0].port_global_index = global_index_counter;
		newDevice.ports[i].nodes[1].type = 'b';
		newDevice.ports[i].nodes[1].network_index = network_index;
		newDevice.ports[i].nodes[1].device_index = device_index;
		newDevice.ports[i].nodes[1].port_local_index = i+1;
		newDevice.ports[i].nodes[1].n_connects = 0;
		newDevice.ports[i].nodes[1].port_global_index = global_index_counter;
		global_index_counter++;
	}
	newDevice.isChecked = 0;
	return newDevice;
}

void addDevice(Network *network, Device device){
	assert(network->network_index==device.network_index);
	assert(device.device_index<=(network->n_devices));
	if (network->isDeviceSet[device.device_index-1]==1){
		printf("ERROR: Device Is Already Set!\n");
		exit(1);
	}
	network->devices[device.device_index-1]=device;
	network->isDeviceSet[device.device_index-1]=1;
}

void connectPorts(Network *network, int port_s_index, 
	int device_s_index, int port_d_index, int device_d_index){
	if (network->devices[device_s_index-1].ports[port_s_index-1].isConnected!=0){
		printf("ERROR: Port Is Already Connected!\n");
		exit(1);
	}
	if (network->devices[device_d_index-1].ports[port_d_index-1].isConnected!=0){
		printf("ERROR: Port Is Already Connected!\n");
		exit(1);
	}
	if (device_s_index<1||device_s_index>network->n_devices){
		printf("ERROR: Wrong Device Index!\n");
		exit(1);
	}
	if (device_d_index<1||device_d_index>network->n_devices){
		printf("ERROR: Wrong Device Index!\n");
		exit(1);
	}
	assert(network->isDeviceSet[device_s_index-1]==1);
	assert(network->isDeviceSet[device_d_index-1]==1);
	assert(network->devices[device_s_index-1].ports[port_s_index-1].isConnected==0);
	assert(network->devices[device_d_index-1].ports[port_d_index-1].isConnected==0);
	int device_s_n_ports = network->devices[device_s_index-1].n_ports;
	int device_d_n_ports = network->devices[device_d_index-1].n_ports;
	network->devices[device_s_index-1].ports[port_s_index-1].nodes[0].inputs=(Connect*)malloc((1)*sizeof(Connect));
	network->devices[device_s_index-1].ports[port_s_index-1].nodes[1].inputs=(Connect*)malloc((device_s_n_ports)*sizeof(Connect));
	network->devices[device_d_index-1].ports[port_d_index-1].nodes[0].inputs=(Connect*)malloc((1)*sizeof(Connect));
	network->devices[device_d_index-1].ports[port_d_index-1].nodes[1].inputs=(Connect*)malloc((device_d_n_ports)*sizeof(Connect));
	network->devices[device_s_index-1].ports[port_s_index-1].nodes[0].n_connects = 1;
	network->devices[device_s_index-1].ports[port_s_index-1].nodes[1].n_connects = device_s_n_ports;
	network->devices[device_d_index-1].ports[port_d_index-1].nodes[0].n_connects = 1;
	network->devices[device_d_index-1].ports[port_d_index-1].nodes[1].n_connects = device_d_n_ports;
	Connect newConnect_a_d={'b', network->devices[device_s_index-1].ports[port_s_index-1].global_index, 1.0};
	Connect newConnect_a_s={'b', network->devices[device_d_index-1].ports[port_d_index-1].global_index, 1.0};
	network->devices[device_d_index-1].ports[port_d_index-1].nodes[0].inputs[0]=newConnect_a_d;
	network->devices[device_s_index-1].ports[port_s_index-1].nodes[0].inputs[0]=newConnect_a_s;
	for (int i=0; i<network->devices[device_s_index-1].n_ports; i++){
		Connect newConnect_b_s={'a', network->devices[device_s_index-1].ports[i].global_index, 
			network->devices[device_s_index-1].S_Matrix.data[network->devices[device_s_index-1].ports[port_s_index-1].local_index-1][i]};
		network->devices[device_s_index-1].ports[port_s_index-1].nodes[1].inputs[i]=newConnect_b_s;
	}
	for (int i=0; i<network->devices[device_d_index-1].n_ports; i++){
		Connect newConnect_b_d={'a', network->devices[device_d_index-1].ports[i].global_index, 
			network->devices[device_d_index-1].S_Matrix.data[network->devices[device_d_index-1].ports[port_d_index-1].local_index-1][i]};
		network->devices[device_d_index-1].ports[port_d_index-1].nodes[1].inputs[i]=newConnect_b_d;
	}
	network->devices[device_s_index-1].ports[port_s_index-1].isConnected=1;
	network->devices[device_d_index-1].ports[port_d_index-1].isConnected=1;
}

int global_excitation_port_index_counter=1;
void addExcitation(Network *network, int device_index, 
	int port_index, int port_global_index){
	assert(network->isDeviceSet[device_index-1]==1);
	assert(network->devices[device_index-1].ports[port_index-1].isConnected==0);
	assert(network->devices[device_index-1].ports[port_index-1].isExcitation==0);
	assert(port_global_index>0);
	if (port_global_index>global_excitation_port_index_counter){
		printf("ERROR: Wrong Excitation Port Index!\n");
		exit(1);
	}
	network->devices[device_index-1].ports[port_index-1].global_excitation_port_index=port_global_index;
	network->devices[device_index-1].ports[port_index-1].isExcitation=1;
	int device_n_ports = network->devices[device_index-1].n_ports;
	network->devices[device_index-1].ports[port_index-1].nodes[0].inputs=NULL;
	network->devices[device_index-1].ports[port_index-1].nodes[1].inputs=(Connect*)malloc((device_n_ports)*sizeof(Connect));
	network->devices[device_index-1].ports[port_index-1].nodes[1].n_connects = device_n_ports;
	for (int i=0; i<network->devices[device_index-1].n_ports; i++){
		Connect newConnect={'a', network->devices[device_index-1].ports[i].global_index, 
			network->devices[device_index-1].S_Matrix.data[network->devices[device_index-1].ports[port_index-1].local_index-1][i]};
		network->devices[device_index-1].ports[port_index-1].nodes[1].inputs[i]=newConnect;
	}
	network->n_excitations++;
	global_excitation_port_index_counter++;
}

typedef struct DataList DataList;
struct DataList{
	int type, device_index, port_index;
};

void solveNetwork(Network *network){
	global_index_counter = 1;
	global_network_index_counter = 0;
	// Checking
	int count=0;
	int count_all=0;
	for (int i=0; i<network->n_devices; i++){
		for (int j=0; j<network->devices[i].n_ports; j++){
			if (network->devices[i].ports[j].isConnected+
				network->devices[i].ports[j].isExcitation!=1){
				printf("ERROR: Missing Excitaions Or Connections!\n");
				exit(1);	
			}
			count_all++;
			if (network->devices[i].ports[j].isConnected){
				count++;	
			}
		}
		network->devices[i].isChecked = 1;
	}
	if (count<1){
		printf("ERROR: Missing Excitaions!\n");
		exit(1);
	}
	network->isChecked = 1;
	//
	int Na=network->n_excitations, Nb=network->n_excitations;
	int Nx=2*(count_all-Na);
	Matrix S_Matrix={.data=NULL};
	allocateMatrix(&S_Matrix, Nb, Na);
	DataList Matrix_a[Na];
	DataList Matrix_b[Nb];
	DataList Matrix_x[Nx];
	int count_a=0, count_b=0, count_x=0;
	for (int i=0; i<network->n_devices; i++){
		for (int j=0; j<network->devices[i].n_ports; j++){
			if (network->devices[i].ports[j].isExcitation){
				Matrix_a[count_a].type = 0;
				Matrix_a[count_a].device_index=network->devices[i].device_index;
				Matrix_a[count_a].port_index=network->devices[i].ports[j].global_index;
				count_a++;
				Matrix_b[count_b].type = 1;
				Matrix_b[count_b].device_index=network->devices[i].device_index;
				Matrix_b[count_b].port_index=network->devices[i].ports[j].global_index;
				count_b++;
			}else{
				Matrix_x[count_x].type = 0;
				Matrix_x[count_x].device_index=network->devices[i].device_index;
				Matrix_x[count_x].port_index=network->devices[i].ports[j].global_index;
				count_x++;
				Matrix_x[count_x].type = 1;
				Matrix_x[count_x].device_index=network->devices[i].device_index;
				Matrix_x[count_x].port_index=network->devices[i].ports[j].global_index;
				count_x++;
			}
		}
	}
	assert(count_a==Na);
	assert(count_b==Nb);
	assert(count_x==Nx);
	Matrix Wba={.data=NULL};
	Matrix Wbx={.data=NULL};
	Matrix Wxa={.data=NULL};
	Matrix Wxx={.data=NULL};
	allocateMatrix(&Wba, Nb, Na);
	allocateMatrix(&Wbx, Nb, Nx);
	allocateMatrix(&Wxa, Nx, Na);
	allocateMatrix(&Wxx, Nx, Nx);
	complex double weight;
	// Wba
	for (int i=0; i<Nb; i++){
		for (int j=0; j<Na; j++){
			weight = 0.0;
			int i_temp=0;
			for (int ii=0; ii<network->devices[Matrix_b[i].device_index-1].n_ports; ii++){
				if (network->devices[Matrix_b[i].device_index-1].ports[ii].global_index==
					Matrix_b[i].port_index){
					i_temp = ii;
				}
			}
			char a='a';
			if (Matrix_a[j].type==0){a='a';};
			if (Matrix_a[j].type==1){a='b';};
			int K=network->devices[Matrix_b[i].device_index-1].ports[i_temp].nodes[Matrix_b[i].type].n_connects;
			for (int k=0; k<K; k++){
				if ((network->devices[Matrix_b[i].device_index-1].ports[i_temp].nodes[Matrix_b[i].type].inputs[k].global_index==Matrix_a[j].port_index)&&
					(network->devices[Matrix_b[i].device_index-1].ports[i_temp].nodes[Matrix_b[i].type].inputs[k].type==a)){
					weight=network->devices[Matrix_b[i].device_index-1].ports[i_temp].nodes[Matrix_b[i].type].inputs[k].weight;
				}
			}
			Wba.data[i][j] = weight;
		}
	}
	// Wbx
	for (int i=0; i<Nb; i++){
		for (int j=0; j<Nx; j++){
			weight = 0.0;
			int i_temp=0;
			for (int ii=0; ii<network->devices[Matrix_b[i].device_index-1].n_ports; ii++){
				if (network->devices[Matrix_b[i].device_index-1].ports[ii].global_index==
					Matrix_b[i].port_index){
					i_temp = ii;
				}
			}
			char x='a';
			if (Matrix_x[j].type==0){x='a';};
			if (Matrix_x[j].type==1){x='b';};
			int K=network->devices[Matrix_b[i].device_index-1].ports[i_temp].nodes[Matrix_b[i].type].n_connects;
			for (int k=0; k<K; k++){
				if ((network->devices[Matrix_b[i].device_index-1].ports[i_temp].nodes[Matrix_b[i].type].inputs[k].global_index==Matrix_x[j].port_index)&&
					(network->devices[Matrix_b[i].device_index-1].ports[i_temp].nodes[Matrix_b[i].type].inputs[k].type==x)){
					weight=network->devices[Matrix_b[i].device_index-1].ports[i_temp].nodes[Matrix_b[i].type].inputs[k].weight;
				}
			}
			Wbx.data[i][j] = weight;
		}
	}
	// Wxa
	for (int i=0; i<Nx; i++){
		for (int j=0; j<Na; j++){
			weight = 0.0;
			int i_temp=0;
			for (int ii=0; ii<network->devices[Matrix_x[i].device_index-1].n_ports; ii++){
				if (network->devices[Matrix_x[i].device_index-1].ports[ii].global_index==
					Matrix_x[i].port_index){
					i_temp = ii;
				}
			}
			char a='a';
			if (Matrix_a[j].type==0){a='a';};
			if (Matrix_a[j].type==1){a='b';};
			int K=network->devices[Matrix_x[i].device_index-1].ports[i_temp].nodes[Matrix_x[i].type].n_connects;
			for (int k=0; k<K; k++){
				if ((network->devices[Matrix_x[i].device_index-1].ports[i_temp].nodes[Matrix_x[i].type].inputs[k].global_index==Matrix_a[j].port_index)&&
					(network->devices[Matrix_x[i].device_index-1].ports[i_temp].nodes[Matrix_x[i].type].inputs[k].type==a)){
					weight=network->devices[Matrix_x[i].device_index-1].ports[i_temp].nodes[Matrix_x[i].type].inputs[k].weight;
				}
			}
			Wxa.data[i][j] = weight;
		}
	}
	// Wxx
	for (int i=0; i<Nx; i++){
		for (int j=0; j<Nx; j++){
			weight = 0.0;
			int i_temp=0;
			for (int ii=0; ii<network->devices[Matrix_x[i].device_index-1].n_ports; ii++){
				if (network->devices[Matrix_x[i].device_index-1].ports[ii].global_index==
					Matrix_x[i].port_index){
					i_temp = ii;
				}
			}
			char x='a';
			if (Matrix_x[j].type==0){x='a';};
			if (Matrix_x[j].type==1){x='b';};
			int K=network->devices[Matrix_x[i].device_index-1].ports[i_temp].nodes[Matrix_x[i].type].n_connects;
			for (int k=0; k<K; k++){
				if ((network->devices[Matrix_x[i].device_index-1].ports[i_temp].nodes[Matrix_x[i].type].inputs[k].global_index==Matrix_x[j].port_index)&&
					(network->devices[Matrix_x[i].device_index-1].ports[i_temp].nodes[Matrix_x[i].type].inputs[k].type==x)){
					weight=network->devices[Matrix_x[i].device_index-1].ports[i_temp].nodes[Matrix_x[i].type].inputs[k].weight;
				}
			}
			Wxx.data[i][j] = weight;
		}
	}
	//
	for (int i=0; i<Nx; i++){
		for (int j=0; j<Nx; j++){
			if (i==j){
				Wxx.data[i][j] = 1.0-Wxx.data[i][j];
			}else{
				Wxx.data[i][j] = -Wxx.data[i][j];
			}
		}
	}
	Matrix temp1={.data=NULL};
	allocateMatrix(&temp1, Nx, Nx);
	inverseMatrix(Wxx, &temp1);
	Matrix temp2={.data=NULL};
	allocateMatrix(&temp2, Nb, Nx);
	multMatrix(Wbx, temp1, &temp2);
	Matrix temp3={.data=NULL};
	allocateMatrix(&temp3, Nb, Na);
	multMatrix(temp2, Wxa, &temp3);
	addMatrix(Wba, temp3, &S_Matrix);
	network->S_Matrix = S_Matrix;
	deallocateMatrix(&temp3);
	deallocateMatrix(&temp2);
	deallocateMatrix(&temp1);
	deallocateMatrix(&Wba);
	deallocateMatrix(&Wbx);
	deallocateMatrix(&Wxa);
	deallocateMatrix(&Wxx);
}
