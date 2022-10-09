//
#include "../include/TestBench.h"

void testNetwork1(){
	
	Matrix SM1={.data=NULL};
	Matrix SM2={.data=NULL};
	
	allocateMatrix(&SM1, 2, 2);
	allocateMatrix(&SM2, 2, 2);
	
	SM1.data[0][0] = 0.111;
	SM1.data[0][1] = 0.112;
	SM1.data[1][0] = 0.121;
	SM1.data[1][1] = 0.122;
	
	SM2.data[0][0] = 0.211;
	SM2.data[0][1] = 0.212;
	SM2.data[1][0] = 0.221;
	SM2.data[1][1] = 0.222;
	
	Network N1=createNetwork(1, 2);
	Device D1=createDevice(1, 1, 2, SM1);
	Device D2=createDevice(1, 2, 2, SM2);
	
	addDevice(&N1, D1);
	addDevice(&N1, D2);
	
	connectPorts(&N1, 2, 1, 1, 2);
	
	addExcitation(&N1, 1, 1, 1);
	addExcitation(&N1, 2, 2, 2);
	
	logNetwork(N1);
	
	solveNetwork(&N1);
	printf("\n");
	printMatrix(N1.S_Matrix);
	
	deallocateMatrix(&SM1);
	deallocateMatrix(&SM2);
	
	deleteNetwork(N1);
}

void testNetwork2(){
	
	Matrix SM1={.data=NULL};
	Matrix SM2={.data=NULL};
	Matrix SM3={.data=NULL};
	Matrix SM4={.data=NULL};
	
	allocateMatrix(&SM1, 2, 2);
	allocateMatrix(&SM2, 2, 2);
	allocateMatrix(&SM3, 2, 2);
	allocateMatrix(&SM4, 3, 3);
	
	SM1.data[0][0] = 0.111;
	SM1.data[0][1] = 0.112;
	SM1.data[1][0] = 0.121;
	SM1.data[1][1] = 0.122;
	
	SM2.data[0][0] = 0.211;
	SM2.data[0][1] = 0.212;
	SM2.data[1][0] = 0.221;
	SM2.data[1][1] = 0.222;
	
	SM3.data[0][0] = 0.311;
	SM3.data[0][1] = 0.312;
	SM3.data[1][0] = 0.321;
	SM3.data[1][1] = 0.322;
	
	SM4.data[0][0] = -1.0/3.0;
	SM4.data[0][1] = +2.0/3.0;
	SM4.data[0][2] = +2.0/3.0;
	SM4.data[1][0] = +2.0/3.0;
	SM4.data[1][1] = -1.0/3.0;
	SM4.data[1][2] = +2.0/3.0;
	SM4.data[2][0] = +2.0/3.0;
	SM4.data[2][1] = +2.0/3.0;
	SM4.data[2][2] = -1.0/3.0;
	
	Network N1=createNetwork(1, 4);
	Device D1=createDevice(1, 1, 2, SM1);
	Device D2=createDevice(1, 2, 2, SM2);
	Device D3=createDevice(1, 3, 2, SM3);
	Device D4=createDevice(1, 4, 3, SM4);
	
	addDevice(&N1, D1);
	addDevice(&N1, D2);
	addDevice(&N1, D3);
	addDevice(&N1, D4);
	
	connectPorts(&N1, 2, 1, 1, 4);
	connectPorts(&N1, 1, 2, 2, 4);
	connectPorts(&N1, 1, 3, 3, 4);
	
	addExcitation(&N1, 1, 1, 1);
	addExcitation(&N1, 2, 2, 2);
	addExcitation(&N1, 3, 2, 3);
	
	logNetwork(N1);
	
	solveNetwork(&N1);
	printMatrix(N1.S_Matrix);
	
	deallocateMatrix(&SM1);
	deallocateMatrix(&SM2);
	deallocateMatrix(&SM3);
	deallocateMatrix(&SM4);
	
	deleteNetwork(N1);
}

void testNetwork3(){
	
	Matrix SM1={.data=NULL};
	Matrix SM2={.data=NULL};
	
	allocateMatrix(&SM1, 2, 2);
	allocateMatrix(&SM2, 1, 1);
	
	SM1.data[0][0] = 0.11;
	SM1.data[0][1] = 0.12;
	SM1.data[1][0] = 0.21;
	SM1.data[1][1] = 0.22;
	
	SM2.data[0][0] = -1.0;
	
	Network N1=createNetwork(1, 2);
	Device D1=createDevice(1, 1, 2, SM1);
	Device D2=createDevice(1, 2, 1, SM2);
	
	addDevice(&N1, D1);
	addDevice(&N1, D2);
	
	connectPorts(&N1, 2, 1, 1, 2);
	
	addExcitation(&N1, 1, 1, 1);
	
	logNetwork(N1);
	
	solveNetwork(&N1);
	printf("\n");
	printMatrix(N1.S_Matrix);
	
	deallocateMatrix(&SM1);
	deallocateMatrix(&SM2);
	
	deleteNetwork(N1);
}

void testNetwork4(){
	
	double GHz=1.0E9;
	int Ns=1601;
	double freq_min=1.0*GHz;
	double freq_max=2.0*GHz;
	double dfreq=(freq_max-freq_min)/(Ns-1.0);
	double freq;
	
	FILE *file=fopen("Data/Test.dat", "w");
	assert(file!=NULL);
	
	Matrix SM1={.data=NULL};
	Matrix SM2={.data=NULL};
	Matrix SM3={.data=NULL};
	Matrix SM4={.data=NULL};
	allocateMatrix(&SM1, 2, 2);
	allocateMatrix(&SM2, 3, 3);
	allocateMatrix(&SM3, 2, 2);
	allocateMatrix(&SM4, 1, 1);
	
	for (int i=0; i<Ns; i++){
		freq = freq_min+i*dfreq;
		
		Network network1 = createNetwork(1, 4);
		
		getSMatrixTL(50.0, 75.0, 1.0, 3.0E8, 0.0, freq, &SM1);
		getSMatrixSplice(3, &SM2);
		getSMatrixTL(50.0, 35.0, 1.5, 3.0E8, 0.0, freq, &SM3);
		getSMatrixShortCircuit(&SM4);
		
		Device device1 = createDevice(1, 1, 2, SM1);
		Device device2 = createDevice(1, 2, 3, SM2);
		Device device3 = createDevice(1, 3, 2, SM3);
		Device device4 = createDevice(1, 4, 1, SM4);
		
		addDevice(&network1, device1);
		addDevice(&network1, device2);
		addDevice(&network1, device3);
		addDevice(&network1, device4);
		
		connectPorts(&network1, 2, 1, 2, 2);
		connectPorts(&network1, 1, 3, 1, 2);
		connectPorts(&network1, 2, 3, 1, 4);
		
		addExcitation(&network1, 1, 1, 1);
		addExcitation(&network1, 2, 3, 2);
		
		// logNetwork(network1);
		
		solveNetwork(&network1);
		
		int n_ports=network1.n_excitations;
		fprintf(file, "%21.14E ", freq);
		for (int ii=0; ii<n_ports; ii++){
			for (int jj=0; jj<n_ports; jj++){
				complex double Sij=network1.S_Matrix.data[ii][jj];
				fprintf(file, "%21.14E %21.14E ", creal(Sij), cimag(Sij));
			}
		}
		fprintf(file, "\n");
		deleteNetwork(network1);
		
	}
	fclose(file);
	
	deallocateMatrix(&SM1);
	deallocateMatrix(&SM2);
	deallocateMatrix(&SM3);
	deallocateMatrix(&SM4);
	
}

void testNetwork5(){
	
	int Ns=1601;
	
	FILE *file=fopen("Data/Test/data.dat", "w");
	assert(file!=NULL);
	
	FILE *file1=fopen("Data/Devices/Device1.dat", "r");
	assert(file1!=NULL);
	FILE *file2=fopen("Data/Devices/Device2.dat", "r");
	assert(file2!=NULL);
	FILE *file3=fopen("Data/Devices/Device3.dat", "r");
	assert(file3!=NULL);
	
	double freq;
	double tempR, tempI;
	complex double j=csqrt(-1.0);
	
	Matrix SM1={.data=NULL};
	Matrix SM2={.data=NULL};
	Matrix SM3={.data=NULL};
	Matrix SM4={.data=NULL};
	Matrix SM5={.data=NULL};
	Matrix SM6={.data=NULL};
	Matrix SM7={.data=NULL};
	Matrix SM8={.data=NULL};
	Matrix SM9={.data=NULL};
	Matrix SM10={.data=NULL};
	Matrix SM11={.data=NULL};
	Matrix SM12={.data=NULL};
	Matrix SM13={.data=NULL};
	allocateMatrix(&SM1, 8, 8);
	allocateMatrix(&SM2, 4, 4);
	allocateMatrix(&SM3, 4, 4);
	allocateMatrix(&SM4, 1, 1);
	allocateMatrix(&SM5, 1, 1);
	allocateMatrix(&SM6, 1, 1);
	allocateMatrix(&SM7, 2, 2);
	allocateMatrix(&SM8, 1, 1);
	allocateMatrix(&SM9, 2, 2);
	allocateMatrix(&SM10, 1, 1);
	allocateMatrix(&SM11, 1, 1);
	allocateMatrix(&SM12, 1, 1);
	allocateMatrix(&SM13, 1, 1);
	
	double Z0=50.0;
	
	Timer T;
	ticTimer(&T);
	
	for (int i=0; i<Ns; i++){
		fscanf(file1, "%lf", &freq);
		for (int ii=0; ii<8; ii++){
			for (int jj=0; jj<8; jj++){
				fscanf(file1, "%lf", &tempR);
				fscanf(file1, "%lf", &tempI);
				SM1.data[ii][jj] = tempR+j*tempI;
			}
		}
		
		fscanf(file2, "%lf", &freq);
		for (int ii=0; ii<4; ii++){
			for (int jj=0; jj<4; jj++){
				fscanf(file2, "%lf", &tempR);
				fscanf(file2, "%lf", &tempI);
				SM2.data[ii][jj] = tempR+j*tempI;
			}
		}
		
		fscanf(file3, "%lf", &freq);
		for (int ii=0; ii<4; ii++){
			for (int jj=0; jj<4; jj++){
				fscanf(file3, "%lf", &tempR);
				fscanf(file3, "%lf", &tempI);
				SM3.data[ii][jj] = tempR+j*tempI;
			}
		}
		
		Network network1 = createNetwork(1, 13);
		
		getSMatrixTerminationImp(Z0, 1000.0, &SM4);
		getSMatrixTerminationImp(Z0, 500.0, &SM5);
		getSMatrixTerminationImp(Z0, 10.0, &SM6);
		getSMatrixSeriesImp(Z0, 100.0, &SM7);
		getSMatrixShortCircuit(&SM8);
		getSMatrixShuntImp(Z0, 1.0/50, &SM9);
		getSMatrixOpenCircuit(&SM10);
		getSMatrixTerminationImp(Z0, 5000.0, &SM11);
		getSMatrixTerminationImp(Z0, 50.0, &SM12);
		getSMatrixTerminationImp(Z0, 10.0, &SM13);
		
		Device device1 = createDevice(1, 1, 8, SM1);
		Device device2 = createDevice(1, 2, 4, SM2);
		Device device3 = createDevice(1, 3, 4, SM3);
		Device device4 = createDevice(1, 4, 1, SM4);
		Device device5 = createDevice(1, 5, 1, SM5);
		Device device6 = createDevice(1, 6, 1, SM6);
		Device device7 = createDevice(1, 7, 2, SM7);
		Device device8 = createDevice(1, 8, 1, SM8);
		Device device9 = createDevice(1, 9, 2, SM9);
		Device device10 = createDevice(1, 10, 1, SM10);
		Device device11 = createDevice(1, 11, 1, SM11);
		Device device12 = createDevice(1, 12, 1, SM12);
		Device device13 = createDevice(1, 13, 1, SM13);
		
		addDevice(&network1, device1);
		addDevice(&network1, device2);
		addDevice(&network1, device3);
		addDevice(&network1, device4);
		addDevice(&network1, device5);
		addDevice(&network1, device6);
		addDevice(&network1, device7);
		addDevice(&network1, device8);
		addDevice(&network1, device9);
		addDevice(&network1, device10);
		addDevice(&network1, device11);
		addDevice(&network1, device12);
		addDevice(&network1, device13);
		
		connectPorts(&network1, 2, 1, 1, 4);
		connectPorts(&network1, 3, 1, 1, 5);
		connectPorts(&network1, 4, 1, 1, 6);
		connectPorts(&network1, 5, 1, 1, 7);
		connectPorts(&network1, 6, 1, 1, 8);
		connectPorts(&network1, 7, 1, 1, 9);
		connectPorts(&network1, 8, 1, 2, 3);
		connectPorts(&network1, 2, 7, 1, 2);
		connectPorts(&network1, 1, 10, 2, 2);
		connectPorts(&network1, 2, 9, 1, 3);
		connectPorts(&network1, 3, 2, 1, 11);
		connectPorts(&network1, 4, 2, 1, 12);
		connectPorts(&network1, 3, 3, 1, 13);
		
		addExcitation(&network1, 1, 1, 1);
		addExcitation(&network1, 3, 4, 2);
		
		// logNetwork(network1);
		
		solveNetwork(&network1);
		
		int n_ports=network1.n_excitations;
		fprintf(file, "%21.14E ", freq);
		for (int ii=0; ii<n_ports; ii++){
			for (int jj=0; jj<n_ports; jj++){
				complex double Sij=network1.S_Matrix.data[ii][jj];
				fprintf(file, "%21.14E %21.14E ", creal(Sij), cimag(Sij));
			}
		}
		fprintf(file, "\n");
		deleteNetwork(network1);
		
	}
	
	tocTimer(&T);
	
	fclose(file1);
	fclose(file2);
	fclose(file3);
	fclose(file);
	
	deallocateMatrix(&SM1);
	deallocateMatrix(&SM2);
	deallocateMatrix(&SM3);
	deallocateMatrix(&SM4);
	deallocateMatrix(&SM5);
	deallocateMatrix(&SM6);
	deallocateMatrix(&SM7);
	deallocateMatrix(&SM8);
	deallocateMatrix(&SM9);
	deallocateMatrix(&SM10);
	deallocateMatrix(&SM11);
	deallocateMatrix(&SM12);
	deallocateMatrix(&SM13);
	
}