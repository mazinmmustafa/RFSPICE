import numpy as np
import os

# Constants
eps0 = 8.8541878128E-12
mu0 = 1.25663706212E-6
c0 = 299792458

class Wire():
    def __init__(self, index=0, x=0.0, y=0.0, r=0.0, rI=0.0,\
                 epsr=0.0, tand=0.0, sigma=0.0):
        self.index = index
        self.x = x
        self.y = y
        self.r = r
        self.rI = rI
        self.epsr = epsr
        self.tand = tand
        self.sigma = sigma
        pass
    def __repr__(self):
        str = "Wire {:2d}:\n".format(self.index)+\
              "x\t: {:21.14E}\n".format(self.x)+\
              "y\t: {:21.14E}\n".format(self.y)+\
              "r\t: {:21.14E}\n".format(self.r)+\
              "rI\t: {:21.14E}\n".format(self.rI)+\
              "epsr\t: {:21.14E}\n".format(self.epsr)+\
              "tand\t: {:21.14E}\n".format(self.tand)+\
              "sigma\t: {:21.14E}\n".format(self.sigma)
        return str

class Config():
    def __init__(self, freq=0.0, N=0, wires=[], Ns=50):
        self.freq = freq
        self.N = N
        self.wires = wires
        self.Ns = Ns
        self.domain = 0.0
        pass
    def __repr__(self):
        str = "Configuration:\n"+\
              "Freq\t: {:21.14E}\n".format(self.freq)+\
              "N\t: {:2d}\n".format(self.N)
        for wire in self.wires:
            str+=wire.__repr__()
            continue
        return str
    def check(self):
        try:
            assert(self.N>1)
            assert(self.N==len(self.wires))
            self.wires = sorted(self.wires, key=lambda x: x.index)
            for i in range(0, self.N):
                assert(self.wires[i].index==i)
                continue
            for i in range(0, self.N):
                if self.wires[i].r>=self.wires[i].rI:
                    exit(1)
                continue
        except:
            print("Error: Incorrect Configuration Parameters!")
            exit(1)
        flag = False
        for wire in self.wires:
            S = np.sqrt(wire.x*wire.x+wire.y*wire.y)
            if S>self.domain:
                flag = True
                self.domain = S
                pass
            continue
        if flag:
            self.domain*=5
            pass
        pass

def solveC0Matrix(config, view=True):
    config.check()
    # Write FF Scripe File
    fileName = "FF/FFScript0.edp"
    file = open(fileName, "w")
    file.write("// Output File\n")
    file.write("ofstream file(\"FF/C0Matrix.dat\");\n")
    file.write("\n")
    file.write("// Constants\n")
    file.write("real eps0={:21.14E};\n".format(eps0))
    file.write("int NCond={};\n".format(config.N-1))
    file.write("\n")
    for wire in config.wires:
        file.write("// Conductor {}\n".format(wire.index))
        file.write("real x{}={:21.14E};\n".format(wire.index, wire.x))
        file.write("real y{}={:21.14E};\n".format(wire.index, wire.y))
        file.write("real r{}={:21.14E};\n".format(wire.index, wire.r))
        file.write("\n")
        continue
    file.write("// Domain\n")
    file.write("real rd={:21.14E};\n".format(config.domain))
    file.write("\n")
    file.write("// Mesh Resolution\n")
    file.write("real NsAir={};\n".format(config.Ns))
    file.write("real NsCond={};\n".format(config.Ns))
    file.write("\n")
    file.write("// Mesh\n")
    for wire in config.wires:
        file.write("border C{}(t=0, 2*pi) {{".format(wire.index))
        file.write("x=x{}+r{}*cos(t); y=y{}+r{}*sin(t);}}\n".format(wire.index,\
                   wire.index, wire.index, wire.index, wire.index))
        continue
    file.write("border Cd(t=0, 2*pi){x=rd*cos(t); y=rd*sin(t);}\n")
    file.write("mesh Th=buildmesh(Cd(NsAir)\n")
    for wire in config.wires:
        file.write("\t+C{}(-NsCond)\n".format(wire.index))
        continue
    file.write("\t);\n")
    file.write("\n")
    file.write("// Solver\n")
    if view:
        file.write("plot(Th, wait=true);\n")
        pass
    file.write("fespace Vh(Th,P2);\n")
    file.write("Vh u, v;\n")
    file.write("Vh eps = eps0;\n")
    file.write("real Q;\n")
    file.write("\n")
    file.write("// BCs Vector\n")
    file.write("matrix<real> Vm;\n")
    file.write("Vm.resize(NCond, 1);\n")
    file.write("\n")
    file.write("// For Loop\n")
    file.write("for(int i=0; i<NCond; i++){\n")
    file.write("\tfor(int j=0; j<NCond; j++){\n")
    file.write("\t\tfor(int k=0; k<NCond; k++){if(k==i){Vm(k,0)=+0.5;}else{Vm(k,0)=-0.5;}}\n")
    file.write("\t\t\tsolve V(u, v)=int2d(Th)(eps*(dx(u)*dx(v)+dy(u)*dy(v)))\n")
    file.write("\t\t\t+on(C0, u=-0.5)")
    for i in range(1, config.N):
        file.write("\n\t\t\t+on(C{}, u=Vm({}, 0))".format(i, i-1))
        continue
    file.write(";\n")
    file.write("\t\t\tif (j==0)\n")
    file.write("\t\t\t\t{Q=int1d(Th, C1)(eps*(N.x*dx(u)+N.y*dy(u)));}\n")
    for i in range(2, config.N):
        file.write("\t\t\telse if (j=={})\n".format(i-1))
        file.write("\t\t\t\t{{Q=int1d(Th, C{})(eps*(N.x*dx(u)+N.y*dy(u)));}}\n".format(i))
        continue
    file.write("\t\t\tfile << Q << \" \";\n")
    file.write("\t\t}\n")
    file.write("\t\tfile << endl;\n")
    file.write("}\n")
    file.write("\nexit(0);")
    file.close()
    # Execute FF Scripe File
    os.system("FreeFEM++ -nw -ne -nc {}".format(fileName))
    pass

def solveCMatrix(config, view=True):
    config.check()
    # Write FF Scripe File
    fileName = "FF/FFScriptI.edp"
    file = open(fileName, "w")
    file.write("// Output File\n")
    file.write("ofstream fileR(\"FF/CMatrixR.dat\");\n")
    file.write("ofstream fileI(\"FF/CMatrixI.dat\");\n")
    file.write("\n")
    file.write("// Constants\n")
    file.write("real eps0={:21.14E};\n".format(eps0))
    file.write("int NCond={};\n".format(config.N-1))
    file.write("\n")
    for wire in config.wires:
        file.write("// Conductor {}\n".format(wire.index))
        file.write("real x{}={:21.14E};\n".format(wire.index, wire.x))
        file.write("real y{}={:21.14E};\n".format(wire.index, wire.y))
        file.write("real r{}={:21.14E};\n".format(wire.index, wire.r))
        file.write("real rI{}={:21.14E};\n".format(wire.index, wire.rI))
        file.write("real EPS{}={:21.14E};\n".format(wire.index, wire.epsr))
        file.write("real tand{}={:21.14E};\n".format(wire.index, wire.tand))
        file.write("\n")
        continue
    file.write("// Domain\n")
    file.write("real rd={:21.14E};\n".format(config.domain))
    file.write("\n")
    file.write("// Mesh Resolution\n")
    file.write("real NsAir={};\n".format(config.Ns))
    file.write("real NsCond={};\n".format(config.Ns))
    file.write("real NsDiel={};\n".format(config.Ns))
    file.write("\n")
    file.write("// Mesh\n")
    for wire in config.wires:
        file.write("border C{}(t=0, 2*pi) {{".format(wire.index))
        file.write("x=x{}+r{}*cos(t); y=y{}+r{}*sin(t);}}\n".format(wire.index,\
                   wire.index, wire.index, wire.index, wire.index))
        file.write("border CI{}(t=0, 2*pi) {{".format(wire.index))
        file.write("x=x{}+rI{}*cos(t); y=y{}+rI{}*sin(t);}}\n".format(wire.index,\
                   wire.index, wire.index, wire.index, wire.index))
        continue
    file.write("border Cd(t=0, 2*pi){x=rd*cos(t); y=rd*sin(t);}\n")
    file.write("mesh Th=buildmesh(Cd(NsAir)\n")
    for wire in config.wires:
        file.write("\t+C{}(-NsCond)+CI{}(NsDiel)\n".format(wire.index, wire.index))
        continue
    file.write("\t);\n")
    file.write("\n")
    file.write("// Solver\n")
    if view:
        file.write("plot(Th, wait=true);\n")
        pass
    file.write("fespace Vh(Th,P2);\n")
    file.write("Vh<complex> u, v;\n")
    file.write("Vh<complex> eps = eps0*(1.0+\n")
    for i in range(0, config.N):
        file.write("\t+(EPS{}*(1.0-1i*tand{})-1.0)*(sqrt((x-x{})^2+(y-y{})^2)<=rI{})\n".format(\
             i, i, i, i, i))
        continue
    file.write(");\n\n")
    file.write("complex Q;\n")
    file.write("\n")
    file.write("// BCs Vector\n")
    file.write("matrix<real> Vm;\n")
    file.write("Vm.resize(NCond, 1);\n")
    file.write("\n")
    file.write("// For Loop\n")
    file.write("for(int i=0; i<NCond; i++){\n")
    file.write("\tfor(int j=0; j<NCond; j++){\n")
    file.write("\t\tfor(int k=0; k<NCond; k++){if(k==i){Vm(k,0)=+0.5;}else{Vm(k,0)=-0.5;}}\n")
    file.write("\t\t\tsolve V(u, v)=int2d(Th)(eps*(dx(u)*dx(v)+dy(u)*dy(v)))\n")
    file.write("\t\t\t+on(C0, u=-0.5)")
    for i in range(1, config.N):
        file.write("\n\t\t\t+on(C{}, u=Vm({}, 0))".format(i, i-1))
        continue
    file.write(";\n")
    file.write("\t\t\tif (j==0)\n")
    file.write("\t\t\t\t{Q=int1d(Th, C1)(eps*(N.x*dx(u)+N.y*dy(u)));}\n")
    for i in range(2, config.N):
        file.write("\t\t\telse if (j=={})\n".format(i-1))
        file.write("\t\t\t\t{{Q=int1d(Th, C{})(eps*(N.x*dx(u)+N.y*dy(u)));}}\n".format(i))
        continue
    file.write("\t\t\tfileR << real(Q) << \" \";\n")
    file.write("\t\t\tfileI << imag(Q) << \" \";\n")
    file.write("\t\t}\n")
    file.write("\t\tfileR << endl;\n")
    file.write("\t\tfileI << endl;\n")
    file.write("}\n")
    file.write("\nexit(0);")
    file.close()
    # Execute FF Scripe File
    os.system("FreeFEM++ -nw -ne -nc {}".format(fileName))
    pass

def solveRMatrix(config, view=True):
    config.check()
    # Write FF Scripe File
    fileName = "FF/FFScriptR.edp"
    file = open(fileName, "w")
    file.write("// Output File\n")
    file.write("ofstream file(\"FF/RMatrix.dat\");\n")
    file.write("\n")
    file.write("// Constants\n")
    file.write("real eps0={:21.14E};\n".format(eps0))
    file.write("real mu0={:21.14E};\n".format(mu0))
    file.write("real freq={:21.14E};\n".format(config.freq))
    file.write("int NCond={};\n".format(config.N-1))
    file.write("\n")
    for wire in config.wires:
        file.write("// Conductor {}\n".format(wire.index))
        file.write("real x{}={:21.14E};\n".format(wire.index, wire.x))
        file.write("real y{}={:21.14E};\n".format(wire.index, wire.y))
        file.write("real r{}={:21.14E};\n".format(wire.index, wire.r))
        file.write("real rI{}={:21.14E};\n".format(wire.index, wire.rI))
        file.write("real EPS{}={:21.14E};\n".format(wire.index, wire.epsr))
        file.write("real tand{}={:21.14E};\n".format(wire.index, wire.tand))
        file.write("real sigma{}={:21.14E};\n".format(wire.index, wire.sigma))
        file.write("\n")
        continue
    file.write("// Domain\n")
    file.write("real rd={:21.14E};\n".format(config.domain))
    file.write("\n")
    file.write("// Mesh Resolution\n")
    file.write("real NsAir={};\n".format(config.Ns))
    file.write("real NsCond={};\n".format(config.Ns))
    file.write("real NsDiel={};\n".format(config.Ns))
    file.write("\n")
    file.write("// Mesh\n")
    for wire in config.wires:
        file.write("border C{}(t=0, 2*pi) {{".format(wire.index))
        file.write("x=x{}+r{}*cos(t); y=y{}+r{}*sin(t);}}\n".format(wire.index,\
                   wire.index, wire.index, wire.index, wire.index))
        file.write("border CI{}(t=0, 2*pi) {{".format(wire.index))
        file.write("x=x{}+rI{}*cos(t); y=y{}+rI{}*sin(t);}}\n".format(wire.index,\
                   wire.index, wire.index, wire.index, wire.index))
        continue
    file.write("border Cd(t=0, 2*pi){x=rd*cos(t); y=rd*sin(t);}\n")
    file.write("mesh Th=buildmesh(Cd(NsAir)\n")
    for wire in config.wires:
        file.write("\t+C{}(-NsCond)+CI{}(NsDiel)\n".format(wire.index, wire.index))
        continue
    file.write("\t);\n")
    file.write("\n")
    file.write("// Solver\n")
    if view:
        file.write("plot(Th, wait=true);\n")
        pass
    file.write("fespace Vh(Th,P2);\n")
    file.write("Vh<complex> u, v;\n")
    file.write("Vh<complex> eps = eps0*(1.0+\n")
    for i in range(0, config.N):
        file.write("\t+(EPS{}*(1.0-1i*tand{})-1.0)*(sqrt((x-x{})^2+(y-y{})^2)<=rI{})\n".format(\
             i, i, i, i, i))
        continue
    file.write(");\n\n")
    file.write("complex Rs0=sqrt(pi*freq*mu0/sigma0);\n")
    file.write("matrix<real> sigma;\n");
    file.write("sigma.resize(NCond, 1);\n");
    for i in range(1, config.N):
        file.write("sigma({}, 0)=sigma{};\n".format(i-1, i));
        continue
    file.write("matrix<real> Rs;\n");
    file.write("Rs.resize(NCond, 1);\n");
    file.write("for(int k=0; k<NCond; k++){Rs(k,0)=sqrt(pi*freq*mu0/sigma(k,0));}\n")
    file.write("complex I0, I, R;\n")
    file.write("\n")
    file.write("// BCs Vector\n")
    file.write("matrix<real> Vm;\n")
    file.write("Vm.resize(NCond, 1);\n")
    file.write("\n")
    file.write("// For Loop\n")
    file.write("for(int i=0; i<NCond; i++){\n")
    file.write("\tfor(int j=0; j<NCond; j++){\n")
    file.write("\t\tfor(int k=0; k<NCond; k++){if(k==i){Vm(k,0)=+0.5;}else{Vm(k,0)=-0.5;}}\n")
    file.write("\t\t\tsolve V(u, v)=int2d(Th)(eps*(dx(u)*dx(v)+dy(u)*dy(v)))\n")
    file.write("\t\t\t+on(C0, u=-0.5)")
    for i in range(1, config.N):
        file.write("\n\t\t\t+on(C{}, u=Vm({}, 0))".format(i, i-1))
        continue
    file.write(";\n")
    file.write("\t\t\tI0=int1d(Th, C0)(sqrt(abs(dx(u)*dx(u)+dy(u)*dy(u))));\n")
    file.write("\t\t\tif (j==0)\n")
    file.write("\t\t\t\t{I=Rs0*int1d(Th, C0)(abs(dx(u)*dx(u)+dy(u)*dy(u)))+\n")
    file.write("\t\t\t\tRs(0,0)*int1d(Th, C1)(abs(dx(u)*dx(u)+dy(u)*dy(u)));\n")
    file.write("\t\t\t\tR=I/(I0*I0);}\n")
    for i in range(2, config.N):
        file.write("\t\t\telse if (j=={})\n".format(i-1))
        file.write("\t\t\t\t{I=Rs0*int1d(Th, C0)(abs(dx(u)*dx(u)+dy(u)*dy(u)))+\n")
        file.write("\t\t\t\tRs({},0)*int1d(Th, C{})(abs(dx(u)*dx(u)+dy(u)*dy(u)));\n".format(i-1, i))
        file.write("\t\t\t\tR=I/(I0*I0);}}\n".format(i-1))
        continue
    file.write("\t\t\tfile << real(R) << \" \";\n")
    file.write("\t\t}\n")
    file.write("\t\tfile << endl;\n")
    file.write("}\n")
    file.write("\nexit(0);")
    file.close()
    # Execute FF Scripe File
    os.system("FreeFEM++ -nw -ne -nc {}".format(fileName))
    pass

def computeParameters(config):
    solveC0Matrix(config, view=False)
    solveCMatrix(config, view=False)
    solveRMatrix(config, view=False)
    C0file = 'FF/C0Matrix.dat'
    C1fileR = 'FF/CMatrixR.dat'
    C1fileI = 'FF/CMatrixI.dat'
    Rfile = 'FF/RMatrix.dat'
    N = config.N
    freq = config.freq
    # Read Raw Data
    C0 = np.ndarray((N-1, N-1))
    count = 0
    with open(C0file, 'r') as file:
        for newLine in file:
            newLine = newLine.split()
            for i in range(0, N-1):
                C0[count, i] = float(newLine[i])
                continue
            count+=1
            continue
        pass
    CR = np.ndarray((N-1, N-1))
    count = 0
    with open(C1fileR, 'r') as file:
        for newLine in file:
            newLine = newLine.split()
            for i in range(0, N-1):
                CR[count, i] = float(newLine[i])
                continue
            count+=1
            continue
        pass
    CI = np.ndarray((N-1, N-1))
    count = 0
    with open(C1fileI, 'r') as file:
        for newLine in file:
            newLine = newLine.split()
            for i in range(0, N-1):
                CI[count, i] = float(newLine[i])
                continue
            count+=1
            continue
        pass
    R = np.ndarray((N-1, N-1))
    count = 0
    with open(Rfile, 'r') as file:
        for newLine in file:
            newLine = newLine.split()
            for i in range(0, N-1):
                R[count, i] = float(newLine[i])
                continue
            count+=1
            continue
        pass
    # Final Results
    C = CR
    L = (1.0/(c0**2))*np.linalg.inv(C0)
    G = -2.0*np.pi*freq*CI
    # Save Results
    np.savetxt('Data/C.dat', C, fmt='%22.14E', delimiter=' ')
    np.savetxt('Data/L.dat', L, fmt='%22.14E', delimiter=' ')
    np.savetxt('Data/G.dat', G, fmt='%22.14E', delimiter=' ')
    np.savetxt('Data/R.dat', R, fmt='%22.14E', delimiter=' ')
    pass