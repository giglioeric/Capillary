###############################################################################
#
#   Author :  GIGLIo ERIC   (CIMAP / CAEN / FRANCE)
#   On ResearchGate:  https://www.researchgate.net/profile/Eric-Giglio-2

#   More detatils in article :   https://doi.org/10.48550/arXiv.2401.13521

#   Reads the input file 'parameter.txt'
#  
#   Generates the output file  'coefficients.txt'
#
#
###############################################################################

# check inputs
def check_inputs(parameters):

        
    print("Checking validity of input list \n")
    traceback_i = 0
    for item in parameters:   
        traceback_i += 1
        

    if (traceback_i !=10) :          # Checks list length
        print("Error !: Input file must have 10 entriees, but has {}\n".format(traceback_i))
        exit()
    
   
    x= parameters[0].split();r1 = float(x[0])
    x= parameters[1].split();r2 = float(x[0])
    x= parameters[2].split();r3 = float(x[0])
    x= parameters[3].split();h = float(x[0])
    x= parameters[4].split();er = float(x[0])
    x= parameters[5].split();kb = float(x[0])
    x= parameters[6].split();ks = float(x[0])
    x= parameters[7].split();mMax = int(x[0])
    x= parameters[8].split();nMax = int(x[0])
    x= parameters[9].split();flag = int(x[0])
    

    # Check additional conditions
    if r1 <= 0:
        print("Error !: R1 <= 0  : R1= {} m".format(r1))
        exit()
    if (r2-r1) <= 0:
        print("Error !: R2 <= R1 : R2 = {} m".format(r2))
        exit()
    if (r3-r2) <= 0:
        print("Error !: R3 <= R2 : R3 = {} m".format(r3))
        exit()
    if h <= 0:
        print("Error !: H <= 0  : H= {} m".format(h))
        exit()

    if er < 1:
        print("Error !: er < 1  : er = {}  ".format(e))
        exit()

    if  (powerof2(mMax)==0) :
        print("Error !: M={} not a power of 2".format(mMax))
        exit()
        
    if  (powerof2(nMax)==0) :
        print("Error !: N={} not a power of 2".format(nMax))
        exit()

    if (flag==1) :
        print("Using Dirichlet (absorbing) boundary conditions at z=H")

    elif (flag==2) :
        print("Using Neumann (zero-flux) boundary conditions at z=H")
    else :
        print("Error ! :flag={} invalid : must be 1 or 2".format(flag))
        exit()
             
    return
 

###############################################################################
# check power of 2

def powerof2(n) :
    if ((n-int(n))!=0) : return False
    if (n == 1) or (n == 2) :  return True
    return (int(bin(n)[3:]) == 0)
            
###############################################################################
# Calculates all the coefficients 0<= m < Mmax, 1<= n <= Nmax and store in 'coefficients.txt'

def build_matrix(parameters):

   
    x= parameters[0].split();r1 = float(x[0])
    x= parameters[1].split();r2 = float(x[0])
    x= parameters[2].split();r3 = float(x[0])
    x= parameters[3].split();h = float(x[0])
    x= parameters[4].split();er = float(x[0])
    x= parameters[5].split();kb = float(x[0])
    x= parameters[6].split();ks = float(x[0])
    x= parameters[7].split();mMax = int(x[0])
    x= parameters[8].split();nMax = int(x[0])
    x= parameters[9].split();flag = int(x[0])

    e0 = 8.8e-12
    
    
    print("mMax = {}  :  nMax = {} \n".format(mMax,nMax))

    matrix = np.array([[0.0,0.0,0.0,0,0],[0,0,0.0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0.0,0,0,0]])
    invmat = np.array([[0.0,0.0,0.0,0,0],[0,0,0.0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0.0,0,0,0]])
    fmn =  np.array([[0.0,0.0],[0.0,0.0]])
    tmn =  np.array([[0.0,0.0],[0.0,0.0]])
    pmn =  np.array([[0.0,0.0],[0.0,0.0]])
    qmn =  np.array([[0.0,0.0],[0.0,0.0]])
    taumn = np.array([0.0,0.0])
    s1mn = np.array([0,0,0,1/e0,0])
    s2mn = np.array([0,0,0,0,1/e0])
    pp1mn = np.array([[1.0,0.0],[0.0,0.0]])
    invtmn = np.array([[0.0,0.0],[0.0,0.0]])

    fout = open("coefficients.txt", "w")
    print("--------------(m,n)-----------------",file=fout)
    print("amn \t \t \t tau^(1)_mn \t \t p11 \t \t \t p12",file=fout)
    print("bmn \t \t \t tau^(2)_mn \t \t p21 \t \t \t p22",file=fout)
   
    
    for m in range(mMax):    
        for n in range(nMax):
            nn=n+1   # because n starts from 0 while we want it to start from 1
            if flag==1 :
                knn = nn*np.pi/h   # absorbing B. C.
            else :
                knn = (nn-0.5)*np.pi/h  # blocking B. C.
                
            i1 = sp.iv(m,knn*r1)     # BesselI()
            i2 = sp.iv(m,knn*r2)
            i3 = sp.iv(m,knn*r3)
            
            k1 = sp.kn(m,knn*r1)     #BesselK()
            k2 = sp.kn(m,knn*r2)
            k3 = sp.kn(m,knn*r3)
                      
            i1p = knn*sp.ivp(m,knn*r1)   # dBesselI()/dr
            i2p = knn*sp.ivp(m,knn*r2)
            
            k1p = knn*sp.kvp(m,knn*r1)   # DBesselK()/dr
            k2p = knn*sp.kvp(m,knn*r2)
            
           
            matrix[0][0]=i1
            matrix[0][1]=-i1
            matrix[0][2]=-k1
            matrix[0][3]=0
            matrix[0][4]=0
            
            matrix[1][0]=0
            matrix[1][1]=i2
            matrix[1][2]=k2
            matrix[1][3]=-i2
            matrix[1][4]=-k2
            
            matrix[2][0]=0
            matrix[2][1]=0
            matrix[2][2]=0
            matrix[2][3]=i3
            matrix[2][4]=k3
            
            matrix[3][0]=i1p
            matrix[3][1]=-i1p*er
            matrix[3][2]=-k1p*er
            matrix[3][3]=0
            matrix[3][4]=0
            
            matrix[4][0]=0
            matrix[4][1]=i2p*er
            matrix[4][2]=k2p*er
            matrix[4][3]=-i2p
            matrix[4][4]=-k2p
            
            
            invmat = np.linalg.inv(matrix)   # inverse matrix
            
            
            a1mn = invmat[0][3]/e0;    # select element of inverse matrix
            a2mn = invmat[0][4]/e0;
            b1mn = invmat[1][3]/e0;
            b2mn = invmat[1][4]/e0;
            c1mn = invmat[2][3]/e0;
            c2mn = invmat[2][4]/e0;
            
            
            u1mn = m*m/(r1*r1)+knn*knn
            u2mn = m*m/(r2*r2)+knn*knn
            
            f11 = (-kb*i1p + ks*u1mn* i1)*b1mn + (-kb*k1p + ks*u1mn*k1)*c1mn
            f12 = (-kb*i1p + ks*u1mn *i1)*b2mn + (-kb*k1p + ks*u1mn*k1)*c2mn
            
            f21 = (kb*i2p + ks*u2mn* i2)*b1mn + (kb*k2p + ks*u2mn*k2)*c1mn
            f22 = (kb*i2p + ks*u2mn* i2)*b2mn + (kb*k2p + ks*u2mn*k2)*c2mn
            
            fmn[0][0]= f11
            fmn[0][1]= f12
            fmn[1][0]= f21
            fmn[1][1]= f22
       
            taumn,tmn=  np.linalg.eig(fmn)
            
            invtmn =  np.linalg.inv(tmn)
            
            zmn = np.dot(pp1mn,invtmn)
            pmn = np.dot(tmn,zmn)
            
            print("--------------({},{})-----------------".format(m,nn),file=fout)
            print("{} \t {}  \t {}  \t {}".format(a1mn,1/taumn[0],pmn[0][0],pmn[0][1]),file=fout)
            print("{} \t {}  \t {}  \t {}".format(a2mn,1/taumn[1],pmn[1][0],pmn[1][1]),file=fout)
            
    print("Done")
    print("File 'coefficients.txt' has been  generated according to the structure")
    print("--------------(m,n)-----------------")
    print("amn \t \t \t tau^(1)_mn \t \t p11 \t \t \t p12")
    print("bmn \t \t \t tau^(2)_mn \t \t p21 \t \t \t p22")
    fout.close()

###############################################################################
#  MAIN

if __name__ == "__main__":
    print("\n**********************************************************************\n")
    print("The Python3 routine 'rates.py' tries to read the input file 'parameters.dat'")
    print("If successful, it checks the validity of the input parameters")
    print("It generates then the output file 'coefficients.txt'\n")
  #  print("python3, numpy , scipy.special must be installed \n")
    print("The dimensions of the capillary (R1,R2,R3,H) are given in meters")
    print("'er' is the relative permittivity of the dielectric")
    print("The bulk conductivity 'kb' is given in S/m")
    print("The surface conductivity 'ks' is given in S")
    print("'M' must be a power of 2 and is the maximum number of angular moments")   
    print("'N' must be a power of 2 and is the maximum number of wave numbers")  
    print("'flag 1' means absorbing boundary condition at z=H , kn=nPi/H")
    print("'flag 2' means blocking  boundary condition at z=H, kn=(n-1/2)Pi/H\n ")
    print("**********************************************************************")
    print("Example of a valid input file 'parameters.txt' :")
    print("1.8e-3   R1")
    print("3.0e-3   R2")
    print("4.0e-3   R3")
    print("45.0e-3  H")
    print("4.6   er")
    print("3.0e-13   kb")
    print("1.0e-12   ks")
    print("16   M")
    print("256  N")
    print("1  flag : 1= Dirichlet; 2 = Neumann")
    print("----------------------\n")
    print("Reading 10 parameters from file 'parameters.txt'\n")
   
    import numpy as np
    import scipy.special as sp

    parameters = []
    try:
        with open("parameters.txt") as file :
            parameters = [line.strip() for line in file]
    except FileNotFoundError:
        print("The input file 'parameters.txt' does not exist.\nPlease create the file 'parameters.txt' with 10 entrees as indicated above")
        exit()
       
    check_inputs(parameters)
    build_matrix(parameters)
    
   
