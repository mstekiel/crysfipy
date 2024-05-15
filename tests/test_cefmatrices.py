from crysfipy.cefmatrices import *

#######################################################################
# Old implementation as separate matrices

# Stevens operators
# Cross checked with the McPhase manual
# https://www2.cpfs.mpg.de/~rotter/homepage_mcphase/manual/node132.html

def O20(J):
    JJ = J*(J+1)
    J2p1 = int(2*J + 1)

    Jz = J_z(J)
    E = eye(J2p1)

    return 3*Jz*Jz - JJ * E

def O22(J):
    Jplus = J_plus(J)
    Jminus = J_minus(J)

    return 0.5 * (mp(Jplus, 2) + mp(Jminus, 2))
    
def O2m2(J):
    Jplus = J_plus(J)
    Jminus = J_minus(J)

    return -0.5j * (mp(Jplus, 2) - mp(Jminus, 2))

def O40(J):
    JJ = J*(J+1)
    J2p1 = int(2*J + 1)

    Jz = J_z(J)
    E = eye(J2p1)

    return 35 * mp(Jz, 4) + (25 - 30 * JJ)*mp(Jz, 2) + E * JJ * (3 * JJ - 6)


def O42(J):
    JJ = J*(J+1)
    J2p1 = int(2*J + 1)

    Jz = J_z(J)
    Jplus = J_plus(J)
    Jminus = J_minus(J)
    E = eye(J2p1)

    # helper matrices:	
    M_1 = 7 * mp(Jz, 2) - E*(JJ + 5)
    M_2 = mp(Jplus, 2) + mp(Jminus, 2)

    return 0.25 * (dot(M_1, M_2) + dot(M_2, M_1))

def O4m2(J):
    JJ = J*(J+1)
    J2p1 = int(2*J + 1)

    Jz = J_z(J)
    Jplus = J_plus(J)
    Jminus = J_minus(J)
    E = eye(J2p1)

    # helper matrices:	
    M_1 = 7 * mp(Jz, 2) - E*(JJ + 5)
    M_2 = mp(Jplus, 2) - mp(Jminus, 2)

    return -0.25j * (dot(M_1, M_2) + dot(M_2, M_1))


def O43(J):
    JJ = J*(J+1)
    J2p1 = int(2*J + 1)

    Jz = J_z(J)
    Jplus = J_plus(J)
    Jminus = J_minus(J)

    # helper matrices:	
    M_1 = Jz
    M_2 = mp(Jplus, 3) + mp(Jminus, 3)
    
    return 0.25 * (dot(M_1, M_2) + dot(M_2, M_1))
    
def O4m3(J):
    JJ = J*(J+1)
    J2p1 = int(2*J + 1)

    Jz = J_z(J)
    Jplus = J_plus(J)
    Jminus = J_minus(J)

    # helper matrices:	
    M_1 = Jz
    M_2 = mp(Jplus, 3) - mp(Jminus, 3)
    
    return -0.25j * (dot(M_1, M_2) + dot(M_2, M_1))


def O44(J):
    JJ = J*(J+1)
    J2p1 = int(2*J + 1)

    Jplus = J_plus(J)
    Jminus = J_minus(J)
    
    return 0.5 * (mp(Jplus, 4) + mp(Jminus, 4))
    
def O4m4(J):
    JJ = J*(J+1)
    J2p1 = int(2*J + 1)

    Jplus = J_plus(J)
    Jminus = J_minus(J)
    
    return -0.5j * (mp(Jplus, 4) - mp(Jminus, 4))


def O60(J):
    JJ = J*(J+1)
    J2p1 = int(2*J + 1)

    Jz = J_z(J)
    E = eye(J2p1)

    return 231 * mp(Jz, 6) + mp(Jz, 4) * (735 - 315 * JJ) + mp(Jz, 2)*(105*JJ**2 - 525*JJ + 294) + E * (-5*JJ**3 + 40 * JJ**2 - 60*JJ)

def O62(J):
    JJ = J*(J+1)
    J2p1 = int(2*J + 1)
    
    E = eye(J2p1)
    Jz = J_z(J)
    Jplus = J_plus(J)
    Jminus = J_minus(J)
    
    # helper matrices:	
    M_1 = 33 * mp(Jz, 4) - mp(Jz, 2) * (18 * JJ + 123) + E * (JJ**2 + 10*JJ + 102)
    M_2 = mp(Jplus, 2) + mp(Jminus, 2)

    return 0.25 * (dot(M_1, M_2) + dot(M_2, M_1))

def O6m2(J):
    JJ = J*(J+1)
    J2p1 = int(2*J + 1)
    
    E = eye(J2p1)
    Jz = J_z(J)
    Jplus = J_plus(J)
    Jminus = J_minus(J)
    
    # helper matrices:	
    M_1 = 33 * mp(Jz, 4) - mp(Jz, 2) * (18 * JJ + 123) + E * (JJ**2 + 10*JJ + 102)
    M_2 = mp(Jplus, 2) - mp(Jminus, 2)

    return -0.25j * (dot(M_1, M_2) + dot(M_2, M_1))
    
def O63(J):
    JJ = J*(J+1)
    J2p1 = int(2*J + 1)

    Jz = J_z(J)
    Jplus = J_plus(J)
    Jminus = J_minus(J)

    # helper matrices:	
    M_1 = 11 * mp(Jz, 3) - Jz * (59 + 3*JJ)
    M_2 = mp(Jplus, 3) + mp(Jminus, 3)
    
    return 0.25 * (dot(M_1, M_2) + dot(M_2, M_1))

def O6m3(J):
    JJ = J*(J+1)
    J2p1 = int(2*J + 1)

    Jz = J_z(J)
    Jplus = J_plus(J)
    Jminus = J_minus(J)

    # helper matrices:	
    M_1 = 11 * mp(Jz, 3) - Jz * (59 + 3*JJ)
    M_2 = mp(Jplus, 3) - mp(Jminus, 3)
    
    return -0.25j * (dot(M_1, M_2) + dot(M_2, M_1))
    
def O64(J):
    JJ = J*(J+1)
    J2p1 = int(2*J + 1)

    E = eye(J2p1)
    Jz = J_z(J)
    Jplus = J_plus(J)
    Jminus = J_minus(J)

    # helper matrices:	
    M_1 = 11 * mp(Jz, 2) - E * (JJ + 38)
    M_2 = mp(Jplus, 4) + mp(Jminus, 4)

    return 0.25 * (dot(M_1, M_2) + dot(M_2, M_1))

def O6m4(J):
    JJ = J*(J+1)
    J2p1 = int(2*J + 1)

    E = eye(J2p1)
    Jz = J_z(J)
    Jplus = J_plus(J)
    Jminus = J_minus(J)

    # helper matrices:	
    M_1 = 11 * mp(Jz, 2) - E * (JJ + 38)
    M_2 = mp(Jplus, 4) - mp(Jminus, 4)

    return -0.25j * (dot(M_1, M_2) + dot(M_2, M_1))

def O66(J):
    JJ = J*(J+1)
    J2p1 = int(2*J + 1)

    Jplus = J_plus(J)
    Jminus = J_minus(J)
    
    return 0.5 * (mp(Jplus, 6) + mp(Jminus, 6))

def O6m6(J):
    JJ = J*(J+1)
    J2p1 = int(2*J + 1)

    Jplus = J_plus(J)
    Jminus = J_minus(J)
    
    return -0.5j * (mp(Jplus, 6) - mp(Jminus, 6))

#######################################################################
# PyCrystalField implementation, very slow
class Operator():
    def __init__(self, J):
        self.O = np.zeros((int(2*J+1), int(2*J+1)))
        self.m = np.arange(-J,J+1,1)
        self.j = J

    @staticmethod
    def Jz(J):
        obj = Operator(J)
        for i in range(len(obj.O)):
            for k in range(len(obj.O)):
                if i == k:
                    obj.O[i,k] = (obj.m[k])
        return obj

    @staticmethod
    def Jplus(J):
        obj = Operator(J)
        for i in range(len(obj.O)):
            for k in range(len(obj.O)):
                if k+1 == i:
                    obj.O[i,k] = np.sqrt((obj.j-obj.m[k])*(obj.j+obj.m[k]+1))
        return obj

    @staticmethod
    def Jminus(J):
        obj = Operator(J)
        for i in range(len(obj.O)):
            for k in range(len(obj.O)):
                if k-1 == i:
                    obj.O[i,k] = np.sqrt((obj.j+obj.m[k])*(obj.j-obj.m[k]+1))
        return obj

    @staticmethod
    def Jx(J):
        objp = Operator.Jplus(J)
        objm = Operator.Jminus(J)
        return 0.5*objp + 0.5*objm

    @staticmethod
    def Jy(J):
        objp = Operator.Jplus(J)
        objm = Operator.Jminus(J)
        return -0.5j*objp + 0.5j*objm

    def __add__(self,other):
        newobj = Operator(self.j)
        if isinstance(other, Operator):
            newobj.O = self.O + other.O
        else:
            newobj.O = self.O + other*np.identity(int(2*self.j+1))
        return newobj

    def __radd__(self,other):
        newobj = Operator(self.j)
        if isinstance(other, Operator):
            newobj.O = self.O + other.O
        else:
            newobj.O = self.O + other*np.identity(int(2*self.j+1))
        return newobj

    def __sub__(self,other):
        newobj = Operator(self.j)
        if isinstance(other, Operator):
            newobj.O = self.O - other.O
        else:
            newobj.O = self.O - other*np.identity(int(2*self.j+1))
        return newobj

    def __mul__(self,other):
        newobj = Operator(self.j)
        if (isinstance(other, int) or isinstance(other, float) or isinstance(other, complex)):
            newobj.O = other * self.O
        else:
            newobj.O = np.dot(self.O, other.O)
        return newobj

    def __rmul__(self,other):
        newobj = Operator(self.j)
        if (isinstance(other, int) or isinstance(other, float)  or isinstance(other, complex)):
            newobj.O = other * self.O
        else:
            newobj.O = np.dot(other.O, self.O)
        return newobj

    def __pow__(self, power):
        newobj = Operator(self.j)
        newobj.O = self.O
        for i in range(power-1):
            newobj.O = np.dot(newobj.O,self.O)
        return newobj

    def __neg__(self):
        newobj = Operator(self.j)
        newobj.O = -self.O
        return newobj

    def __repr__(self):
        return repr(self.O)

def StevensOpPCF(J,n,m):
    """generate stevens operator for a given total angular momentum
    and a given n and m state"""
    Jz = Operator.Jz(J=J)
    Jp = Operator.Jplus(J=J)
    Jm = Operator.Jminus(J=J)
    X = J*(J+1.)

    if [n,m] == [0,0]:
        return np.zeros((int(2*J+1), int(2*J+1)))
    elif [n,m] == [1,0]:
        matrix = Jz
    elif [n,m] == [1,1]:
        matrix = 0.5 *(Jp + Jm)
    elif [n,m] == [1,-1]:
        matrix = -0.5j *(Jp - Jm)

    elif [n,m] == [2,2]:
        matrix = 0.5 *(Jp**2 + Jm**2)
    elif [n,m] == [2,1]:
        matrix = 0.25*(Jz*(Jp + Jm) + (Jp + Jm)*Jz)
    elif [n,m] == [2,0]:
        matrix = 3*Jz**2 - X
    elif [n,m] == [2,-1]:
        matrix = -0.25j*(Jz*(Jp - Jm) + (Jp - Jm)*Jz)
    elif [n,m] == [2,-2]:
        matrix = -0.5j *(Jp**2 - Jm**2)

    elif [n,m] == [3,3]:
        matrix = 0.5 *(Jp**3 + Jm**3)
    elif [n,m] == [3,2]:
        matrix = 0.25 *((Jp**2 + Jm**2)*Jz + Jz*(Jp**2 + Jm**2))
    elif [n,m] == [3,1]:
        matrix = 0.25*((Jp + Jm)*(5*Jz**2 - X - 0.5) + (5*Jz**2 - X - 0.5)*(Jp + Jm))
    elif [n,m] == [3,0]:
        matrix = 5*Jz**3 - (3*X-1)*Jz
    elif [n,m] == [3,-1]:
        matrix = -0.25j*((Jp - Jm)*(5*Jz**2 - X - 0.5) + (5*Jz**2 - X - 0.5)*(Jp - Jm))
    elif [n,m] == [3,-2]:
        matrix = -0.25j*(Jz*(Jp**2 - Jm**2) + (Jp**2 - Jm**2)*Jz)
    elif [n,m] == [3,-3]:
        matrix = -0.5j *(Jp**3 - Jm**3)

    elif [n,m] == [4,4]:
        matrix = 0.5 *(Jp**4 + Jm**4)
    elif [n,m] == [4,3]:
        matrix = 0.25 *((Jp**3 + Jm**3)*Jz + Jz*(Jp**3 + Jm**3))
    elif [n,m] == [4,2]:
        matrix = 0.25 *((Jp**2 + Jm**2)*(7*Jz**2 -X -5) + (7*Jz**2 -X -5)*(Jp**2 + Jm**2))
    elif [n,m] == [4,1]:
        matrix = 0.25 *((Jp + Jm)*(7*Jz**3 -(3*X+1)*Jz) + (7*Jz**3 -(3*X+1)*Jz)*(Jp + Jm))
    elif [n,m] == [4,0]:
        matrix = 35*Jz**4 - (30*X -25)*Jz**2 + 3*X**2 - 6*X
    elif [n,m] == [4,-4]:
        matrix = -0.5j *(Jp**4 - Jm**4)
    elif [n,m] == [4,-3]:
        matrix = -0.25j *((Jp**3 - Jm**3)*Jz + Jz*(Jp**3 - Jm**3))
    elif [n,m] == [4,-2]:
        matrix = -0.25j *((Jp**2 - Jm**2)*(7*Jz**2 -X -5) + (7*Jz**2 -X -5)*(Jp**2 - Jm**2))
    elif [n,m] == [4,-1]:
        matrix = -0.25j *((Jp - Jm)*(7*Jz**3 -(3*X+1)*Jz) + (7*Jz**3 -(3*X+1)*Jz)*(Jp - Jm))

    elif [n,m] == [6,6]:
        matrix = 0.5 *(Jp**6 + Jm**6)
    elif [n,m] == [6,5]:
        matrix = 0.25*((Jp**5 + Jm**5)*Jz + Jz*(Jp**5 + Jm**5))
    elif [n,m] == [6,4]:
        matrix = 0.25*((Jp**4 + Jm**4)*(11*Jz**2 -X -38) + (11*Jz**2 -X -38)*(Jp**4 + Jm**4))
    elif [n,m] == [6,3]:
        matrix = 0.25*((Jp**3 + Jm**3)*(11*Jz**3 -(3*X+59)*Jz) + (11*Jz**3 -(3*X+59)*Jz)*(Jp**3 + Jm**3))
    elif [n,m] == [6,2]:
        matrix = 0.25*((Jp**2 + Jm**2)*(33*Jz**4 -(18*X+123)*Jz**2 +X**2 +10*X +102) +\
                    (33*Jz**4 -(18*X+123)*Jz**2 +X**2 +10*X +102)*(Jp**2 + Jm**2))
    elif [n,m] == [6,1]:
        matrix = 0.25*((Jp +Jm)*(33*Jz**5 -(30*X-15)*Jz**3 +(5*X**2 -10*X +12)*Jz) +\
                    (33*Jz**5 -(30*X-15)*Jz**3 +(5*X**2 -10*X +12)*Jz)*(Jp+ Jm))
    elif [n,m] == [6,0]:
        matrix = 231*Jz**6 - (315*X-735)*Jz**4 + (105*X**2 -525*X +294)*Jz**2 -\
                5*X**3 + 40*X**2 - 60*X
    elif [n,m] == [6,-6]:
        matrix = -0.5j *(Jp**6 - Jm**6)
    elif [n,m] == [6,-5]:
        matrix = -0.25j*((Jp**5 - Jm**5)*Jz + Jz*(Jp**5 - Jm**5))
    elif [n,m] == [6,-4]:
        matrix = -0.25j*((Jp**4 - Jm**4)*(11*Jz**2 -X -38) + (11*Jz**2 -X -38)*(Jp**4 - Jm**4))
    elif [n,m] == [6,-3]:
        matrix = -0.25j*((Jp**3 - Jm**3)*(11*Jz**3 -(3*X+59)*Jz) + (11*Jz**3 -(3*X+59)*Jz)*(Jp**3 - Jm**3))
    elif [n,m] == [6,-2]:
        matrix = -0.25j*((Jp**2 - Jm**2)*(33*Jz**4 -(18*X+123)*Jz**2 +X**2 +10*X +102) +\
                    (33*Jz**4 -(18*X+123)*Jz**2 +X**2 +10*X +102)*(Jp**2 - Jm**2))
    elif [n,m] == [6,-1]:
        matrix = -0.25j*((Jp - Jm)*(33*Jz**5 -(30*X-15)*Jz**3 +(5*X**2 -10*X +12)*Jz) +\
                (33*Jz**5 -(30*X-15)*Jz**3 +(5*X**2 -10*X +12)*Jz)*(Jp - Jm))
    else:
        raise KeyError(f'Index {(n,m)} not implemented')

    return matrix.O.T

###########################
def test_matrices(J):
    O = StevensBase(J)
    equivalent_ops = {
        'O11':  (J_x(J), O.O11(),  StevensOpPCF(J,1, 1)) ,
        'O10':  (J_z(J), O.O10(),  StevensOpPCF(J,1, 0)) ,
        'O1m1': (J_y(J), O.O1m1(), StevensOpPCF(J,1,-1)) ,

        'O22': (O22(J), O.O22(), StevensOpPCF(J,2,2)) ,
        'O21': (None, O.O21(), StevensOpPCF(J,2,1)) ,
        'O20': (O20(J), O.O20(), StevensOpPCF(J,2,0)) ,
        'O2m1': (None, O.O2m1(), StevensOpPCF(J,2,-1)) ,
        'O2m2': (O2m2(J), O.O2m2(), StevensOpPCF(J,2,-2)) ,

        'O33': (None, O.O33(), StevensOpPCF(J,3,3)) ,
        'O32': (None, O.O32(), StevensOpPCF(J,3,2)) ,
        'O31': (None, O.O31(), StevensOpPCF(J,3,1)) ,
        'O30': (None, O.O30(), StevensOpPCF(J,3,0)) ,
        'O3m1': (None, O.O3m1(), StevensOpPCF(J,3,-1)) ,
        'O3m2': (None, O.O3m2(), StevensOpPCF(J,3,-2)) ,
        'O3m3': (None, O.O3m3(), StevensOpPCF(J,3,-3)) ,

        'O44': (O44(J), O.O44(), StevensOpPCF(J,4,4)) ,
        'O43': (O43(J), O.O43(), StevensOpPCF(J,4,3)) ,
        'O42': (O42(J), O.O42(), StevensOpPCF(J,4,2)) ,
        'O41': (None,    O.O41(), StevensOpPCF(J,4,1)) ,
        'O40': (O40(J), O.O40(), StevensOpPCF(J,4,0)) ,
        'O4m1': (None,     O.O4m1(), StevensOpPCF(J,4,-1)) ,
        'O4m2': (O4m2(J), O.O4m2(), StevensOpPCF(J,4,-2)) ,
        'O4m3': (O4m3(J), O.O4m3(), StevensOpPCF(J,4,-3)) ,
        'O4m4': (O4m4(J), O.O4m4(), StevensOpPCF(J,4,-4)) ,


        'O66':  (O66(J), O.O66(), StevensOpPCF(J,6,6)) ,
        'O65':  (None,    O.O65(), StevensOpPCF(J,6,5)) ,
        'O64':  (O64(J), O.O64(), StevensOpPCF(J,6,4)) ,
        'O63':  (O63(J), O.O63(), StevensOpPCF(J,6,3)) ,
        'O62':  (O62(J), O.O62(), StevensOpPCF(J,6,2)) ,
        'O61':  (None,    O.O61(), StevensOpPCF(J,6,1)) ,
        'O60':  (O60(J), O.O60(), StevensOpPCF(J,6,0)) ,
        'O6m1': (None,     O.O6m1(), StevensOpPCF(J,6,-1)) ,
        'O6m2': (O6m2(J), O.O6m2(), StevensOpPCF(J,6,-2)) ,
        'O6m3': (O6m3(J), O.O6m3(), StevensOpPCF(J,6,-3)) ,
        'O6m4': (O6m4(J), O.O6m4(), StevensOpPCF(J,6,-4)) ,
        'O6m5': (None,     O.O6m5(), StevensOpPCF(J,6,-5)) ,
        'O6m6': (O6m6(J), O.O6m6(), StevensOpPCF(J,6,-6)) ,
    }

    results = {}
    flag = True
    for Oname, (Oold, Onew, Opcf) in equivalent_ops.items():
        atol = np.max(Onew)*1e-10	# effectively makes it into rtol
        # print(J,Oname)
        # print(Oold, Onew, Opcf)
        if Oold is not None:
            if not np.allclose(Onew, Oold, atol=atol):
                print(f'\tOnew does not match Oold for {Oname} J={J}')
                flag = False

        if Opcf is not None:
            tp = np.allclose(Onew, Opcf, atol=atol)
            tm = np.allclose(Onew,-Opcf, atol=atol)
            if not (tp or tm):
                print(f'\tOnew does not match Opcf for {Oname} J={J}')
                flag = False
            # if (tp and tm):
            # 	print(f'\tOnew matches Ocf +- for {Oname} J={J}')

    return flag
    
#########################
def test_execution(Ntot):
    J=5
    O = StevensBase(J)

    N = 0
    t_start = time()
    while N < Ntot:
        _ = O['O42']
        N += 1

    t_end = time()
    print(f't1 took {t_end-t_start} s')

    N = 0
    t_start = time()
    while N < Ntot:
        _ = StevensOpPCF(J,4,2)
        N += 1

    t_end = time()
    print(f't2 took {t_end-t_start} s')

def test_implementation():
    StevensBase.implementedOps()

if __name__ == '__main__':
    O = StevensBase(2)
    print(O['O21'])
    print(O.O21())

    from time import time
    import numpy as np

    # for Jh in np.arange(1, 101):
    # 	result = test_implementation(Jh/2)
    # 	print(f'J={Jh/2} test: {result}')


    test_execution(10_000)

