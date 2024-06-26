from typing import Tuple

from .ion import  Ion
from .cefpars import CEFpars
from . import constants as C
from .cefmatrices import StevensBase, J_x, J_y, J_z

import numpy as np
from scipy.linalg import schur
# from numpy import conj, transpose, dot, diag
# import numbers



class CEFion:
    r'''
    Object representing a rare-earth ion in CF potential. It is internally calculated in the meV units.
    
    Parameters:
        ion : string
            Name of the ion. They are tabulated in :obj:`ion` with their parameters.
        Hfield : 1D array of floats
            External magnetic field in *T* units.
        cfp : ``crysfipy.CEFpars``
            Crystal field parameters
        diagonalize : bool, optional
            If true (default) then it automatically diagonalizes Hamiltonian, calculates energy levels and sorts all matrices so that the first eigenvector corresponds to the lowest energy level and the last one to the highest.

    Examples:
        
        TODO
        
    Attributes:
        ion : ``crysfipy:ion``
            The ``ion`` object that represents an isolated rare-earth ion.
        Jval : int/2
            The J value corresponding to the L+S quantum numbers. Not to be confused with the ``J`` operator.
        Hfield : array_like
            Vector in real space corresponding to the external magnetic field in T units.
        cfp : ``crysfipy.reion.cfp``
            Crystal field parameters
        hamiltonian : ndarray
            Hamiltonian operator. :math:`\hat{\mathcal{H}} = \sum_{ij} B_i^j \hat{O}_i^j + g_J (H_x \hat{J}_x + H_y \hat{J}_y + H_z \hat{J}_z)`
        Jx, Jy, Jz : ndarray
            Matrices representing the prinicpal quantum operators :math:`\hat{J}_\\alpha`` with :math:`\\alpha=[x,y,z]`.
        J : ndarray
            Total angular momentum vector operator :math:`\hat{J}=[\hat{J}_x, \hat{J}_y, \hat{J}_z]`. It's multidimensionality takes some time to get used to.
        energies : ndarray
            Eigenenergies of the Hamiltonian.
        eigenvectors : ndarray
            Eigenvectors of the Hamiltonian.
        degeneracies : list
            List containing entries like [energy, degeneracy], which describes how degenerated is the given energy level.
        freeionkets : list
            List of kets corresponding to the J basis of the free ion problem.
    '''

    def __init__(self, ion: Ion, cfp: CEFpars,  Hfield: Tuple[float, float, float]=[0,0,0], diagonalize: bool=True):
        self.ion = ion
        self.Jval = self.ion.J
        Jval = self.Jval
        
        self.Hfield = np.array(Hfield)
        self.cfp = cfp

        # Prepare the rest of the fields based on the main input parameters

        # Main Ji matrices
        self.Jx = J_x(Jval) 
        self.Jy = J_y(Jval)
        self.Jz = J_z(Jval)
        self.J = [self.Jx, self.Jy, self.Jz]
        
        # Matrix with projection of moments to x, y, z directions for all J2p1 levels
        mu = self.ion.gJ * np.einsum('ijj->ji', self.J)
        self.moment = np.real(mu + mu.conj())/2

        
        # Prepare a list os kets that form the basis
        if Jval%1==0:
            self.freeionkets = [f'|{int(x)}>' for x in np.linspace(Jval,-Jval,int(2*Jval+1))]
        elif Jval%1==0.5:
            self.freeionkets = [f'|{int(2*x):d}/2>' for x in np.linspace(Jval,-Jval,int(2*Jval+1))]
        
        
        # THE HAMILTONIAN
        H = -C.uB * self.ion.gJ * np.einsum('ijk,i', self.J, self.Hfield)
        
        # Store Stevens operators in the dictionary containing pointers to functions
        O = StevensBase(self.Jval)
        for Bname, Bvalue in self.cfp.Bpars.items():
            H += Bvalue * O[Bname[1:]]
            

            
        # Might take cc of everything to get rid of complex eigenvalues.
        # However, it should not be a problem
        self.hamiltonian = H
        # self.hamiltonian = (H + H.conj())/2

        # Uncomputed, empty fields
        self.eigenvectors = []
        self.energies = []
        
        # Diagonalize the Hamiltonian
        if (diagonalize):
            self.diagonalize()

    def get_energies(self, sortWithE: bool=True, shiftToZero: bool=True, Eimag_threshold: float=1e-10):
        '''
        Determine energies of the CEF hamiltonian.

        Updates
        -------
        self.energies
        '''

        # Diagonalize the Hamiltonian
        # Other functions than scipy.linalg.schur do not produce perpendicular subspaces for degenerated eigenvalues
        E, _ = schur(self.hamiltonian)
        E = np.diag(E)

        # Check if energies are largely imaginary
        if np.any(np.abs(np.imag(E)/np.real(E))>Eimag_threshold):
            raise ValueError('Final energies are complex! Check `Eimag_threshold` option.')
        
        self.energies = np.real(E) - int(shiftToZero)*min(np.real(E))     # shift to zero level     

        # Sorting
        if sortWithE:
            sortedIndices = self.energies.argsort()
        else:
            sortedIndices = np.arange(self.Jval)
            
        self.energies =  self.energies[sortedIndices]
        
    
    def diagonalize(self, sortWithE: bool=True, shiftToZero: bool=True, Eimag_precision: int=10):
        """
        Diagonalize the Hamiltonian, and change to the sorted eigenvector base. The default sorting is done according to eigenenergies, so that the first vector [1,0,...,0] is the lowest eigenstate, and the last one [0,...,0,1] is the highest one. Changing the base greatly faiclitates further calculations based on the evaluation of matrix elements.
        
        Updates
        -------
        ``Jx``,``Jy``,``Jz``,``J``,``energies``,``eigenvectors``,``degeneracies``
        
        Returns
        -------
            E, U
                Products of diagonalization, H = U E U^T.
            
        Raises
        ------
            ValueError 
                If the calculated eigenenergies are not real.
                
        """
    
        # Diagonalize the Hamiltonian
        # Other functions than scipy.linalg.schur do not produce perpendicular subspaces for degenerated eigenvalues
        E, U = schur(self.hamiltonian)
        eigenvalues = np.diag(E)

        # Check if energies are largely imaginary
        if np.any( np.abs(np.imag(np.around(eigenvalues, Eimag_precision))) > 0):
            raise ValueError(f'Final energies are complex! {eigenvalues} Check `Eimag_precision` option.')
        
        self.energies = np.real(eigenvalues) - int(shiftToZero)*min(np.real(eigenvalues))     # shift to zero level

        # TODO check if U is orthogonal based on comparison with QR decomposition
        self.eigenvectors = U
        
        # Sorting
        if sortWithE:
            sortedIndices = self.energies.argsort()
        else:
            sortedIndices = np.arange(self.Jval)

        self.eigenvectors = self.eigenvectors[:,sortedIndices]
        self.energies =  self.energies[sortedIndices]
        
        # Change the basis of principal operators to the eigenstate basis with specified sorting scheme
        self.Jx = np.dot(np.dot(self.eigenvectors.conj().transpose(), self.Jx), self.eigenvectors)
        self.Jy = np.dot(np.dot(self.eigenvectors.conj().transpose(), self.Jy), self.eigenvectors)
        self.Jz = np.dot(np.dot(self.eigenvectors.conj().transpose(), self.Jz), self.eigenvectors)
       
        self.J = np.array([self.Jx, self.Jy, self.Jz])

        # Calculate magnetic moment
        mu = self.ion.gJ * np.einsum('ijj->ji', self.J)
        self.moment = np.real(mu + mu.conj())/2

        # Degeneracy
        deg_e = []
        levels = 0
        deg_e.append([self.energies[0], 0])
        for x in self.energies:
            if not np.isclose(deg_e[levels][0], x):
                levels+=1
                deg_e.append([x, 1])
            else:
                deg_e[levels][1] += 1

        self.degeneracies = np.array(deg_e)

        return E, U
        
    def __str__(self, precision: float=4):
        """
        Nice printout of calculated parameters
        
        Parameters:
            precision : int, optional
                How many significant digits should be shown for printout of energy and eigenvector coefficients.
        """
        if not len(self.energies):
            return "Undiagonalized Hamiltonian"

        ret = ""
        ret += "Energy levels and corresponding eigenfunctions:\n"
               
        n = 0
        levels = []
        for level,x in enumerate(self.degeneracies):
            level_str = ''
            energy = x[0]
            degeneracy = int(x[1])
            level_str += f'E({level:d}) =\t{energy:.{precision}} meV\t{degeneracy:2d}fold-degenerated\n'
            
            # List degenerated eigenvectors
            for ev in self.eigenvectors.T[n:n+degeneracy]:
                ev_components = []
                for c,ket in zip(ev,self.freeionkets):
                    if np.abs(c) > 1/np.power(10,precision):    # Arbitrary tolerance
                        ev_components.append(f'({c:.{precision}f}){ket}')
                        
                level_str += f'ev_{level}: ' + ' + '.join(ev_components) + '\n'

            levels.append(level_str)
            n += degeneracy
        
        ret += '\n'.join(levels)
        return ret