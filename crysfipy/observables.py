from . import constants as C
# from .cefmatrices import *
from .cefion import CEFion

import numpy as np

from typing import Tuple
    

def boltzman_population(energies: np.ndarray, temperature: float) -> np.ndarray:
    r'''
    Calculate the population of energy levels at given temperature based on the Boltzmann statistic.
    :math:`p_i = \frac{1}{Z} e^{-E_i/k_B T}`
    :math:`Z = \sum_i e^{-Ei/k_BT}`
    One important distinction, is that this function works with eigenvalues (energies) from the whole Hilbert space,
    as it needs to evaluate :math:`Z` on its own. This works well for the total angular momentum Hilbert space, and does care about degeneracies.
    
    Args:
        energies: List of energy levels in meV units
        temperature: Temperature at which to evaluate the statistic

    Returns:
        List of occupation probabilities for each energy level.
    '''
    
    p = np.exp(-np.array(energies - min(energies))*C.meV2K/temperature)
    Z = sum(p)
    return p / Z


def neutronint(cefion: CEFion, temperature: float, Q: Tuple[float, float, float] , scheme: str, Ei: float=1e+6, B_iso: float=0) -> Tuple[np.ndarray, np.ndarray]:
    r"""
    Returns matrix of energies and inelastic neutron scattering spectral weights for all possible transitions at given temperature.
    
    The spectral weight is calculated by equation from Enderle book following Stephane Raymond article.
    
    | :math:`S(\vec{Q},\omega) = N (\gamma r_0)^2 f^2_m(Q) e^{-W(Q)} \sum_{if} \frac{k_f}{k_i} p_i |<\lambda_f|J_\perp|\lambda_i>|^2 \delta(E_i - E_f - \hbar \omega)`
    
    where:

    | :math:`N (\gamma r_0)^2` : ignored, acts as units.
    | :math:`f^2_m(Q)` : magnetic form factor, taken from internal tables in ``mikibox.crysfipy.Ion`` class.
    | :math:`e^{-W(Q)}` : Debye-Waller factor, where :math:`W(Q) = \frac{1}{2} B Q^2`. :math:`B = B_{iso}` is the isotropic mean-square displacement obtained from crystal structure refinement.
    | :math:`\frac{k_f}{k_i}` : scaling factor calculated from energy, which is used more widely :math:`\frac{k_f}{k_i} = \sqrt{1-\frac{\Delta E}{E_i}}`. there is a minus under the square root, because positive energy transfer corresponds to neutron energy loss.
    | :math:`p_i` : Boltzmann population factor.
    | :math:`|<\lambda_f|J_\perp|\lambda_i>|^2` : matrix elements, exact description depends on ``Q``, see below.

    
    The intensities are evaluated based on the :math:`|<\lambda_f|J_\perp|\lambda_i>|^2` matrix elements, which form a matrix :math:`|J_\perp|^2`.
    Two main cases are implemented and encoded in the :data:`Q` parameter.
    
    Args:
        cefion: :class:`CEFion` Rare-earth ion in the crystal field
        tempreature:  Temperature in *K*
        Q :  List of Q vectors used to evaluate the spectral weight
        
        scheme : 'powder', 'single-crystal':
            Scheme according to which :math:`|J_\perp|^2` is calculated.
            
            * powder :  :math:`|<\lambda_f|J_\perp|\lambda_i>|^2 = 2/3\sum_\alpha |<\lambda_f|J_\alpha|\lambda_i>|^2` (default).
            * (single-crystal : :math:`|<\lambda_f|J_\perp|\lambda_i>|^2 = \sum_\alpha (1-\frac{Q_\alpha}{Q})|<\lambda_f|J_\alpha|\lambda_i>|^2`. Q is a vector representing a direction in the reciprocal space in respect to which a perpendicular projection of :math:`J` will be calculated.

            
    Returns:
        energies: Array containing energies of the transitions
        intensities: Array containing intensities of the transitions
            
    Raises:
        ValueError 
            When an invalid ``Q`` parameter is chosen, or the dimension of the ``Q`` vector is not 3. 
        RuntimeWarning
            When Q=[0,0,0], where the spectral weight is ill defined.        
    """
    
    # The way it is calculated is that the factor within the sum is calculated as a matrix, which is flattened at the end. Energy degeneracies are not taken into account, as it is easier to handle.
    
    # Magnetic form factor
    f2m = cefion.ion.mff(np.linalg.norm(Q))**2
    
    # Debye-Waller factor
    eDW = np.exp(-1/3 * np.linalg.norm(Q)**2 * B_iso)
    
    # Tricky way to create a 2D array of energies associated with transitions between levels
    jumps = cefion.energies - cefion.energies[np.newaxis].T    
    
    # Calculate the |<\Gamma_f|J_\perp|\Gamma_i>|^2 matrix
    if scheme=='single-crystal':
        if Q.shape[0] != 3:
            raise ValueError('Dimension of the `Q` vector is not 3')
            
        # First implementation does not seem to work well
        # Qperp_projectCEFion = ms.perp_matrix(Q)
        # J_perp = np.einsum('ij,jkl',Qperp_projectCEFion, cefion.J)
        # J2_perp = np.einsum('ijk->jk', np.square(np.abs(J_perp)))
        
        J2 = np.square(np.abs(cefion.J)) # J is already diagonalized and in the eigenvector basis
        projection = 1-(Q/np.linalg.norm(Q))**2
        
        J2_perp = np.einsum('i,ijk',projection, J2)
    elif scheme=='powder':
        J2_perp = 2.0 / 3 * np.einsum('ijk->jk',np.square(np.abs(cefion.J)))    
    else:
        raise ValueError('Invalid ``Q`` parameter')
        
        
    # kf/ki factor, which is actually a matrix
    kfki = np.sqrt(1-jumps/Ei)
    
    # Occupation
    transition_probs = boltzman_population(cefion.energies, temperature)
        
    # Multiply the factors, vectors and matrices properly to get spectral weight.
    Sqw = f2m * eDW * kfki * J2_perp * transition_probs[:, np.newaxis]
    
    Denergies = jumps.flatten()
    sorting = Denergies.argsort()
    return (Denergies[sorting], Sqw.flatten()[sorting])


def magnetization(cefion: CEFion, temperature: float, Hfield: Tuple[float,float,float]) -> Tuple[float, float, float]:
    r'''
    Calculate the magnetization of the single ion in the crystal field.
    Returned value is in :math:`\mu_B` units.

    Calculation follows from reevaluating the Hamiltonian with the Zeeman term :math:`\vec{\mu} \cdot \vec{H}`
    and determining the expectation value of magneitzation:

    :math:`M_\alpha = g_J \sum_n p_n |<\lambda_n | \hat{J}_\alpha | \lambda_n>|`

    Parameters:
        cefion: 
            Rare-earth ion in the crystal field
        temperature: 
            Temperature at which to calculate magnetization   
        Hfield: 
            Applied magnetic field. Both direction and value are important.
    '''
    
    # Solve the Hamiltonian of a copy of the given ion
    cefion_inH = CEFion(ion=cefion.ion, cfp=cefion.cfp, Hfield=Hfield, diagonalize=True)

    # The diagonalized Hamiltonians' operators are already transformed into the sorted eigenvector base
    p = boltzman_population(cefion_inH.energies, temperature)
    # M = cefion_inH.ion.gJ * np.abs( np.einsum('ijj,j', cefion_inH.J, p) )
    M = np.einsum('ij,i', cefion_inH.moment, p)
    
    return tuple(M)


def magnetization_exchange(cefion: CEFion, temperature: float, Hfield: Tuple[float,float,float], lam: float, precision: float=1e-5) -> Tuple[float, float, float]:
    r'''
    Calculate the magnetization of the single ion in the crystal field, with correction of exchange/molecular fields.
    Returned value is in :math:`\mu_B` units.

    As for standard magnetization, `observables.magnetization()` calculation follows from reevaluating the Hamiltonian with the Zeeman term,
    but with an effective field :math:`\vec{\mu} \cdot \vec{H}_eff`, where the effective field is corrected for magnetization
    :math:`\vec{H}_eff = \vec{H}_{external} + \vec{H}_exchange = \vec{H}_{external} + \lambda \vec{M}`.
    Exchange field is calulated in a self-consistent manner.

    Parameters:
        cefion: 
            Rare-earth ion in the crystal field.
        temperature: 
            Temperature at which to calculate magnetization.
        Hfield: 
            Applied magnetic field. Both direction and value are important.
        lam:
            Exchange field parameter.
        precision:
            Self-consistent loop will be performed as long as abs(dM/M)>precision.
    '''
    
    # Solve the Hamiltonian of a copy of the given ion
    M = [magnetization(cefion, temperature, Hfield)]
    residual = 42

    while residual > precision:
        Hexchange = lam*np.array(M[-1])
        M_next = magnetization(cefion, temperature, Hfield+Hexchange)
        M.append(M_next)
        dM = [m1-m2 for m1,m2 in zip(M[-1], M[-2])]
        residual = np.linalg.norm(dM)/np.linalg.norm(M[-1])
    
    return M[-1]


def susceptibility(cefion: CEFion, temperatures: np.ndarray, Hfield_direction: Tuple[float,float,float], method: str='perturbation') -> np.ndarray:
    r"""
    Calculate the magnetic susceptibility at listed temperatures, based on one of the implemented methods.
    
    
    `perturbation`:
    Based on assuming an infinitezimal applied field that perturbes the Hamiltonian, the eigenstates
    are calculated by means od perturbation theory and given by formula:

    | :math:`\chi_{CEF} = (g_J \mu_B)^2 \left[ \sum_{n,m \neq n}  p_n \frac{1-exp(-\Delta_{m,n}/k_B T)}{\Delta_{m,n}} |<\lambda_m|J|\lambda_n>|^2  +  \frac{1}{k_B T} \sum_{n} p_n |<\lambda_n|J|\lambda_n>|^2 \right]`
    
    where
    
    | :math:`g_J \mu_B` : Lande factor, Bohr magneton.
    | :math:`p_n` : Ptobabitlity of occupying the energy level :math:`\lambda_n` at given temperature.
    | :math:`\Delta_{m,n}` : Energy of transition between the :math:`\lambda_n` and :math:`\lambda_m` levels.
    | :math:`<\lambda_m|J|\lambda_n>` : J matrix elements.


    

    `magnetization`:
    Based on calculating a numerical derivative of the magnetization, with a very small magnetic field

    | :math:`\chi_{CEF} = \left. \frac{\partial \vec{M}}{\partial \vec{H}} \right|_{\vec{H}=\epsilon}`

    Magnetization is calculated internally by :func:`magnetization`. :math:`\epsilon = 10^{-8}`

  

    
    Parameters:
        cefion : :class:`crysfipy.CEFion`
            Rare-earth ion in crystal field\
        temperatures : ndarray
            Array ocntaining temperatures at which to compute the susceptibility.
        Hfield_direction:
            Direction of the applied magnetic field. Value can be arbitrary, it is normalized in the code.
        method: optional, 'perturbation', 'magnetization'
            Method by which to calculate susceptibility. Old implementation 
            
    Returns:
        List of susceptibility values calculated at given temperatures. In the units of :math:`\mu_B^2`.
    """
    


    if method=='magnetization':
        susceptibility = np.zeros(len(temperatures))
        eps = 1e-5 
        
        for it, temperature in enumerate(temperatures):
            Hfield = eps * np.array(Hfield_direction)/np.linalg.norm(Hfield_direction)
            M = magnetization(cefion, temperature, Hfield)
            susceptibility[it] = np.linalg.norm(M)/eps
    elif method=='perturbation':
        susceptibility = np.empty(len(temperatures))
        #susceptibility.fill(np.nan)

        # Tricky way to create a 2D array of energies associated with transitions between levels
        jumps = cefion.energies - cefion.energies[np.newaxis].T

        # Clean up
        jumps_with_zero_energy = np.where(np.abs(jumps)< C.numerical_zero)
        jumps_with_positive_energy = np.where(jumps>0)
        # jumps_with_negative_energy = np.where(jumps<0) # obsolete but left for clarifty, these will be set to zero

        # Calculate the J^2 matrix
        J2 = np.square(np.abs(cefion.J))
        J2_directed = np.einsum('i,ijk',Hfield_direction, J2)/np.linalg.norm(Hfield_direction)  # TODO WTF is this projection?

        print(jumps)
        for it, temperature in enumerate(temperatures):
            # Define the transition matrix taking into account the diagonal and negative energy jumps
            Tmx = np.zeros(np.shape(jumps))
            Tmx[jumps_with_zero_energy] = 1*C.meV2K/temperature
            Tmx[jumps_with_positive_energy] = -np.expm1(-jumps[jumps_with_positive_energy]*C.meV2K/temperature)/jumps[jumps_with_positive_energy]

            # Include the Boltzmann factor
            transition_probs = boltzman_population(cefion.energies, temperature)
            Tmx = Tmx * transition_probs[:, np.newaxis]

            
            susceptibility[it] = (cefion.ion.gJ)**2 * np.sum(Tmx * J2_directed)
    else:
        raise ValueError('Unknown method to calculate magnetization.')
        
    return susceptibility
        

def thermodynamics(cefion: CEFion, T: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r"""
    Calculate the fundamental thermodynamic values as a function of temperature.
    
    These functions are calculated altogether taking advantage of the fact 
    that thermodynamics can be determined from the partition function :math:`Z` 
    taken as the weighing function and calculating moments of average energy <E^n>.
    
    | Partition function: :math:`Z = \sum_n e^{-\beta E_n}`
    | Average (internal) energy: :math:`\langle E \rangle = U = \frac{1}{Z} \sum_n E_n e^{-\beta E_n}`
    | Free energy: :math:`F = -T ln(Z)`
    | Entropy: :math:`S = (U-F)/T`
    | Specific heat: :math:`C_V = (<E^2> - <E>)/T^2`
    
    Parameters:
        cefion:
            Rare-earth ion in crystal field
        T : ndarray
            Temperature in Kelvin

            
    Returns:
        Z, E, F, S Cv : The partition function, average energy (internal energy), free energy, entropy, and heat capacity, respectively.
    """
    Z = np.zeros(len(T))
    E = np.zeros(len(T))
    E2 = np.zeros(len(T))
    
    for En in cefion.energies:
        print( np.exp(-En*C.meV2K/T) )
        Z += np.exp(-En*C.meV2K/T)
        E += En*np.exp(-En*C.meV2K/T)
        E2 += En*En*np.exp(-En*C.meV2K/T)

    E /= Z
    E2 /= Z
    F = -T*np.log(Z)
    S = (E-F)/T
    Cv = (E2-E**2)/T**2
            
    return Z, E, F, S, Cv
