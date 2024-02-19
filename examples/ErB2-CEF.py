import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.optimize import curve_fit

from timeit import default_timer as timer

import sys
sys.path.append('C:/Users/Stekiel/Documents/GitHub/crysfipy')
import crysfipy


# Analysis of the crystal field effects in CePdAl3
calculateMagnetism = True

Bij = [1.26987, -0.00062433, -2.86435e-07, 1.01125e-05]
# Bij = [1.26987, 0, 0, 0]

CEFpars = crysfipy.CEFpars('D6h', Bij, 'meV')

cefion_ErB2 = crysfipy.CEFion(crysfipy.Ion('Er'),[0,0,0], CEFpars)
cefion = cefion_ErB2
print(CEFpars)
print(cefion_ErB2.cfp)
print(cefion_ErB2)








if calculateMagnetism:
    temperature=0.1
    temperatures = np.linspace(1,300,601)
    kB = 8.61733e-2 # in meV/K


    fig, axs = plt.subplots(figsize=(5,13), nrows=4, tight_layout=True)

    # Specific heat
    ax = axs[0]

    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Magnetic susceptibility ($\mu_B^2$/ion)')

    _, _, _, _, cp = crysfipy.thermodynamics(cefion, temperatures)
    for ecef in cefion.energies:
        Tcef =ecef/kB
        if Tcef < temperatures.max():
            ax.axvline(Tcef, 0, 0.5, color='red')

    ax.plot(temperatures, cp, marker='d', color='black', label='cp')

    ax.legend(title='Specific heat')


    # Set up the susceptibility calculated from magnetization
    ax = axs[1]

    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Magnetic susceptibility ($\mu_B^2$/ion)')

    chi_100 = crysfipy.susceptibility(cefion, temperatures , (1,0,0), method='magnetization')
    chi_001 = crysfipy.susceptibility(cefion, temperatures , (0,0,1), method='magnetization')

    ax.plot(temperatures, chi_100, marker='^', color='red', label='$\chi$ || [100]', mfc=None)
    ax.plot(temperatures, chi_001, marker='o', color='blue', label='$\chi$ || [001]', mfc=None)

    ax.set_yscale('log')
    ax.legend(title=f'ErB$_2$ magnetization')

    # Set up the susceptibility calculated from perturbation
    ax = axs[2]

    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Magnetic susceptibility ($\mu_B^2$/ion)')
    #ax.set_ylim(0,4)

    chi_100 = crysfipy.susceptibility(cefion, temperatures , (1,0,0), method='perturbation')
    chi_010 = crysfipy.susceptibility(cefion, temperatures , (0,1,0), method='perturbation')
    chi_001 = crysfipy.susceptibility(cefion, temperatures , (0,0,1), method='perturbation')

#    print(chi_100)

    ax.plot(temperatures, chi_100, marker='^', color='red', label='$\chi$ || [100]', mfc=None)
    ax.plot(temperatures, chi_010, marker='s', color='green', label='$\chi$ || [010]', mfc=None)
    ax.plot(temperatures, chi_001, marker='o', color='blue', label='$\chi$ || [001]', mfc=None)

    ax.set_yscale('log')
    ax.legend(title=f'ErB$_2$ perturbation')

    # Set up the magnetization
    ax = axs[3]

    ax.set_xlabel('Applied field (T)')
    ax.set_ylabel('Magnetic moment ($\mu_B$/ion)')

    fields = np.linspace(1e-5, 1000, 101)
    Ma, Mb, Mc = [], [], []
    for field in fields:
        Ma.append(crysfipy.magnetization(cefion, temperature, (field,0,0))[0])
        Mb.append(crysfipy.magnetization(cefion, temperature, (0,field,0))[1])
        Mc.append(crysfipy.magnetization(cefion, temperature, (0,0,field))[2])



    ax.plot(fields, Ma, marker='^', color='red', label='H || [100]')
    ax.plot(fields, Mb, marker='s', color='green', label='H || [010]')
    ax.plot(fields, Mc, marker='o', color='blue',label='H || [001]')

    ax.legend(title=f'ErB$_2$, T={temperature} K')

    fig.savefig(r'C:\Users\Stekiel\Documents\GitHub\crysfipy\examples\ErB2-magnetism-mikibox.png',dpi=200)