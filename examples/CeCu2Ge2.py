import numpy as np
import matplotlib.pyplot as plt
import time

import crysfipy

calculateTAS = True
calculateMagnetism = True

# Following calculation is to cross check the calculations for the orthorhombic case with Ce
# Source: Klicpera et al, PHYSICAL REVIEW B 95, 085107 (2017)
#cefion_CePd2Ga2 = ms.crysfipy.CEFion(ms.crysfipy.Ion('Ce'),[0,0,0], ["o", 0.33, 0.472, -0.009, 0.111, 0.055])

# Published results:  E1=7.2 meV,   E2=12.2 meV
# Calculated results: E1=7.091 meV, E2=11.772 meV



CEFpars = crysfipy.CEFpars('D4', [-10.26, -0.056, 2.67], 'K')
# CEFpars = crysfipy.CEFpars('C2', [0.33, 0.472, -0.009, 0.111, 0.055], 'meV')
# CEFpars = crysfipy.CEFpars('D6h', [1.26987, -0.00062433, -2.86435e-07, 1.01125e-05], 'meV')
cefion = crysfipy.CEFion(crysfipy.Ion('Pr'), CEFpars, diagonalize=True)

temperature=0.1
temperatures = np.linspace(0.1,300,100)


if True:
    fig, axs = plt.subplots(figsize=(5,10), nrows=3, tight_layout=True)


    # Set up the susceptibility calculated from magnetization
    ax = axs[0]

    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Magnetic susceptibility ($\mu_B^2$/Ce)')
    ax.set_yscale('log')


    sc = 17.5

    chi_100 = sc*crysfipy.susceptibility(cefion, temperatures , [1,0,0], method='magnetization')
    chi_110 = sc*crysfipy.susceptibility(cefion, temperatures , [1,1,0], method='magnetization')
    chi_001 = sc*crysfipy.susceptibility(cefion, temperatures , [0,0,1], method='magnetization')

    ax.plot(temperatures, chi_100, marker='^', color='red', label='$\chi$ || [100]', mfc=None)
    ax.plot(temperatures, chi_110, marker='s', color='green', label='$\chi$ || [110]', mfc=None)
    ax.plot(temperatures, chi_001, marker='o', color='blue', label='$\chi$ || [001]', mfc=None)

    ax.legend(title=f'CeCu$_2$Ge$_2$ magnetization')

    # Set up the susceptibility calculated from perturbation
    ax = axs[1]

    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Magnetic susceptibility ($\mu_B^2$/Ce)')
    ax.set_yscale('log')
    #ax.set_ylim(0,4)

    chi_100 = crysfipy.susceptibility(cefion, temperatures , [1,0,0], method='perturbation')
    chi_110 = crysfipy.susceptibility(cefion, temperatures , [1,1,0], method='perturbation')
    chi_001 = crysfipy.susceptibility(cefion, temperatures , [0,0,1], method='perturbation')

    #    print(chi_100)

    ax.plot(temperatures, chi_100, marker='^', color='red', label='$\chi$ || [100]', mfc=None)
    ax.plot(temperatures, chi_110, marker='s', color='green', label='$\chi$ || [110]', mfc=None)
    ax.plot(temperatures, chi_001, marker='o', color='blue', label='$\chi$ || [001]', mfc=None)

    ax.legend(title=f'CeCu$_2$Ge$_2$ perturbation')

    # Set up the magnetization
    ax = axs[2]

    ax.set_xlabel('Applied field (T)')
    ax.set_ylabel('Magnetic moment ($\mu_B$/Ce)')

    fields = np.linspace(1e-5, 5, 101)
    Ma, Mab, Mc = [], [], []
    for field in fields:
        Ma.append(crysfipy.magnetization(cefion, temperature, [field,0,0])[0])
        Mab.append(crysfipy.magnetization(cefion, temperature, [field/np.sqrt(2),field/np.sqrt(2),0])[0])
        Mc.append(crysfipy.magnetization(cefion, temperature, [0,0,field])[2])


    ax.plot(fields, Ma, marker='^', color='red', label='H || [100]')
    ax.plot(fields, Mab, marker='s', color='green', label='H || [110]')
    ax.plot(fields, Mc, marker='o', color='blue',label='H || [001]')

    ax.legend(title=f'CeCu$_2$Ge$_2$, T={temperature} K')

    fig.savefig(r'C:\Users\Stekiel\Documents\GitHub\crysfipy\examples\CeCu2Ge2-magnetism-mikibox.png',dpi=200)
