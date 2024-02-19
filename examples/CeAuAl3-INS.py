import numpy as np
import matplotlib.pyplot as plt

import crysfipy


# Analysis of the crystal field effects in CeAuAl3
#
# Based on:
# https://doi.org/10.1088/1361-648X/aac408

calculateTAS = False
calculateTOF = False
calculateFormFactor = False
calculateMagnetism = True


CEFpars = crysfipy.CEFpars('D4', [1.203, -0.001, 0.244], 'meV')
cefion = crysfipy.CEFion(crysfipy.Ion('Ce'),[0,0,0], CEFpars)

a, c = 4.35, 10.6
Q1 = np.array([0,0,4])*2*np.pi/c
Q2 = np.array([1,1,0])*2*np.pi/a

print(cefion.cfp)
print(cefion)



fig, ax = plt.subplots()

e = np.linspace(-10,40,1000)
temperature = 10

for Q in [Q1, Q2]:
    De, Dint = crysfipy.neutronint(cefion,temperature, Q, scheme='single-crystal', Ei=1.1*np.max(e))
    mainTransitions = (np.where(Dint> 0.1*np.max(Dint)))
    
    spectrum = np.zeros(e.shape)
    for Etr, Itr in zip(De, Dint):
        spectrum += Itr*np.exp(-(Etr-e)**2/0.2)

    label = f'Q=({Q[0]} {Q[1]} {Q[2]})'
    #label = f'T={temperature} K'
    ax.plot(e, spectrum, label=label)

ax.legend()
fig.savefig(r'C:\Users\Stekiel\Documents\GitHub\crysfipy\examples\CeAuAl3-TAS-spectra.png')