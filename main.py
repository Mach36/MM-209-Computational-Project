import numpy as np
import matplotlib.pyplot as plt
from math import *
# from numpy.core.function_base import linspace

T = float(input('Temperature(K) = '))
R = 8.314
### Literature
# delG = delH - T*delS   J/mol
# Ni + 0.5O2 = NiO                 K1 -> AB
# NiS + 2O2 = NiSO4                K2 -> DG
# NiO + SO2 + 0.5O2 = NiSO4        K3 -> DH
# Ni3S2 + 3.5O2 = 3NiO + 2SO2      K4 -> BC
# Ni3S2 + SO2 = 3NiS + O2          K5 -> CF
# NiS + 1.5O2 = NiO + SO2          K6 -> CD
# 3Ni + 2SO2 = Ni3S2 + 2O2         K7 -> EB

# In order to calculate delG of the above reactions, we will use the equation:
# delG_equation = Summation ui*delG_fi(product) - Summation vi*delG_fi(reactants)

delG1 = -239743.2-37.99072*T-0.5*(-T*205.028552)-(-T*29.87376)
delG2 = -872907.92-97.07*T-(-82006.4-T*52.96944)-2*(-T*205.028552)
delG3 = -872907.92-97.07*T-(-239743.2-37.99072*T)-(-296829.696-T*248.1112)-0.5*(-T*205.028552)
delG4 = 3*(-239743.2-37.99072*T)+2*(-296829.696-T*248.1112)-(-202924-T*133.888)-3.5*(-T*205.028552)
delG5 = 3*(-82006.4-T*52.96944)+(-T*205.028552)-(-202924-T*133.888)-(-296829.696-T*248.1112)
delG6 = -239743.2-37.99072*T-296829.696-T*248.1112-(-82006.4-T*52.96944)-1.5*(-T*205.028552)
delG7 = -202924-T*133.888+2*(-T*205.028552)-3*(-T*29.87376)-2*(-296829.696-T*248.1112)

logK1 = log10(exp(-delG1/(R*T)))
logK2 = log10(exp(-delG2/(R*T)))
logK3 = log10(exp(-delG3/(R*T)))
logK4 = log10(exp(-delG4/(R*T)))
logK5 = log10(exp(-delG5/(R*T)))
logK6 = log10(exp(-delG6/(R*T)))
logK7 = log10(exp(-delG7/(R*T)))

### Graph
# AB
plt.vlines(x = -2*logK1, ymin= -150, ymax = -3.5*logK1 + 0.5*logK4)
# DG
plt.vlines(x=-0.5*logK2, ymin=0.25*logK2 - logK3, ymax=200)
# DH
logpo2_3 = np.linspace(-0.5*logK2, 100, 500)
logpso2_3 = -0.5*logpo2_3 - logK3
plt.plot(logpo2_3, logpso2_3)
# BD
logpo2_4 = np.linspace(-2*logK1, -0.5*logK2, 500)
logpso2_4 = 1.75*logpo2_4 + 0.5*logK4
plt.plot(logpo2_4, logpso2_4)
# CF
logpo2_5 = np.linspace(-150, (-logK5 -0.5*logK4)/0.75, 500)
logpso2_5 = logpo2_5 - logK5
plt.plot(logpo2_5, logpso2_5)
# EB
logpo2_6 = np.linspace(-150, -2*logK1, 500)
logpso2_6 = logpo2_6 - 0.5*logK7
plt.plot(logpo2_6, logpso2_6)

plt.xlabel('log(p$_O$$_2$)')
plt.ylabel('log(p$_S$$_O$$_2$)')
plt.title('Predominance Diagram of Ni-S-O system')
plt.ylim((-150,100))

### Labelling of points in the graph:
plt.annotate('A', (-2*logK1, -150))
plt.annotate('B', (-2*logK1, -3.5*logK1 + 0.5*logK4), textcoords="offset points",xytext=(1,-5))
plt.annotate('C', ((-logK5 -0.5*logK4)/0.75, (-logK5 -0.5*logK4)/0.75 - logK5), textcoords="offset points",xytext=(-10, -2))
plt.annotate('D', (-0.5*logK2, 0.25*logK2 - logK3))
plt.annotate('E', (-150, -150 - 0.5*logK7), textcoords="offset points",xytext=(-2,-9))
plt.annotate('F', (-150, -150 - logK5), textcoords="offset points",xytext=(-3,3))
plt.annotate('G', (-0.5*logK2, 100), textcoords="offset points",xytext=(0,-9))
plt.annotate('H', (100, -50 - logK3))
### Labelling of areas in the graph:
plt.text(-logK1 - 75, -150 - 0.5*logK7, 'Ni')
plt.text(-logK1 + 50, -50 - logK3, 'NiO')
plt.text(-0.25*logK2 - 75, 0.125*logK2 - 0.5*logK3 + 50, 'NiS')
plt.text(-0.25*logK2 + 50, 0.125*logK2 - 0.5*logK3 + 50, 'NiSO$_4$')
plt.text(-150, - 0.25*logK7 - 0.5*logK5 - 150, 'Ni$_3$S$_2$', rotation = 45)

plt.show()

