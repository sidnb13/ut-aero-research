from requirements import *
from NumericalMethods import CurveFit, RootApproximation, Plotter

# output data in desmos-readable format

with open('phi0.txt','w') as f:
    for i in range(N):
        f.write(f'{delta[i]}\t{phi_0[i]}\n')

with open('phimin.txt','w') as f:
    for i in range(N):
        f.write(f'{delta[i]}\t{phi_min[i]}\n')