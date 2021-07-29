from requirements.imports import *
from requirements.numericalMethods import *
from requirements.variables import *

from interpolationPhi import latexify

import numpy.polynomial.chebyshev as chebyshev
from scipy.interpolate import lagrange

# define r range

rm = 0

wr = lambda r : r_0/np.sqrt(r*(2*r_0-r)) * (1 - (a_T*r_0**2/mu)*(3*np.arcsin(np.sqrt(r/(2*r_0)))-3*np.pi/4+2))+(a_T*r_0/(2*mu))*(r+3*r_0)

yr = lambda r : np.sqrt(mu*(2*r_0-r)/(r*r_0)) * np.sqrt(1-np.power(wr(r),2))

for i in np.linspace(1.5,1.8,N):
    if np.isclose(wr(i), 1, EPSILON):
        rm = i

r,rstep = np.linspace(1,rm,retstep=True)

def get_cheby_nodes(n):
    for i in range(1,100):
        temp = chebyshev.chebpts1(i)
        cheb_nodes = temp[temp >= 0] + np.amin(r) # min guaranteed to be 1, axis translation
        cheb_nodes = cheb_nodes[cheb_nodes <= np.amax(r)] # restrict to domain but retain num pts

        if len(cheb_nodes) is n: 
            return cheb_nodes

rNode = get_cheby_nodes(9) # optimum is 9 nodes
yNode = yr(rNode)

y_approx_poly = lagrange(rNode, yNode)

print(latexify(y_approx_poly.coeffs, 'r'))

# for rN,yN in zip(rNode,yNode):
#     plt.plot(rN,yN,'ro')

# plt.plot(r,yr(r),label='$y_r(r)$')
# plt.plot(r,y_approx_poly(r),label='Lagrange polynomial')

# plt.xlabel('$r$')
# plt.ylabel('$y_r(r)$')

# plt.legend()
# plt.show()

# print(y_approx_poly)