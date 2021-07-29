from requirements.imports import *
from requirements.numericalMethods import *
from requirements.variables import *

import numpy.polynomial.chebyshev as chebyshev
from scipy.interpolate import lagrange

NUM_PTS = 9 # deg = NUM_PTS - 1

cheb_n = 2*NUM_PTS# due to even symmetry
cheb_nodes = chebyshev.chebpts1(cheb_n)
nodes = cheb_nodes[cheb_nodes >= 0] # cannot use x-axis < 0

# find correspondings points from both phi_0,phi_min, y(r)
phi0_pts = np.zeros(NUM_PTS)
phimin_pts = np.zeros(NUM_PTS)
y_pts = np.zeros(NUM_PTS)

# collect corresponding interpolation points
for i in range(N):
    for j in range(NUM_PTS):
        if np.abs(delta[i] - nodes[j]) <= 1e-3:
            phi0_pts[j] = phi_0[i]
            phimin_pts[j] = phi_min[i]

# perform lagrange interpolation
poly_phi0 = lagrange(nodes, phi0_pts)
poly_phimin = lagrange(nodes, phimin_pts)

def latexify(c,var):
    s = f'$f({var})='
    for i in range(c.shape[0]):
        ind_var = f'{var}^{{{c.shape[0]-i-1}}}' if i < c.shape[0]-1 else ''
        s += f'{int(np.round(c[i],0))}{ind_var}+'
    return s.strip('+').replace('+-','-') + '$'

# plt.figure(1)

# plt.xlim(*limits['xlim'])
# plt.ylim(*limits['ylim'])

# for i in range(NUM_PTS):
#     plt.plot(nodes[i], phi0_pts[i], 'ro')
# plt.plot(delta,poly_phi0(delta),label='polynomial')
# plt.plot(delta, phi_0, label='$\\phi_0$')

# plt.xlabel('$\\delta$')

# plt.legend()
# plt.show()

# plt.figure(2)

# plt.xlim(*limits['xlim'])
# plt.ylim(*limits['ylim'])

# for i in range(NUM_PTS):
#     plt.plot(nodes[i], phimin_pts[i], 'ro')
# plt.plot(delta,poly_phimin(delta),label='polynomial')
# plt.plot(delta, phi_min, label='$\\phi_\\mathrm{min}$')

# plt.xlabel('$\\delta$')

# plt.legend()
# plt.show()