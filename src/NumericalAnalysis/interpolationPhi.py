from requirements.imports import *
from requirements.numericalMethods import *
from requirements.variables import *

import numpy.polynomial.chebyshev as chebyshev
from scipy.interpolate import lagrange

NUM_PTS = 7 # deg = NUM_PTS - 1

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
        if (np.abs(nodes[j] - delta[i]) <= EPSILON):
            phi0_pts[j] = phi_0[i]
            phimin_pts[j] = phi_min[i]

# perform lagrange interpolation
poly_phi0 = lagrange(nodes, phi0_pts)
poly_phimin = lagrange(nodes, phimin_pts)

plt.figure(1)

plt.xlim(*limits['xlim'])
plt.ylim(*limits['ylim'])

for i in range(NUM_PTS):
    plt.plot(nodes[i], phi0_pts[i], 'ro')
plt.plot(delta,poly_phi0(delta),label='polynomial')
plt.plot(delta, phi_0, label='$\\phi_0$')

plt.legend()
plt.show()

plt.figure(2)

plt.xlim(*limits['xlim'])
plt.ylim(*limits['ylim'])

for i in range(NUM_PTS):
    plt.plot(nodes[i], phimin_pts[i], 'ro')
plt.plot(delta,poly_phimin(delta),label='polynomial')
plt.plot(delta, phi_min, label='$\\phi_\\mathrm{min}$')

plt.legend()
plt.show()