from requirements.imports import *
from requirements.numericalMethods import *
from requirements.variables import *

import numpy.polynomial.chebyshev as cheb


# most accurate chebyshev approximation, 
# return chebyshev polynomial coefficient vector
def chebyshevApproximation(delta, phi, deg):
    chebCoeffs = cheb.chebfit(delta,phi,deg)
    return (cheb.chebval(delta,chebCoeffs,tensor=False), cheb.cheb2poly(chebCoeffs))

cheb0, cheb0Coeffs = chebyshevApproximation(delta, phi_0,5)
chebmin, chebminCoeffs = chebyshevApproximation(delta, phi_min,5)

chebPlots = [
    {'xvar': delta, 'yvar': cheb0, 'label': f'Approx for $\\phi_0$, 1st deg Chebyshev'},
    {'xvar': delta, 'yvar': chebmin, 'label': f'Approx for $\\phi_0$, 1st deg Chebyshev'},
    {'xvar': delta, 'yvar': phi_0, 'label': '$\\phi_0$'},
    {'xvar': delta, 'yvar': phi_min, 'label': '$\\phi_\\mathrm{min}$'}
]
chebFile = 'plots/cheby-approx.pdf'

chebPlotObj = Plotter(chebPlots, limits, chebFile, {'x': '$\\delta$', 'y': '$\\phi$'})