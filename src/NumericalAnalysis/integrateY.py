from requirements.imports import *
from requirements.numericalMethods import *
from requirements.variables import *

from interpolationY import y_approx_poly, rm

from scipy.integrate import quad
from scipy.optimize import fsolve

# inverted polynomial for use in integration
def inv_poly(x):
    return 1/y_approx_poly(x)

value = quad(inv_poly, 1, rm)[0] # integ. the LHS of equation dr/y(r)=dt

# essentially int_0^x(t) - value
def unit_integral_solve(x):
    return quad(lambda t: t, 0, x)[0] - value

sol = fsolve(unit_integral_solve, 1)

print(f'Finding upper limit of t in solving dr/y(r)=dt:\n t = {sol}')