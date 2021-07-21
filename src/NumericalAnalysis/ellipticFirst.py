from requirements import *
from NumericalMethods import CurveFit, RootApproximation, Plotter
from approximations import *

# ELLIPTIC APPROXIMATION

def elliptic(k, n): # unfortunately cannot be vectorized, do not want to use a power series
    intvals = np.zeros(N)
    for i in range(0,N):
        intvals[i], temp = integrate.quad(lambda x: -np.pi/(2*n) + 1/(n*np.sqrt((1-np.power(x,2))*(1-k[i]**2*x**2))),0,1)
    return intvals

ell0 = CurveFit(delta,phi_min,delta,f=elliptic)
ell0Coeffs = ell0.getCoeffs()
ellmin = CurveFit(delta,phi_0,delta,f=elliptic)
ellminCoeffs = ellmin.getCoeffs()

ellPlots = [
    {'xvar': delta, 'yvar': elliptic(delta,*ell0Coeffs), 'label': f'Approx for $\\phi_0$ with $n={ell0Coeffs[0]}$'},
    {'xvar': delta, 'yvar': elliptic(delta,*ellminCoeffs), 'label': f'Approx for $\\phi_\\mathrm{{min}}$ with $n={ellminCoeffs[0]}$'},
    {'xvar': delta, 'yvar': phi_0, 'label': '$\\phi_0$'},
    {'xvar': delta, 'yvar': phi_min, 'label': '$\\phi_\\mathrm{min}$'}]
ellFile = 'plots/elliptical-approx.pdf'

ellPlotObj = Plotter(ellPlots,limits,ellFile,{'x': '$\\delta$', 'y': '$\\phi$'})