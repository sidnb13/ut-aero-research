from requirements import *
from NumericalMethods import CurveFit, RootApproximation, Plotter
from approximations import *

# POWER FIT

def Power(x, a, b, c, d):
    return b*x/(a*np.sqrt(-(x-d)/c))

g = np.asarray([0.84,0.007,150,1])

pow0 = CurveFit(delta,phi_0,f=Power)
pow0Coeffs = pow0.getCoeffs()
powmin = CurveFit(delta,phi_min,f=Power,guess=g)
powminCoeffs = powmin.getCoeffs()

powPlots = [
    {'xvar': delta, 'yvar': Power(delta,*pow0Coeffs), 'label': f'Approx for $\\phi_0$ with ${pow0Coeffs[0]}\\exp({pow0Coeffs[1]}(\\delta-{pow0Coeffs[2]}))+{pow0Coeffs[3]}$'},
    {'xvar': delta, 'yvar': Power(delta,*powminCoeffs), 'label': f'Approx for $\\phi_0$ with ${powminCoeffs[0]}\\exp({powminCoeffs[1]}(\\delta-{powminCoeffs[2]}))+{powminCoeffs[3]}$'},
    {'xvar': delta, 'yvar': phi_0, 'label': '$\\phi_0$'},
    {'xvar': delta, 'yvar': phi_min, 'label': '$\\phi_\\mathrm{min}$'}
]
powFile = 'plots/exponential-approx.pdf'

powPlotObj = Plotter(powPlots,limits,powFile,{'x': '$\\delta$', 'y': '$\\phi$'})

powPlotObj.plot(show=True)