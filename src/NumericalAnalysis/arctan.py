from requirements.imports import *
from requirements.numericalMethods import *
from requirements.variables import *

# ARCTAN FIT

def trig(x, a, b, c, d):
    return a*np.arctan(b*(x-c))+d

trig0 = CurveFit(delta,phi_0,f=trig)
trig0Coeffs = trig0.getCoeffs()
trigmin = CurveFit(delta,phi_min,f=trig)
trigminCoeffs = trigmin.getCoeffs()

trigPlots = [
    {'xvar': delta, 'yvar': trig(delta,*trig0Coeffs), 'label': f'Approx for $\\phi_0$'},
    {'xvar': delta, 'yvar': trig(delta,*trigminCoeffs), 'label': f'Approx for $\\phi_0$'},
    {'xvar': delta, 'yvar': phi_0, 'label': '$\\phi_0$'},
    {'xvar': delta, 'yvar': phi_min, 'label': '$\\phi_\\mathrm{min}$'}
]
trigFile = 'plots/trig-approx.pdf'

trigPlotObj = Plotter(trigPlots,limits,trigFile,{'x': '$\\delta$', 'y': '$\\phi$'})

trigPlotObj.plot(show=True)