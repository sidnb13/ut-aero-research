from matplotlib.pyplot import plot
from requirements import *
from NumericalMethods import CurveFit, RootApproximation, Plotter

phi_0 = np.zeros(N)
phi_min = np.zeros(N)

for i in range(N):
    # populate phi_0
    phi_0[i] = RootApproximation(yPhi,phiReal(beta[i]),beta[i]).bisectionRoots()
    # populate phi_min
    phi_min[i] = RootApproximation(wPhiDeriv,phi,beta[i],func2=wPhiDeriv2Inst,numIter=100).newtonRoots()


# common plot limits
limits = {'xlim': (np.amin(delta),np.amax(delta)), 'ylim': (0,1)}

#TODO: refactor everything below this to interface with the classes
# control panel to enable/disable plots

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

# POWER FIT

def Power(x, a, b, c, d):
    return a * np.exp(b*(x-c)) + d

guess = np.asarray([1,5,1.05,0])

pow0 = CurveFit(delta,phi_0,f=Power)
pow0Coeffs = pow0.getCoeffs()
powmin = CurveFit(delta,phi_min,f=Power)
powminCoeffs = powmin.getCoeffs()

powPlots = [
    {'xvar': delta, 'yvar': Power(delta,*pow0Coeffs), 'label': f'Approx for $\\phi_0$ with ${pow0Coeffs[0]}\\exp({pow0Coeffs[1]}(\\delta-{pow0Coeffs[2]}))+{pow0Coeffs[3]}$'},
    {'xvar': delta, 'yvar': Power(delta,*powminCoeffs), 'label': f'Approx for $\\phi_0$ with ${powminCoeffs[0]}\\exp({powminCoeffs[1]}(\\delta-{powminCoeffs[2]}))+{powminCoeffs[3]}$'},
    {'xvar': delta, 'yvar': phi_0, 'label': '$\\phi_0$'},
    {'xvar': delta, 'yvar': phi_min, 'label': '$\\phi_\\mathrm{min}$'}
]
powFile = 'plots/exponential-approx.pdf'

powPlotObj = Plotter(powPlots,limits,powFile,{'x': '$\\delta$', 'y': '$\\phi$'})

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

# normal polynomials

if 0:
    func1 = fitCurve(delta,phi_0,5)
    func2 = fitCurve(delta,phi_min,5)

    plt.subplot(211)

    plt.plot(delta,func1(delta),label='$\\phi_\\mathrm{min}(\\delta);\mathrm{deg}=5$')
    plt.plot(delta,func2(delta),label='$\\phi_0(\\delta);\mathrm{deg}=5$')

    plt.plot(delta,phi_min, label='$\\phi_\\mathrm{min}$')
    plt.plot(delta,phi_0,label='$\\phi_0$')

    plt.legend(loc='upper left')
    plt.xlim(np.amin(delta),np.amax(delta))
    plt.ylim(0,1)

    plt.ylabel('$\\phi$')

    plt.subplot(212)

    plt.plot(delta,np.abs(func1(delta)-phi_0),label='$\\mathrm{err}(\phi_\\mathrm{min}(\\delta))$')
    plt.plot(delta,np.abs(func2(delta)-phi_min), label='$\\mathrm{err}(\phi_0(\\delta))$')

    plt.legend(loc='upper left')
    plt.ylim(0,0.1)
    plt.xlim(np.amin(delta),np.amax(delta))

    plt.xlabel('$\\delta$')
    plt.ylabel('approximation error')

    plt.savefig('plots/phi-curve-fit-poly.pdf',backend='pgf',bbox_inches='tight')

# regplot

plist = [
    {'xvar': delta, 'yvar': phi_0, 'label': '$\\phi_0$'},
    {'xvar': delta, 'yvar': phi_min, 'label': '$\\phi_\\mathrm{min}$'}
]

regular_plot = Plotter(plist,limits,'plots/phi-functions.pdf', {'x': '$\\delta$', 'y': '$\\phi$'})
