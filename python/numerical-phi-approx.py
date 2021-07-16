from initialize import *

class RootApproximation:
    def __init__(func1, phi, betaVal, func2=None, numIter=200):
        self.f1 = func1
        self.f2 = func2
        self.phi = phi
        self.betaVal = betaVal
        self.numIter = numIter
    
    def bisectionRoots():
        try:
            # iterate through until root is found
            a,b = np.amin(self.phi),np.amax(self.phi)
            c,n = 0,1
            # calculating iteration upper bound
            iterMax = int(np.log((b-a)/EPSILON)/np.log(2))
        except:
            return 0

        while n <= iterMax:
            c = (a+b)/2
            if (b-a)/2 <= EPSILON or np.isclose(self.f1(c,self.betaVal),0,atol=EPSILON):
                return c
            if np.sign(self.f1(a,self.betaVal)) == np.sign(self.f1(c,self.betaVal)):
                a = c
            else: 
                b = c
            n = n+1
        return c

    def newtonRoots():
        # guess the midpoint of interval
        x0 = np.median(self.phi)
        # iterate
        for _ in range(self.numIter):
            y = self.f1(x0,self.betaVal)
            yd = self.f2(x0,self.betaVal)

            x1 = x0 - y/yd
            if np.abs(x1 - x0) <= EPSILON:
                return x0
            x0 = x1

class Plotter:
    def __init__(plotList, limDict, fileName):
        self.plotList = plotList
        self.limDict = limDict
        self.fileName = fileName

    def plot(show=True,save=False):
        for plot in self.plotList:
            plt.plot(plot['xvar'],plot['yvar'],label=plot['label'])
        
        plt.xlim(*self.limDict['xlim'])
        plt.ylim(*self.limDict['ylim'])
        plt.legend()
        plt.savefig(self.fileName,backend='pgf',bbox_inches='tight')
        
class CurveFit:
    def __init__(f, delta, phi, k):
        self.delta = delta
        self.phi = phi
        self.k = k
        self.f = f

    def coeffs():
        coeffs, pcov = curve_fit(f, self.delta, self.phi)
        return coeffs

    # generating a closed-form expression by solving Ac=b for a vector c
    def polyApprox():
        m = self.k + 1 # num of coeff. and k is degree
        # define vectors and matrices constant
        b = np.zeros(m)
        Amat = np.zeros((m,m))
        # populate b and A
        for j in range(0,m):
            b[j] = np.sum(np.power(self.delta,self.k-j)*self.phi)
            for i in range(0,m):
                Amat[j,i] = np.sum(np.power(self.delta,2*self.k-i-j))
        # solve equation Ac=b for c using numpy routine
        c = np.linalg.solve(Amat,b)
        
        # generate a lambda function of the curve
        return polyLambdaFromCoeffs(c)
    
    # risky function that converts a polynomial coefficient vector -> lambda
    def polyLambdaFromCoeffs(c):
        polyString = ''
        for i in range(0,np.size(c)):
            polyString = polyString + f'{c[i]}*np.power(d,{np.size(c)-i-1}) + ' # need the -1 for a constant term on end
        polyString = polyString.rstrip(' + ')

        print(polyString)

        return eval(f'lambda d: {polyString}')

# generating vectors and plots for delta/beta
delta = np.linspace(0.001,1,N)
beta = np.sqrt(1/delta)

phi_0 = np.zeros(N)
phi_min = np.zeros(N)

for i in range(N):
    # populate phi_0
    phi_0[i] = RootApproximation(yPhi,phiIntReal(beta[i]),beta[i]).bisectionRoots()
    # populate phi_min
    phi_min[i] = RootApproximation(wPhiDeriv,phiInt,beta[i],func2=wPhiDeriv2Inst,numIter=100).newtonRoots()


#TODO: refactor everything below this to interface with the classes
# control panel to enable/disable plots
ELLIPTIC = True
EXP = False
CHEBYSHEV = False
PLOT_CURVE_FITS = False
PLOT_NUMERICS = False

if ELLIPTIC:
    def elliptic(k, n): # unfortunately cannot be vectorized, do not want to use a power series
        intvals = np.zeros(N)
        for i in range(0,N):
            intvals[i], temp = integrate.quad(lambda x: -np.pi/(2*n) + 1/(n*np.sqrt((1-np.power(x,2))*(1-k[i]**2*x**2))),0,1)
        return intvals
    
    phi0_coeffs, pcov = curve_fit(elliptic, delta, phi_0)
    phimin_coeffs, pcov = curve_fit(elliptic, delta, phi_min)

    plt.plot(delta,elliptic(delta,*phi0_coeffs),label=f'Approx for $\\phi_0$ with $n={phi0_coeffs[0]}$')
    plt.plot(delta,elliptic(delta,*phimin_coeffs),label=f'Approx for $\\phi_\\mathrm{{min}}$ with $n={phimin_coeffs[0]}$')

    plt.plot(delta,phi_0,label='Original $\\phi_0$')
    plt.plot(delta,phi_min,label='Original $\\phi_\\mathrm{min}$')

    plt.xlim(np.amin(delta),np.amax(delta))
    plt.ylim(0,1)

    plt.legend()
    plt.savefig('plots/elliptical-approx.pdf',backend='pgf',bbox_inches='tight')

if EXP:
    # power fit
    def func(x, a, b, c, d):
        return a * np.exp(b*(x-c)) + d

    guess = np.asarray([1,5,1.05,0])

    phi0_coeffs, pcov = curve_fit(func, delta, phi_0,p0=guess)
    phimin_coeffs, pcov = curve_fit(func, delta, phi_min,p0=guess)

    # approximations
    plt.plot(delta,func(delta,*phi0_coeffs),label='Approx $\\phi_0$')
    plt.plot(delta,func(delta,*phimin_coeffs),label='Approx $\\phi_\\mathrm{min}$')

    phi0_texstring = f'$\\phi_0={phi0_coeffs[0]}e^{{{phi0_coeffs[1]}(\\delta-{phi0_coeffs[2]})}}+{phi0_coeffs[3]}$'
    phimin_texstring = f'$\\phi_\\mathrm{{min}}={phimin_coeffs[0]}e^{{{phimin_coeffs[1]}(\\delta-{phimin_coeffs[2]})}}+{phimin_coeffs[3]}$'

    print('%s\n%s' % (phi0_texstring,phimin_texstring))

    # original
    plt.plot(delta,phi_0,label='Original $\\phi_0$')
    plt.plot(delta,phi_min,label='Original $\\phi_\\mathrm{min}$')

    plt.legend()
    plt.savefig('plots/exponential-approx.pdf',backend='pgf',bbox_inches='tight')

    plt.show()

if CHEBYSHEV:
    # most accurate chebyshev approximation, 
    # return chebyshev polynomial coefficient vector
    def chebyshevApproximation(delta, phi, deg):
        chebCoeffs = cheb.chebfit(delta,phi,deg)
        return (cheb.chebval(delta,chebCoeffs,tensor=False), cheb.cheb2poly(chebCoeffs))

    phi0, phi0_c = chebyshevApproximation(delta,phi_0,2)
    phimin, phimin_c = chebyshevApproximation(delta,phi_min,2)

    # chebyshev functions
    # flip because of return order
    phi0_f = polyLambdaFromCoeffs(np.flip(phi0_c))
    phimin_f = polyLambdaFromCoeffs(np.flip(phimin_c))

    plt.plot(delta,phi0_f(delta),label='Chebyshev approximation for $\\phi_0$')
    plt.plot(delta,phimin_f(delta),label='Chebyshev approximation for $\\phi_\\mathrm{min}$')

    plt.plot(delta,phi_0,label='Original $\\phi_0$')
    plt.plot(delta,phi_min,label='Original $\\phi_\\mathrm{min}$')

    plt.xlim(np.amin(delta),np.amax(delta))
    plt.ylim(0,1)

    plt.legend()
    plt.savefig('plots/cheby-approx.pdf',backend='pgf',bbox_inches='tight')

    plt.show()

if PLOT_CURVE_FITS:
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

if PLOT_NUMERICS:
    plt.plot(delta,phi_min, label='$\\phi_\\mathrm{min}$')
    plt.plot(delta,phi_0,label='$\\phi_0$')

    plt.xlim(np.amin(delta),np.amax(delta))
    plt.ylim(0,1)

    plt.xlabel('$\\delta$')
    plt.ylabel('$\\phi$')

    plt.legend(loc='upper left')
    plt.savefig('plots/phi-functions.pdf',backend='pgf',bbox_inches='tight')