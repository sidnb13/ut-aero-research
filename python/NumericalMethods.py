from requirements import *

class RootApproximation:
    def __init__(self, func1, phi, betaVal, func2=None, numIter=200):
        self.f1 = func1
        self.f2 = func2
        self.phi = phi
        self.betaVal = betaVal
        self.numIter = numIter
    
    def bisectionRoots(self):
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

    def newtonRoots(self):
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
    def __init__(self, plotList, limDict, fileName, axes):
        self.plotList1 = plotList
        self.limDict = limDict
        self.fileName = fileName
        self.axes = axes

    def plot(self, show=True,save=False):
        for plot in self.plotList1:
            plt.plot(plot['xvar'],plot['yvar'],label=plot['label'])

        plt.xlim(*self.limDict['xlim'])
        plt.ylim(*self.limDict['ylim'])

        plt.ylabel(self.axes['y'])
        plt.xlabel(self.axes['x'])

        plt.legend()

        if show:
            plt.show()
        if save:    
            plt.savefig(self.fileName,backend='pgf',bbox_inches='tight')
        
class CurveFit:
    def __init__(self, delta, phi, k=None, guess=None, f=None):
        self.delta = delta
        self.phi = phi
        self.k = k
        self.f = f
        self.guess = guess

    def getCoeffs(self):
        c, pcov = curve_fit(self.f, self.delta, self.phi, p0=self.guess)
        return c

    # generating a closed-form expression by solving Ac=b for a vector c
    def polyApprox(self):
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

        def polyLambdaFromCoeffs(c):
            polyString = ''
            for i in range(0,np.size(c)):
                polyString = polyString + f'{c[i]}*np.power(d,{np.size(c)-i-1}) + ' # need the -1 for a constant term on end
            polyString = polyString.rstrip(' + ')
            return eval(f'lambda d: {polyString}')
        
        # generate a lambda function of the curve
        return c, polyLambdaFromCoeffs(c)