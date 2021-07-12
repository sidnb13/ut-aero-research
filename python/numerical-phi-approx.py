import numpy as np
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings('ignore')

# set domain and consts
N = 1000
PLOT = False
START = 0
END = np.pi/4

phiInt, step = np.linspace(START,END,N, retstep=True)
# set constants
EPSILON = 1e-6
a_T = 0.2
r_0 = 1
mu = 1
Beta = np.sqrt(4*mu/((3*np.pi+8)*a_T*np.power(r_0,2)))

# define functions
wPhi = lambda phi, Beta : 1/np.cos(2*phi) * (1-(3*phi+2)*4/(np.power(Beta,2)*(3*np.pi+8))) + (np.sin(2*phi)+4) * 2/(np.power(Beta,2)*(3*np.pi+8))

gPhi = lambda phi, Beta : np.power((3*np.pi+8),2)*np.power(Beta,4)*np.power(np.sin(2*phi+np.pi/2),2)\
        -(np.power((3*np.pi+8),2)*np.power(Beta,4)-8*(3*phi+2)*(3*np.pi+8)*np.power(Beta,2)+16*np.power((3*phi+2),2))\
        -2*np.sin(2*phi+np.pi/2)*(4*np.power(np.sin(phi+np.pi/4),2)+6)*((3*np.pi+8)*np.power(Beta,2)-4*(3*phi+2))\
        -np.power(np.sin(2*phi+np.pi/2),2)*(16*np.power(np.sin(phi+np.pi/4),4)+48*np.power(np.sin(phi+np.pi/4),2)+36)
yPhi = lambda phi, Beta : np.sqrt(gPhi(phi,Beta))/(np.power(Beta,2)*(3*np.pi+8)*np.sin(2*phi+np.pi/2))*np.sqrt(mu/r_0)*1/np.tan(phi+np.pi/4)

# derivatives
wPhiDeriv = lambda phi, beta : 2*(np.power(np.sin(2*phi+np.pi/2),2)-3)*np.sin(2*phi+np.pi/2)/((3*np.pi+8)*np.power(beta,2)*np.power(np.sin(2*phi+np.pi/2),3))\
                                - np.cos(2*phi+np.pi/2)*(1-(3*phi+2)*4/((3*np.pi+8)*np.power(beta,2)))/(r_0*np.power(np.sin(2*phi+np.pi/2),3))

wPhiDeriv2Inst = lambda phi, beta: (wPhiDeriv(phi+step,beta)-wPhiDeriv(phi-step,beta))/(2*step)

# calculate either overall or inst. deriv
phiIntReal = lambda beta : phiInt[~np.isnan(yPhi(phiInt,beta))]

yPhiDerivInst = lambda phi, beta : (yPhi(phi+step,beta)-yPhi(phi-step,beta))/(2*step)
yPhiDeriv = lambda beta : np.diff(yPhi(phiIntReal,beta))/np.diff(phiIntReal)

# bisection algorithm to find roots of y(phi) and w'(phi)
def bisectionAlgorithm(f,phiInt,Beta):
    try:
        # iterate through until root is found
        a,b = np.amin(phiInt),np.amax(phiInt)
        c,n = 0,1
        # calculating iteration upper bound
        iterMax = int(np.log((b-a)/EPSILON)/np.log(2))
    except:
        return 0

    while n <= iterMax:
        c = (a+b)/2
        if (b-a)/2 <= EPSILON or np.isclose(f(c,Beta),0,atol=EPSILON):
            return c
        if np.sign(f(a,Beta)) == np.sign(f(c,Beta)):
            a = c
        else: 
            b = c
        n = n+1
    return c

# Newton-Raphson algorithm to find roots
def newtonRaphson(f, f2, phiInt, Beta, numIter):
    # guess the midpoint of interval
    x0 = np.median(phiInt)
    # iterate
    for _ in range(numIter):
        y = f(x0,Beta)
        yd = f2(x0,Beta)

        x1 = x0 - y/yd
        if np.abs(x1 - x0) <= EPSILON:
            return x0
        x0 = x1

# generating vectors and plots for delta/beta
delta = np.linspace(0.001,1,N)
beta = np.sqrt(1/delta)

phi_0 = np.zeros(N)
phi_min = np.zeros(N)

for i in range(N):
    # populate phi_0
    phi_0[i] = bisectionAlgorithm(yPhi,phiIntReal(beta[i]),beta[i])
    # populate phi_min
    phi_min[i] = newtonRaphson(wPhiDeriv,wPhiDeriv2Inst,phiInt,beta[i],100)

# generating a closed-form expression by solving Ac=b for a vector c
def fitCurve(delta,phi,k):
    m = k + 1 # num of coeff. and k is degree
    # define vectors and matrices constant
    b = np.zeros(m)
    Amat = np.zeros((m,m))
    # populate b and A
    for j in range(0,m):
        b[j] = np.sum(np.power(delta,k-j)*phi)
        for i in range(0,m):
            Amat[j,i] = np.sum(np.power(delta,2*k-i-j))
    # solve equation Ac=b for c using numpy routine
    c = np.linalg.solve(Amat,b)
    # generate a lambda function of the curve
    polyString = ''
    for i in range(0,m):
        polyString = polyString + f'{c[i]}*np.power(d,{k-i}) + '
    polyString = polyString.rstrip(' + ')

    return (eval(f'lambda d: {polyString}'),polyString)

func1,f1 = fitCurve(delta,phi_0,7)
func2,f2 = fitCurve(delta,phi_min,40)

plt.plot(delta,func1(delta))
plt.plot(delta,func2(delta))
plt.xlim(np.amin(delta),np.amax(delta))
plt.ylim(0,1)
plt.show()

if PLOT:
    plt.rcParams.update({
        "text.usetex": True,
    })

    plt.plot(delta,phi_min, label='$\\phi_\\mathrm{min}$')
    plt.plot(delta,phi_0,label='$\\phi_0$')

    plt.xlim(np.amin(delta),np.amax(delta))
    plt.ylim(0,1)

    plt.legend(loc='upper left')
    plt.savefig('plots/phi-functions.pdf',backend='pgf',bbox_inches='tight')