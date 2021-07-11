import numpy as np
import matplotlib.pyplot as plt

# set domain and consts
N = 1000
PLOT = False
START = 0
END = np.pi/4

phiInt, step = np.linspace(START,END,N, retstep=True)
# set constants
EPSILON = 1e-4
EPSILON_0 = np.argmax(phiInt) - np.argmin(phiInt)
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
    # calculating iteration upper bound
    iterMax = int(np.log(EPSILON_0/EPSILON)/np.log(2))
    # iterate through until root is found
    a,b = phiInt[np.argmin(phiInt)],phiInt[np.argmax(phiInt)]
    c,n = 0,1

    while n <= iterMax:
        c = (a+b)/2
        if (b-a)/2 <= EPSILON or np.isclose(f(c,Beta),0,atol=EPSILON):
            return c
        if np.sign(f(a,Beta)) == np.sign(f(c,Beta)):
            a = c
        else: 
            b = c
        n = n+1
    return None

# Newton-Raphson algorithm to find roots
def newtonRaphson(f, f2, phiInt, Beta, numIter):
    # guess the midpoint of interval
    x0 = f(np.median(phiInt),Beta)
    tol = EPSILON
    epsilon = np.power(EPSILON,2)

    for _ in range(numIter):
        y = f(x0,Beta)
        yd = f2(x0,Beta)

        if np.abs(yd) <= epsilon:
            return x0
        x1 = x0 - y/yd
        if np.abs(x1 - x0) <= tol:
            return x0
        x0 = x1
    return None

# generating vectors and plots for delta/beta
delta = np.linspace(0.001,1,N)
beta = np.sqrt(1/delta)



if PLOT:
    plt.rcParams.update({
        "text.usetex": True,
    })

    plt.xlim(START,END)
    plt.ylim(-1,5)

    plt.legend(loc='upper left')
    plt.savefig('bruh.pdf',backend='pgf',bbox_inches='tight')