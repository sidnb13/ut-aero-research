from matplotlib import lines
import numpy as np
import matplotlib.pyplot as plt

# set domain
ARGLEN = 1000
START = 0
END = np.pi/4

phiInt = np.linspace(START,END,ARGLEN)
# set constants
EPSILON = 1e-6
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
wPhiDeriv = lambda phi, beta : (a_T*np.power(r_0,2)*(np.power(np.sin(2*phi+np.pi/2),3)-3)*np.sin(2*phi+np.pi/2)\
            -2*mu*np.cos(2*phi+np.pi/2)*(1-a_T*np.power(r_0,2)*(3*phi+2)/mu))/(2*mu*r_0*np.power(np.sin(2*phi+np.pi/2),3))

# bisection algorithm to find roots of y(phi) and w'(phi)
def bisectionAlgorithm(f,phiInt):
    # calculating iteration upper bound
    iterMax = int(np.log2(EPSILON_0/EPSILON))
    # iterate through until root is found
    a,b = phiInt[np.argmin(phiInt)],phiInt[np.argmax(phiInt)]
    c,n = 0,0

    while n <= iterMax:
        c = (a+b)/2
        if (b-a)/2 <= EPSILON or f(c,Beta) == 0:
            return c
        if np.sign(f(a,Beta)) == np.sign(f(c,Beta)):
            a = c
        else: b = c
        n = n+1
    return None

# Newton-Raphson algorithm to find roots
def NewtonRaphson():
    pass

# root = bisectionAlgorithm(wPhiDeriv,phiInt)
# print(root)

# plt.plot(phiInt,wPhiDeriv(phiInt,Beta))
# plt.plot(phiInt,np.zeros(ARGLEN),linestyle='dashdot')
# plt.plot(root,wPhiDeriv(root,Beta),marker='o', markersize=5, color="black")

# plt.xlim(START,END)
# plt.ylim(-1,5)
# plt.show()