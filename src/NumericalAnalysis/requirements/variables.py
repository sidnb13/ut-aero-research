import numpy as np
import warnings

warnings.filterwarnings('ignore')

# set domain and consts
N = 1000
START = 0
END = np.pi/4

phi, step = np.linspace(START,END,N, retstep=True)
# set constants
EPSILON = 1e-6
a_T = 0.2
r_0 = 1
mu = 1
betaVal = np.sqrt(4*mu/((3*np.pi+8)*a_T*np.power(r_0,2)))

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
phiReal = lambda beta : phi[~np.isnan(yPhi(phi,beta))]

yPhiDerivInst = lambda phi, beta : (yPhi(phi+step,beta)-yPhi(phi-step,beta))/(2*step)
yPhiDeriv = lambda beta : np.diff(yPhi(phiIntReal,beta))/np.diff(phiIntReal)

# generating vectors and plots for delta/beta
delta = np.linspace(0.001,1,N)
beta = np.sqrt(1/delta)