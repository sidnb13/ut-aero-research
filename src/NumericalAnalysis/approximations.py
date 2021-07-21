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

# regplot

plist = [
    {'xvar': delta, 'yvar': phi_0, 'label': '$\\phi_0$'},
    {'xvar': delta, 'yvar': phi_min, 'label': '$\\phi_\\mathrm{min}$'}
]

regular_plot = Plotter(plist,limits,'plots/phi-functions.pdf', {'x': '$\\delta$', 'y': '$\\phi$'})
