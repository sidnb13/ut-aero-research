import requirements
from NumericalMethods import CurveFit, RootApproximation, Plotter
from approximations import *

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