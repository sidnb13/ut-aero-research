from requirements import *
from NumericalMethods import Plotter
from approximations import *

from scipy.interpolate import UnivariateSpline

SMOOTH_FACTOR = 0.001

splinePhi0 = UnivariateSpline(delta, phi_0, s=SMOOTH_FACTOR)
splinePhiMin = UnivariateSpline(delta,phi_min,s=SMOOTH_FACTOR)

plt1List = [
    {'xvar': delta, 'yvar': splinePhi0(delta), 'label': 'Spline $\\phi_0(\\delta)$'},
    {'xvar': delta, 'yvar': splinePhiMin(delta), 'label': 'Spline $\\phi_\\mathrm{min}(\\delta)$'},
    # {'xvar': delta, 'yvar': phi_0, 'label': '$\\phi_0$'},
    # {'xvar': delta, 'yvar': phi_min, 'label': '$\\phi_\\mathrm{min}$'}
]

plt.subplot(211)

plt.title(f'Spline approximation with $s={SMOOTH_FACTOR}$')

for p in plt1List:
    plt.plot(p['xvar'], p['yvar'], label=p['label'])

plt.ylabel('$\\phi$')

plt.xlim(*limits['xlim'])
plt.ylim(*limits['ylim'])

plt.legend()

plt.subplot(212)

plt.plot(delta, np.abs(phi_0-splinePhi0(delta)), label='err for spl($\\phi_0$)')
plt.plot(delta, np.abs(phi_min-splinePhiMin(delta)), label='err for spl($\\phi_\\mathrm{min}$)')

plt.ylabel('error')
plt.xlabel('$\\delta$')

plt.xlim(*limits['xlim'])
plt.ylim(0,0.01)

plt.legend()

print('---------- phi_0 spline approx ----------')
print(f'Error = \t{splinePhi0.get_residual()}')
print(f'Knots = \t{splinePhi0.get_knots()}')
print(f'Coeffs = \t{splinePhi0.get_coeffs()}\n')

print('---------- phi_min spline approx ----------')
print(f'Error = \t{splinePhiMin.get_residual()}')
print(f'Knots = \t{splinePhiMin.get_knots()}')
print(f'Coeffs = \t{splinePhiMin.get_coeffs()}')

plt.show()
plt.savefig('plots/spline.pdf',backend='pgf',bbox_inches='tight')