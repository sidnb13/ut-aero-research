from requirements.imports import *
from requirements.numericalMethods import *
from requirements.variables import *

# regplot
plist = [
    {'xvar': delta, 'yvar': phi_0, 'label': '$\\phi_0$'},
    {'xvar': delta, 'yvar': phi_min, 'label': '$\\phi_\\mathrm{min}$'}
]

regular_plot = Plotter(plist,limits,'plots/phi-functions.pdf', {'x': '$\\delta$', 'y': '$\\phi$'})