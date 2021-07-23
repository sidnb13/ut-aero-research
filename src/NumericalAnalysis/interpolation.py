from requirements.imports import *
from requirements.numericalMethods import *
from requirements.variables import *

import numpy.polynomial.chebyshev as chebyshev

num_nodes = 5

# get n chebyshev nodes of 1st kind
nodes = chebyshev.chebpts1(num_nodes)

