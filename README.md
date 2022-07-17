# UT Aerospace Engineering Internship repository

My work is primarily focused on applying numerical methods to closed-form equations involving spaceflight parameters. 
This repository is composed of a writeup and its associated code and static files.

**[See project writeup](./output/writeup.pdf)**

## Tree

```
.
├── output          CPDF of writeup, compiled frequently
├── plots           Plots included in the writeup, generated using matplotlib or GNU Octave
├── presentation    Initial writeup TeX code and PDF
├── src             Python module (NumericalAnalysis) + some Octave codes
├── static          Static PDFs
├── tex             TeX source files for the writeup
└── writeup.tex     Master TeX source file for writeup
```

## Dependencies

You will need the following to play around with the code and writeup sources:

- GNU Octave
- Python 3 and the following modules
  - `numpy`
  - `scipy`
  - `matplotlib`
- A LaTeX distribution

## Numerical methods used

Below are some resources and methods I found useful.

### Approximation methodologies

- Curve fitting: https://www.dam.brown.edu/people/alcyew/handouts/leastsq.pdf
- Bisection method: https://en.wikipedia.org/wiki/Bisection_method
- Newton's method: https://en.wikipedia.org/wiki/Newton%27s_method
- Lagrange interpolation: https://www.neelsomani.com/blog/chebyshev-points-for-lagrange-interpolation.php

### Interpolation

- Chebyshev nodes: https://numpy.org/doc/stable/reference/generated/numpy.polynomial.chebyshev.chebpts1.html#numpy.polynomial.chebyshev.chebpts1
- Lagrange polynomial: https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.lagrange.html
