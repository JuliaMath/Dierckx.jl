Dierckx.jl
==========

*Julia library for 1-d and 2-d splines*

[![.github/workflows/ci.yml](https://github.com/JuliaMath/Dierckx.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaMath/Dierckx.jl/actions/workflows/ci.yml)
[![Coverage Status](http://img.shields.io/coveralls/JuliaMath/Dierckx.jl.svg?style=flat-square)](https://coveralls.io/github/JuliaMath/Dierckx.jl)

This is a Julia wrapper for the
[dierckx](http://www.netlib.org/dierckx/index.html) Fortran library,
the same library underlying the spline classes in scipy.interpolate.
Some of the functionality here overlaps with
[Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl),
a pure-Julia interpolation package.  Take a look at it if you have a
use case not covered here.

All new development on `Dierckx.jl` will be for Julia v1.3 and above.
The `master` branch is therefore incompatible with earlier versions
of Julia.

### Features

- Implements B-splines (basis splines).
- Splines from first order to fifth order; default is third order (cubic).
- Fit and evaluate 1-d and 2-d splines on irregular grids.
- Fit and evaluate 2-d splines at unstructured points.
- Fit "smooth" (non-interpolating) splines with adjustable smoothing factor s.
- Derivatives, integrals and roots of 1-d splines.
- Parametric B-splines.

Install (Julia 1.6 and later)
-----------------------------

```julia
(v1.6) pkg> add Dierckx
```

(Type `]` to enter package mode.)

Example Usage
-------------

```julia
using Dierckx
```

Fit a 1-d spline to some input data (points can be unevenly spaced):

```julia
x = [0., 1., 2., 3., 4.]
y = [-1., 0., 7., 26., 63.]  # x.^3 - 1.
spl = Spline1D(x, y)
```

Evaluate the spline at some new points:

```julia
spl([1.5, 2.5])  # result = [2.375, 14.625]
spl(1.5)  # result = 2.375
```

Equivalent to the above:

```julia
evaluate(spl, [1.5, 2.5])
evaluate(spl, 1.5)
```

Evaluate derivative, integral, or roots:

```julia
derivative(spl, 1.5)  # derivate at x=1.5; result is 5.75
integrate(spl, 0., 4.)  # integrate from x=0 to x=4; result is 60.0
roots(spl)  # result is [1.0]
```

*Note that `roots()` only works for cubic splines (k=3).*

Fit a 2-d spline to data on a (possibly irregular) grid:

```julia
x = [0.5, 2., 3., 4., 5.5, 8.]
y = [0.5, 2., 3., 4.]
z = [1. 2. 1. 2.;  # size is (length(x), length(y))
     1. 2. 1. 2.;
     1. 2. 3. 2.;
     1. 2. 2. 2.;
     1. 2. 1. 2.;
     1. 2. 3. 1.]

spline = Spline2D(x, y, z)
```

*Note that if you consider `z` as a matrix, `x` refers to row
 coordinates and `y` refers to column coordinates.*

Evaluate at element-wise points:

```julia
xi = [1., 1.5, 2.3, 4.5, 3.3, 3.2, 3.]
yi = [1., 2.3, 5.3, 0.5, 3.3, 1.2, 3.]
zi = spline(xi, yi)  # 1-d array of length 7
zi = evaluate(spline, xi, yi)  # equivalent to previous line
```

Evaluate at grid spanned by input arrays:

```julia
xi = [1., 1.5, 2.3, 4.5]
yi = [1., 2.3, 5.3]
zi = evalgrid(spline, xi, yi)  # 2-d array of size (4, 3)
```


Reference
---------

### 1-d Splines

```julia
Spline1D(x, y; w=ones(length(x)), k=3, bc="nearest", s=0.0)
Spline1D(x, y, xknots; w=ones(length(x)), k=3, bc="nearest")
```

- Create a spline of degree `k` (1 = linear, 2 = quadratic, 3 = cubic,
  up to 5) from 1-d arrays `x` and `y`. The parameter `bc` specifies
  the behavior when evaluating the spline outside the support domain,
  which is `(minimum(x), maximum(x))`. The allowed values are
  `"nearest"`, `"zero"`, `"extrapolate"`, `"error"`.

  In the first form, the number and positions of knots are chosen
  automatically. The smoothness of the spline is then achieved by
  minimizing the discontinuity jumps of the `k`th derivative of the
  spline at the knots. The amount of smoothness is determined by the
  condition that `sum((w[i]*(y[i]-spline(x[i])))**2) <= s`, with `s` a
  given non-negative constant, called the smoothing factor. The number
  of knots is increased until the condition is satisfied. By means of
  this parameter, the user can control the tradeoff between closeness
  of fit and smoothness of fit of the approximation.  if `s` is too
  large, the spline will be too smooth and signal will be lost ; if
  `s` is too small the spline will pick up too much noise. in the
  extreme cases the program will return an interpolating spline if
  `s=0.0` and the weighted least-squares polynomial of degree `k` if
  `s` is very large.

  In the second form, the knots are supplied by the user. There is no
  smoothing parameter in this form. The program simply minimizes the
  discontinuity jumps of the `k`th derivative of the spline at the
  given knots.

```julia
evaluate(spl, x)
```

- Evaluate the 1-d spline `spl` at points given in `x`, which can be a
  1-d array or scalar. If a 1-d array, the values must be monotonically
  increasing.

```julia
derivative(spl, x; nu=1)
```

- Evaluate the `nu`-th derivative of the spline at points in `x`.

```julia
integrate(spl, a, b)
```

-  Definite integral of spline between `x=a` and `x=b`.


```julia
roots(spl; maxn=8)
```

- For cubic splines (`k=3`) only, find roots. Only up to `maxn` roots
  are returned. A warning is issued if the spline has more roots than
  the number returned.


### Parametric Splines

These are the B-spline representation of a curve through N-dimensional space.

```julia
ParametricSpline(X; s=0.0, ...)
ParametricSpline(u, X; s=0.0, ...)
ParametricSpline(X, knots, ...)
ParametricSpline(u, X, knots, ...)
```

- `X` is a 2-d array with size `(N, m)`: `N` is the number of dimensions
  of the space (must be between 1 and 10) and `m` is the number of points.
  `X[:, i]` gives the coordinates of the `i`th point.

- `u` is a 1-d array giving parameter values at each of the `m` points. If not
  given, values are calculated automatically.

- `knots` is a 1-d array giving user-specified knots, if desired.

Keyword arguments common to all constructor methods:

- `w`: weight applied to each point (length `m` 1-d array).
- `k`: Spline order (between 1 and 5; default 3).
- `bc`: Boundary condition (default `'nearest'`).
- `periodic`: Treat curve as periodic? (Default is `false`).


### 2-d Splines

```julia
Spline2D(x, y, z; w=ones(length(x)), kx=3, ky=3, s=0.0)
Spline2D(x, y, z; kx=3, ky=3, s=0.0)
```

- Fit a 2-d spline to the input data. `x` and `y` must be 1-d arrays.

  If `z` is also a 1-d array, the inputs are assumed to represent
  unstructured data, with `z[i]` being the function value at point
  `(x[i], y[i])`. In this case, the lengths of all inputs must match.

  If `z` is a 2-d array, the data are assumed to be gridded: `z[i, j]`
  is the function value at `(x[i], y[j])`. In this case, it is
  required that `size(z) == (length(x), length(y))`. (Note that when
  interpreting `z` as a matrix, `x` gives the row coordinates and `y`
  gives the column coordinates.)

```julia
evaluate(spl, x, y)
```

- Evaluate the 2-d spline `spl` at points `(x[i], y[i])`. Inputs can be
  Vectors or scalars. Points outside the domain of the spline are set to
  the values at the boundary.

```julia
evalgrid(spl, x, y)
```

- Evaluate the 2-d spline `spl` at the grid points spanned by the
  coordinate arrays `x` and `y`. The input arrays must be
  monotonically increasing. The output is a 2-d array with size
  `(length(x), length(y))`: `output[i, j]` is the function value at
  `(x[i], y[j])`. In other words, when interpreting the result as a
  matrix, `x` gives the row coordinates and `y` gives the column
  coordinates.

```julia
derivative(spl, x, y; nux = 1, nuy = 1)
```

- Evaluate the partial derivative of the 2-d spline `spl` at the
  grid points spanned by the coordinate arrays `x` and `y`.
  Note that `x` refers to the row coordinates, and `y` to the column
  ones. The order of the derivatives may be between `0` and `kx-1` for `x`
  and `0` and `ky-1` for `y`.

- integral of a 2-d spline over the domain `[x0, x1]*[y0, y1]`

```julia
integrate(spl, x0, x1, y0, y1)
```


Translation from scipy.interpolate
----------------------------------

The `Spline` classes in scipy.interpolate are also thin wrappers
for the Dierckx Fortran library. The performance of Dierckx.jl should
be similar or better than the scipy.interpolate classes. (Better for
small arrays where Python overhead is more significant.) The
equivalent of a specific classes in scipy.interpolate:

| scipy.interpolate class      | Dierckx.jl constructor method              |
| ---------------------------- | ------------------------------------------ |
| UnivariateSpline             | `Spline1D(x, y; s=length(x))`              |
| InterpolatedUnivariateSpline | `Spline1D(x, y; s=0.0)`                    |
| LSQUnivariateSpline          | `Spline1D(x, y, xknots)`                   |
| SmoothBivariateSpline        | `Spline2D(x, y, z; s=length(x))`           |
| LSQBivariateSpline           |                                            |
| RectBivariateSpline          | `Spline2D(x, y, z; s=0.0)` (z = 2-d array) |
| SmoothSphereBivariateSpline  |                                            |
| LSQSphereBivariateSpline     |                                            |
| RectSphereBivariateSpline    |                                            |

Parametric splines:

| scipy.interpolate function   | Dierckx.jl constructor method              |
| ---------------------------- | ------------------------------------------ |
| `splprep(X)`                 | `ParametricSpline(X)`                      |
| `splprep(X, u=...)`          | `ParametricSpline(u, X)`                   |
| `splprep(X, t=...)`          | `ParametricSpline(X, t)`  (t = knots)      |
| `splprep(X, u=..., t=...)`   | `ParametricSpline(u, X, t)`                |


License
-------

Dierckx.jl is distributed under a 3-clause BSD license. See LICENSE.md
for details. The real*8 version of the Dierckx Fortran library as well as
some test cases and error messages are copied from the scipy package,
which is distributed under this license.


Citation
--------

If you use this package in a pulication and wish to cite it, you may want to cite the original [Fortran dierckx library](http://www.netlib.org/dierckx/):

> Paul Dierckx, Curve and Surface Fitting with Splines, Oxford University Press, 1993

If convenient, you can also include a link to the github repository for this package: https://github.com/JuliaMath/Dierckx.jl.
