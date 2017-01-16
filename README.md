# Volterra equation
## Volterra equation of the second kind
## R language
## Numeric solution by methods - *finite sums* and *approximation*

---

## Parameters

|          **Parameter**         |      **Name**     |             **Example**             |
|:------------------------------:|:-----------------:|:-----------------------------------:|
| Free element of the equation   |     **FreeEl**    |            *function(g) exp( - g )* |
| Kernel                         |      **Core**     | *function(g, s) exp( - ( g - s ) )* |
| Count of steps in a grid       |      **step**     |                               *100* |
| Interval start point           |    **intstart**   |                                 *0* |
| Interval end point             |     **intend**    |                                *10* |
| X-axis limits                  |     **limits**    |                       *c(0.2, 0.5)* |
| Quadrature for **finite sums** |     **method**    |                       *"rectangle"* |
| Method for calculation         | **method_global** |                             *"sum"* |
| Error for **approximation**    |      **eps**      |                              *0.05* |

Supported quadrature formulas for **finite sums**:

1. Trapezoidal
2. Parabolic (Simpson's rule)
3. Rectangle (right)

---

## Usage

#### Finite sums

```R
f <- function(g) dnorm(g, 3, 0.5)
v <- volt(
	FreeEl = f, 
	Core = function(g, s) f(g-s), 
	step = 1000, 
	intend = 20, 
	method = "rectangle"
)
v$plot = v$plot + geom_hline(yintercept = 1/3, color = "red")
print(v$plot)
```

#### Approximation

```R
f <- function(g) dnorm(g, 3, 0.5)
v <- volt(
    FreeEl = f, 
    Core = function(g, s) f(g-s), 
    step = 1000, 
    intend = 20, 
    method_global = "approx",
    eps = 0.07
)
v$plot = v$plot + geom_hline(yintercept = 1/3, color = "red")
print(v$plot)
```

**Plot**:

<p align="center">
	<img src="https://raw.githubusercontent.com/hexeh/volterra/master/plot.png" alt = "Plot Image">
</p>

---

## Dependencies

[R-3.3.2](https://cran.r-project.org/)

### Packages:

- [ggplot2](https://cran.r-project.org/web/packages/ggplot2/)
- [Deriv](https://cran.r-project.org/web/packages/Deriv/)

---

### Todo

1. Fix error definition in finite sums
2. Time tracker
