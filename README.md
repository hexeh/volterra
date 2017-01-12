# Volterra equation
## Volterra equation of the second kind
## R language

---

Supported quadrature formulas:

1. Trapezoidal
2. Parabolic (Simpson's rule)
3. Rectangle (right)

---

## Usage

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
Plot:
<p align="center">
	<img src="https://raw.githubusercontent.com/hexeh/volterra/master/plot.png" alt = "Plot Image">
</p>

## Dependencies

[R-3.3.2](https://cran.r-project.org/)

### Packages:

- [ggplot2](https://cran.r-project.org/web/packages/ggplot2/)
