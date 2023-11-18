# "Numerical Integration Methods and Adaptive Techniques"

In computational physics, numerical integration is essential because it offers ways to approximate definite integrals in situations where analytical solutions are difficult or unfeasible. Three numerical integration methods are examined in this project: the adaptive trapezoidal method, the trapezoidal rule, and Simpson's rule. These techniques advance the study of numerical analysis by providing computational tools for approximating definite integrals.

## Objectives:
**1. Applying Numerical Integration Techniques:**
- Use Python to implement the trapezoidal rule, the Simpson's rule, and an adaptive trapezoidal technique.
- Recognize the underlying mathematical ideas of each technique and how numerical integration uses them.

**2. Accuracy and Efficiency Comparison:**
- Use Python to implement the trapezoidal rule, the Simpson's rule, and an adaptive trapezoidal technique.
- Recognize the underlying mathematical ideas of each technique and how numerical integration uses them.

**3. Utilize NumPy for Integration Comparison:**
- Incorporate the numerical integration functions that come with NumPy so that the custom implementations can be compared.
- Examine and talk about the results' parallels and discrepancies.

**4. Documentation and Testing:**
- Document the code by adding docstrings to functions for clarity and understanding.
- Implement unit tests to ensure the correctness of each numerical integration function.

**5. GitHub Actions Integration:**
- Integrate GitHub Actions for linting and unit testing to maintain code quality.
- Implement continuous integration practices to streamline development.


# Integral Functions in `calculus.py`:

## 1. Simpson's Rule (`simpson` function):

The function uses Simpson's rule to approximate the definite integral of `f` from `a` to `b`. Simpson's rule is a numerical integration technique that uses parabolic arcs to approximate the area under a curve.

```python
def simpson(f, a, b, n):
    h = (b - a) / n
    i = np.arange(0, n)

    s = f(a) + f(b) 
    s += 4 * np.sum(f(a + i[1::2] * h))
    s += 2 * np.sum(f(a + i[2:-1:2] * h))
    
    return s * h / 3
```
- `f`: This is the function to be integrated.
- `a` and `b`: These are the lower and upper bounds of integration.
- `n`: This is the number of subintervals used in the approximation.
- `h`: Represents the width of each subinterval.
- `i`: Creates an array of indices [0, 1, ..., n-1].
- `s`: Initializes the sum with the values of f(a) and f(b).
The main computation involves summing values of the function at different points with appropriate weights to get the integral approximation.


## 2. Trapezoidal Rule (`trapezoid` function):

```python
def trapezoid(f, a, b, n):
    h = (b - a) / n
    s = f(a) + f(b)
    i = np.arange(0, n)
    s += 2 * np.sum(f(a + i[1:] * h))
    return s * h / 2
```
Similar to Simpson's rule, this function approximates the definite integral of `f` from `a` to `b` using the trapezoidal rule.
The computation involves summing values of the function at different points, but here, the weights are different, reflecting the trapezoidal shape of each subinterval.

## 3. Adaptive Trapezoidal Method (`adaptive_trapezoid` function):

```python
def adaptive_trapezoid(f, a, b, acc, output=False):
    old_s = np.inf
    h = b - a
    n = 1
    s = (f(a) + f(b)) * 0.5
    if output:
        print("N = " + str(n + 1) + ",  Integral = " + str(h * s))
    while abs(h * (old_s - s * 0.5)) > acc:
        old_s = s
        for i in np.arange(n):
            s += f(a + (i + 0.5) * h)
        n *= 2.
        h *= 0.5
        if output:
            print("N = " + str(n) + ",  Integral = " + str(h * s))
    return h * s
```

This function uses an adaptive trapezoidal method to compute the definite integral of `f` from `a` to `v` to a desired accuracy acc.
The algorithm starts with a single trapezoid and iteratively doubles the number of subintervals until the desired accuracy is achieved.
The variable output controls whether intermediate results are printed.
The main loop keeps refining the estimate until the difference between the old and new estimates is within the specified accuracy.


