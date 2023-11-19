import math
import sys
import cmath
import numpy as np


def simpson(f, a, b, n):
    """Approximates the definite integral of the given function using composite Simpson's rule.

    Parameters:
    - f (function): The integrand function.
    - a (float): The lower limit of integration.
    - b (float): The upper limit of integration.
    - n (int): The number of subintervals for the composite Simpson's rule.

    Returns:
    float: The approximate integral value.
    """
    h = (b - a) / n
    i = np.arange(0, n)

    s = f(a) + f(b)
    s += 4 * np.sum(f(a + i[1::2] * h))
    s += 2 * np.sum(f(a + i[2:-1:2] * h))

    return s * h / 3

def trapezoid(f, a, b, n):
    """Approximates the definite integral of the given function using composite trapezoidal rule.

    Parameters:
    - f (function): The integrand function.
    - a (float): The lower limit of integration.
    - b (float): The upper limit of integration.
    - n (int): The number of subintervals for the composite trapezoidal rule.

    Returns:
    float: The approximate integral value.
    """
    h = (b - a) / n
    s = f(a) + f(b)
    i = np.arange(0, n)
    s += 2 * np.sum(f(a + i[1:] * h))
    return s * h / 2

def adaptive_trapezoid(f, a, b, acc, output=False):
    """Uses the adaptive trapezoidal method to compute the definite integral of the given function.

    Parameters:
    - f (function): The integrand function.
    - a (float): The lower limit of integration.
    - b (float): The upper limit of integration.
    - acc (float): The desired accuracy for stopping the integration.
    - output (bool): If True, print intermediate results during computation.

    Returns:
    float: The approximate integral value.
    """
    old_s = np.inf
    h = b - a
    n = 1
    s = (f(a) + f(b)) * 0.5
    if output:
        print("Number of Subintervals (N) = " + str(n + 1) + ", Approximate Integral = " + str(h * s))
    while abs(h * (old_s - s * 0.5)) > acc:
        old_s = s
        for i in np.arange(n):
            s += f(a + (i + 0.5) * h)
        n *= 2.
        h *= 0.5
        if output:
            print("Number of Subintervals (N) = " + str(n) + ", Approximate Integral = " + str(h * s))
    return h * s

def root_print_header(algorithm, accuracy):
    sys.stdout.write("\n ROOT FINDING using " + algorithm +
                     "\n Requested accuracy = " +repr(accuracy) +
                     "\n Step     Guess For Root          Step Size      " +
                     "     Function Value" +
                     "\n ----  --------------------  --------------------" +
                     "  --------------------" + "\n")
    
def root_print_step(step, x, dx, f_of_x):
    sys.stdout.write(repr(step).rjust(5))
    for val in [x, dx, f_of_x]:
        sys.stdout.write("  " + repr(val).ljust(20))
    sys.stdout.write("\n")

def root_max_steps(algorithm, max_steps):
    raise Exception(" " + algorithm + ": maximum number of steps " +
                    repr(max_steps) + " exceeded\n")

def root_simple(f, x, dx, accuracy=1.0e-6, max_steps=1000, root_debug=False):
    """Return root of f(x) given guess x and step dx with specified accuracy.
    Step must be in direction of root: dx must have same sign as (root - x).
    """
    f0 = f(x)
    fx = f0
    step = 0
    iterations = []
    if root_debug:        
        root_print_header("Simple Search with Step Halving", accuracy)
        root_print_step(step, x, dx, f0)
        iterations.append([x,f0])
    while abs(dx) > abs(accuracy) and f0 != 0.0:
        x += dx
        fx = f(x)
        if f0 * fx < 0.0:   # stepped past root
            x -= dx         # step back
            dx /= 2.0       # use smaller step
        step += 1
        if step > max_steps:
            root_max_steps("root_simple", max_steps)
        if root_debug:
            root_print_step(step, x, dx, fx)
            iterations.append([x,fx])
    return x,np.array(iterations)

def root_secant(f, x0, x1, accuracy=1.0e-6, max_steps=20, root_debug=False):
    """Return root of f(x) given guesses x0 and x1 with specified accuracy.
    Uses secant root-finding algorithm.
    """
    iterations = []
    f0 = f(x0)
    dx = x1 - x0
    step = 0
    if root_debug:
        iterations.append([x0, f0])
    if f0 == 0:
        return x0
    while abs(dx) > abs(accuracy):
        f1 = f(x1)
        if f1 == 0:
            return x1
        if f1 == f0:
            raise Exception("Secant horizontal f(x0) = f(x1) algorithm fails")
        dx *= -f1 / (f1 - f0)
        x0 = x1
        f0 = f1
        x1 += dx
        step += 1
        if step > max_steps:
            root_max_steps("root_secant", max_steps)
        if root_debug:
            iterations.append([x1, f1])
    return x1, np.array(iterations)


