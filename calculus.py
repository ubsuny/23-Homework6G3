calculusHW6
# calculus.py
=======
"""
calculus.py

This library provides functions for numerical integration and root-finding algorithms.

Functions:
- simpson(f, a, b, n): Approximates the definite integral of a function using composite Simpson's rule.
- trapezoid(f, a, b, n): Approximates the definite integral of a function using composite trapezoidal rule.
- adaptive_trapezoid(f, a, b, acc, output=False): Uses adaptive trapezoidal method to compute the definite integral.
- root_simple(f, x, dx, accuracy=1.0e-6, max_steps=1000, root_debug=False): Returns the root of a function using simple search with step halving.
- root_bisection(f, x1, x2, accuracy=1.0e-6, max_steps=1000, root_debug=False): Returns the root of a function using bisection search.
- root_secant(f, x0, x1, accuracy=1.0e-6, max_steps=20, root_debug=False): Returns the root of a function using the secant algorithm.
- root_tangent(f, fp, x0, accuracy=1.0e-6, max_steps=20, root_debug=False): Returns the root of a function using the Newton-Raphson (tangent) algorithm.
- root_print_header(algorithm, accuracy): Prints the header for root-finding algorithms.
- root_print_step(step, x, dx, f_of_x): Prints the details of each iteration step in root-finding algorithms.
- root_max_steps(algorithm, max_steps): Raises an exception when the maximum number of steps is exceeded.

"""
import math
import sys
import cmath
import numpy as np
main

def trapezoid_rule(func, a, b, n):
    h = (b - a) / n
    result = 0.5 * (func(a) + func(b))
    for i in range(1, n):
        result += func(a + i * h)
    result *= h
    return result

def adaptive_trapezoid_rule(func, a, b, epsilon):
    n = 1
calculusHW6
    old_result = 0
    result = trapezoid_rule(func, a, b, n)
    while abs(result - old_result) > epsilon:
        n *= 2
        old_result = result
        result = trapezoid_rule(func, a, b, n)
    return result

def simpsons_rule(func, a, b, n):
    h = (b - a) / n
    result = func(a) + func(b)
    for i in range(1, n, 2):
        result += 4 * func(a + i * h)
    for i in range(2, n-1, 2):
        result += 2 * func(a + i * h)
    result *= h / 3
    return result
=======
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
    """Prints the header for the root-finding algorithm.

    Parameters:
    - algorithm (str): The name of the root-finding algorithm.
    - accuracy (float): The specified accuracy for the root.
    """
    sys.stdout.write("\n ROOT FINDING using " + algorithm +
                     "\n Requested accuracy = " +repr(accuracy) +
                     "\n Step     Guess For Root          Step Size      " +
                     "     Function Value" +
                     "\n ----  --------------------  --------------------" +
                     "  --------------------" + "\n")
    
def root_print_step(step, x, dx, f_of_x):
    """Prints the details of each iteration step in the root-finding algorithm.

    Parameters:
    - step (int): The current iteration step.
    - x (float): The guess for the root.
    - dx (float): The step size.
    - f_of_x (float): The function value at the current guess.
    """
    sys.stdout.write(repr(step).rjust(5))
    for val in [x, dx, f_of_x]:
        sys.stdout.write("  " + repr(val).ljust(20))
    sys.stdout.write("\n")

def root_max_steps(algorithm, max_steps):
    """Raises an exception when the maximum number of steps is exceeded.

    Parameters:
    - algorithm (str): The name of the root-finding algorithm.
    - max_steps (int): The maximum number of allowed steps.
    """
    raise Exception(" " + algorithm + ": maximum number of steps " +
                    repr(max_steps) + " exceeded\n")

def root_simple(f, x, dx, accuracy=1.0e-6, max_steps=1000, root_debug=False):
    """Return root of f(x) given guess x and step dx with specified accuracy.
    Step must be in direction of root: dx must have the same sign as (root - x).
    
    Parameters:
    - f (function): The function for which to find the root.
    - x (float): The initial guess for the root.
    - dx (float): The step size.
    - accuracy (float): The desired accuracy for the root.
    - max_steps (int): The maximum number of allowed steps.
    - root_debug (bool): If True, print intermediate results during computation.

    Returns:
    tuple: The root value and an array containing iteration details.
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

def root_bisection(f, x1, x2, accuracy=1.0e-6, max_steps=1000, root_debug=False):
    """
    Return the root of the function f(x) within the interval [x1, x2] with the specified accuracy.

    Parameters:
    - f: function
        The function for which the root is to be found.
    - x1, x2: float
        The interval [x1, x2] in which the root is to be found. Assumes that f(x) changes sign within this interval.
    - accuracy: float, optional (default=1.0e-6)
        The desired accuracy for the root.
    - max_steps: int, optional (default=1000)
        The maximum number of iterations allowed.
    - root_debug: bool, optional (default=False)
        If True, print debugging information for each iteration.

    Returns:
    - x_root: float
        The estimated root of the function f(x) within the specified interval.
    - iterations: numpy.ndarray, shape (n, 2)
        An array containing the iterations during the root-finding process. Each row represents [x, f(x)].

    Raises:
    - Exception:
        If f(x1) * f(x2) > 0.0, indicating that there is no root within the specified interval.
        If the maximum number of steps is reached.

    """
    iterations = []
    f1 = f(x1)
    f2 = f(x2)
    if f1 * f2 > 0.0:
        raise Exception("f(x1) * f(x2) > 0.0")
    x_mid = (x1 + x2) / 2.0
    f_mid = f(x_mid)
    dx = x2 - x1
    step = 0
    if root_debug:
        iterations = []
        root_print_header("Bisection Search", accuracy)
        root_print_step(step, x_mid, dx, f_mid)
        iterations.append([x_mid,f_mid])
    while abs(dx) > accuracy:
        if f_mid == 0.0:
            dx = 0.0
        else:
            if f1 * f_mid > 0:
                x1 = x_mid
                f1 = f_mid
            else:
                x2 = x_mid
                f2 = f_mid
            x_mid = (x1 + x2) / 2.0
            f_mid = f(x_mid)
            dx = x2 - x1
        step += 1
        if step > max_steps:
            warning = "Too many steps (" + repr(step) + ") in root_bisection"
            raise Exception(warning)
        if root_debug:
            root_print_step(step, x_mid, dx, f_mid)
            iterations.append([x_mid,f_mid])
    return x_mid, np.array(iterations)


def root_secant(f, x0, x1, accuracy=1.0e-6, max_steps=20, root_debug=False):
    """Return root of f(x) given guesses x0 and x1 with specified accuracy.
    Uses secant root-finding algorithm.

    Parameters:
    - f (function): The function for which to find the root.
    - x0 (float): The first initial guess for the root.
    - x1 (float): The second initial guess for the root.
    - accuracy (float): The desired accuracy for the root.
    - max_steps (int): The maximum number of allowed steps.
    - root_debug (bool): If True, print intermediate results during computation.

    Returns:
    tuple: The root value and an array containing iteration details.
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

def root_tangent(f, fp, x0, accuracy=1.0e-6, max_steps=20, root_debug=False):
    """Return root of f(x) with derivative fp = df(x)/dx
    given initial guess x0, with specified accuracy.
    Uses Newton-Raphson (tangent) root-finding algorithm.

    Parameters:
    - f (function): The function for which to find the root.
    - fp (function): The derivative of the function.
    - x0 (float): The initial guess for the root.
    - accuracy (float): The desired accuracy for the root.
    - max_steps (int): The maximum number of allowed steps.
    - root_debug (bool): If True, print intermediate results during computation.

    Returns:
    tuple: The root value and an array containing iteration details.
    """
    iterations = []
    f0 = f(x0)
    fp0 = fp(x0)
    if fp0 == 0.0:
        raise Exception(" root_tangent df/dx = 0 algorithm fails")
    dx = -f0 / fp0
    step = 0
    if root_debug:
        iterations.append([x0, f0])
    if f0 == 0.0:
        return x0
    while True:
        fp0 = fp(x0)
        if fp0 == 0.0:
            raise Exception(" root_tangent df/dx = 0 algorithm fails")
        dx = -f0 / fp0
        x0 += dx
        f0 = f(x0)
        if abs(dx) <= accuracy or f0 == 0.0:
            return x0
        step += 1
        if step > max_steps:
            root_max_steps("root_tangent", max_steps)
        if root_debug:
            iterations.append([x0, f0])
    return x0
 main
