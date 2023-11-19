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
