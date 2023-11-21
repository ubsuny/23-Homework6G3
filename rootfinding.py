#!/usr/bin/env python
import math
import sys
import cmath
import numpy as np

def root_print_header(algorithm, accuracy):
    """Prints the header for the root finding process.

    Parameters
    ----------
    algorithm : str
        The name of the root finding algorithm used.
    accuracy : float
        The requested accuracy for the root finding process.
    """
    sys.stdout.write("\n ROOT FINDING using " + algorithm +
                     "\n Requested accuracy = " +repr(accuracy) +
                     "\n Step     Guess For Root          Step Size      " +
                     "     Function Value" +
                     "\n ----  --------------------  --------------------" +
                     "  --------------------" + "\n")

def root_print_step(step, x, dx, f_of_x):
    """Prints the information for each step of the root finding process.

    Parameters
    ----------
    step : int
        The number of the current step.
    x : float
        The current guess for the root.
    dx : float
        The current step size.
    f_of_x : float
        The value of the function at the current guess.
    """
    sys.stdout.write(repr(step).rjust(5))
    for val in [x, dx, f_of_x]:
        sys.stdout.write("  " + repr(val).ljust(20))
    sys.stdout.write("\n")

def root_max_steps(algorithm, max_steps):
    """Raises an exception when the maximum number of steps is exceeded.

    Parameters
    ----------
    algorithm : str
        The name of the root finding algorithm used.
    max_steps : int
        The maximum number of steps allowed.

    Raises
    ------
    Exception
        When the maximum number of steps is exceeded.
    """
    raise Exception(" " + algorithm + ": maximum number of steps " +
                    repr(max_steps) + " exceeded\n")

def root_simple(f, x, dx, accuracy=1.0e-6, max_steps=1000, root_debug=False):
    """Returns the root of f(x) given a guess x and a step dx with specified accuracy.

    Parameters
    ----------
    f : function
        The function to find the root of.
    x : float
        The initial guess for the root.
    dx : float
        The initial step size. Must have the same sign as (root - x).
    accuracy : float, optional
        The requested accuracy for the root finding process. Default is 1.0e-6.
    max_steps : int, optional
        The maximum number of steps allowed. Default is 1000.
    root_debug : bool, optional
        Whether to print the information for each step. Default is False.

    Returns
    -------
    x : float
        The final guess for the root.
    iterations : numpy.array
        An array of the guesses and function values for each step.
    step : int
        The number of steps taken.

    Raises
    ------
    Exception
        When the maximum number of steps is exceeded.
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
    return x,np.array(iterations), step

def root_bisection(f, x1, x2, accuracy=1.0e-6, max_steps=1000, root_debug=False):
    """Returns the root of f(x) in the interval bracketed by x1 and x2 with specified accuracy.

    Parameters
    ----------
    f : function
        The function to find the root of.
    x1 : float
        The lower bound of the interval. Must have opposite sign of f(x2).
    x2 : float
        The upper bound of the interval. Must have opposite sign of f(x1).
    accuracy : float, optional
        The requested accuracy for the root finding process. Default is 1.0e-6.
    max_steps : int, optional
        The maximum number of steps allowed. Default is 1000.
    root_debug : bool, optional
        Whether to print the information for each step. Default is False.

    Returns
    -------
    x_mid : float
        The final guess for the root.
    iterations : numpy.array
        An array of the guesses and function values for each step.
    step : int
        The number of steps taken.

    Raises
    ------
    Exception
        When f(x1) * f(x2) > 0.0 or when the maximum number of steps is exceeded.
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
    return x_mid,np.array(iterations), step

def root_secant(f, x0, x1, accuracy=1.0e-6, max_steps=20,root_debug=False):
    """Returns the root of f(x) given two guesses x0 and x1 with specified accuracy.

    Parameters
    ----------
    f : function
        The function to find the root of.
    x0 : float
        The first guess for the root.
    x1 : float
        The second guess for the root.
    accuracy : float, optional
        The requested accuracy for the root finding process. Default is 1.0e-6.
    max_steps : int, optional
        The maximum number of steps allowed. Default is 20.
    root_debug : bool, optional
        Whether to print the information for each step. Default is False.

    Returns
    -------
    x1 : float
        The final guess for the root.
    iterations : numpy.array
        An array of the guesses and function values for each step.
    num_steps : int
        The number of steps taken.

    Raises
    ------
    Exception
        When f(x0) = f(x1) or when the maximum number of steps is exceeded.
    """
    iterations=[]
    f0 = f(x0)
    dx = x1 - x0
    step = 0
    num_steps = 0 # Add a variable to store the number of steps
    if root_debug:        
        root_print_header("Secant Search", accuracy)
        root_print_step(step, x0, dx, f0)
        iterations.append([x0,f0])
    if f0 == 0:
        return x0, np.array(iterations), num_steps # Return x0 as the root and the number of steps
    while abs(dx) > abs(accuracy):
        f1 = f(x1)
        if f1 == 0:
            return x1, np.array(iterations), num_steps # Return x1 as the root and the number of steps
        if f1 == f0:
            raise Exception("Secant horizontal f(x0) = f(x1) algorithm fails")
        dx *= - f1 / (f1 - f0)
        x0 = x1
        f0 = f1
        x1 += dx
        step += 1
        num_steps += 1 # Increment the number of steps by one
        if step > max_steps:
            root_max_steps("root_secant", max_steps)
        if root_debug:
            root_print_step(step, x1, dx, f1)
            iterations.append([x1,f1])
    return x1, np.array(iterations), num_steps # Return the final outputs

def root_tangent(f, fp, x0, accuracy=1.0e-6, max_steps=20, root_debug=False):
    """Returns the root of f(x) with derivative fp = df(x)/dx
    given an initial guess x0, with specified accuracy.
    Uses Newton-Raphson (tangent) root-finding algorithm.

    Parameters
    ----------
    f : function
        The function to find the root of.
    fp : function
        The derivative of the function f.
    x0 : float
        The initial guess for the root.
    accuracy : float, optional
        The requested accuracy for the root finding process. Default is 1.0e-6.
    max_steps : int, optional
        The maximum number of steps allowed. Default is 20.
    root_debug : bool, optional
        Whether to print the information for each step. Default is False.

    Returns
    -------
    x0 : float
        The final guess for the root.
    iterations : numpy.array
        An array of the guesses and function values for each step.
    num_steps : int
        The number of steps taken.

    Raises
    ------
    Exception
        When fp(x0) = 0 or when the maximum number of steps is exceeded.
    """
    iterations = []
    f0 = f(x0)
    fp0 = fp(x0)
    if fp0 == 0.0:
        raise Exception(" root_tangent df/dx = 0 algorithm fails")
    dx = - f0 / fp0
    step = 0
    num_steps = 0 # Add a variable to store the number of steps
    if root_debug:
        root_print_header("Tangent Search", accuracy)
        root_print_step(step, x0, dx, f0)
        iterations.append([x0,f0])
    # Check if f(x0) is zero
    if f0 == 0.0:
        return x0, np.array(iterations), num_steps # Return x0 as the root and the number of steps
    while True:
        fp0 = fp(x0)
        if fp0 == 0.0:
            raise Exception(" root_tangent df/dx = 0 algorithm fails")
        dx = - f0 / fp0
        x0 += dx
        f0 = f(x0)
        if abs(dx) <= accuracy or f0 == 0.0:
            return x0, np.array(iterations), num_steps # Return the root, the iterations, and the number of steps
        step += 1
        num_steps += 1 # Increment the number of steps by one
        if step > max_steps:
            root_max_steps("root_tangent", max_steps)
        if root_debug:
            root_print_step(step, x0, dx, f0)
            iterations.append([x0,f0])
    return x0, np.array(iterations), num_steps # Return the final outputs
