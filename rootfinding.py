#!/usr/bin/env python
import math
import sys
import cmath
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

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

# Define a function to plot a function
def plot_function(f, x1, x2, name):
    """Plots a function f(x) on a given interval [x1, x2] and saves it as a .png file.

    Args:
        f (function): The function to be plotted.
        x1 (float): The lower bound of the interval.
        x2 (float): The upper bound of the interval.
        name (str): The name of the function and the file.

    Returns:
        matplotlib.pyplot: The pyplot object that contains the plot.

    """
    # Create an array of x values
    xvals = np.linspace(x1, x2, 1000)

    # Compute the corresponding y values
    yvals = f(xvals)
    y3vals = np.zeros(len(xvals)) # Create an array of zeros

    # Plot the function
    plt.plot(xvals, yvals)
    plt.plot(xvals, y3vals, color='red') # Plot the line y=0 in red

    # Set the labels for x and y axes
    plt.xlabel("x")
    plt.ylabel(f'{name}') # Change the y label to match the function

    # Create an array of x tick values
    xticks = np.arange(-5*np.pi, 5*np.pi + np.pi/2, np.pi/2)

    # Create an array of x tick labels
    xticklabels = [r'$-5\pi$', r'$-9\pi/2$', r'$-4\pi$', r'$-7\pi/2$', r'$-3\pi$', r'$-5\pi/2$', r'$-2\pi$', r'$-3\pi/2$', r'$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$', r'$5\pi/2$', r'$3\pi$', r'$7\pi/2$', r'$4\pi$', r'$9\pi/2$', r'$5\pi$']

    # Set the ticks for x axis in terms of pi
    plt.xticks(xticks, xticklabels)

    # Set the limits for x and y axis
    plt.xlim(x1, x2)
    plt.ylim(-10, 10) # Add this line

    # Save the plot as a .png file
    plt.savefig(f'{name}.png'); # Use f-string syntax and add a semicolon
    # Return the plt object
    return plt

def find_roots(f, a, b, tol):
    # Find the number of roots in the interval using the sign change method
    n = 0 # Initialize the number of roots
    x = np.linspace(a, b, 100) # Create an array of 100 points in the interval
    y = f(x) # Evaluate the function at the points
    for i in range(len(y)-1):
        if y[i] * y[i+1] < 0: # If there is a sign change between two consecutive points
            n += 1 # Increment the number of roots

    # Find the roots using the brentq function
    roots = [] # Initialize the list of roots
    for i in range(len(y)-1):
        if y[i] * y[i+1] < 0: # If there is a sign change between two consecutive points
            root = opt.brentq(f, x[i], x[i+1], xtol=tol, rtol=tol) # Find the root in the subinterval
            # Check if the function is tan(x)
            if np.isclose(f(root), math.tan(root), atol=tol):
                # Check if the root is close to an odd number times pi/2
                if not np.isclose(root % np.pi, np.pi/2, atol=tol):
                    roots.append(root) # Append the root to the list
                else:
                    n -= 1 # Decrement the number of roots
            else:
                roots.append(root) # Append the root to the list
    # Return the number and the list of roots
    return n, roots

def print_roots(n, root_deg, root, name):
  print(f'      Algorithms for the root of {name}:')
  print("------------------------------------------------")
  print(f'There are {n} actual roots (in degrees), \nthey are: {(root_deg)} \n')
  print(f'There are {n} actual roots (in radians), \nthey are: {(root)} \n')

def print_results(algorithms, functions, arguments):
  # Initialize the lists
  answers_in_deg = []
  answers_in_rad = []
  steps = []

  # Loop over the algorithms, functions, and arguments
  for i, (alg, func, args) in enumerate(zip(algorithms, functions, arguments)):
    # Print the algorithm name, answer, iterations, and step
    print(f"{i+1}.{alg}")
    print("------------------")

    # Call the function with the arguments and unpack the results
    answer, iterations, step = func(*args)

    # Format the answer as a string with 20 decimal places
    # answer = format(float(answer), ".20f")
    answers_in_rad.append(answer)
    answer_to_deg = np.degrees(answer)
    print(f'root in deg: {answer_to_deg}')
    print(f'root in rad: {answer:.20f} \nsteps required: {step} \n')
    answers_in_deg.append(answer_to_deg)
    steps.append(step)

  # Return the lists
  return answers_in_deg, answers_in_rad, steps

# Define a function named efficiency with five parameters
def efficiency(answers_in_rad, root, algorithms, steps, tolerance):
  # Create an empty list to store the digits values
  digits_list = []
  # Loop through each element of answers_in_rad and get their indices
  for i, answer in enumerate(answers_in_rad):
    # Get the algorithm name and the number of steps for the current answer
    alg = algorithms[i]
    stp = steps[i]
    # Check if the answer is close to any element in the root list
    if any(abs(answer - r) < tolerance for r in root):
      # Print valid and the solution with the algorithm name
      print(f"The root by {alg} is {answer} is valid with steps {stp}.") 
      # Find the closest element in the root list to the answer
      r = min(root, key=lambda x: abs(x - answer))
      # Count the digits that are accurate in the answer
      digits = 0
      for i in range(len(str(answer))):
        # Compare each character in the answer with the corresponding character in the root element
        if str(answer)[i] == str(r)[i]:
          # Increment the digits count
          digits += 1
        else:
          # Break the loop if there is a mismatch
          break
      # Print the number of accurate digits
      if digits == 0:
        # If correct digits is 0, print all digits are correct
        print("The number of correct digits is 0 -- almost of the digits are correct.\n") # Merge the two print statements
      else:
        # Otherwise, print the number of correct digits
        print(f"The number of correct digits is {digits}.\n")
    else:
      # Print invalid and the solution
      print(f"The root by {alg} is {answer} is invalid with steps {stp}.\n")
      # Set the digits value to a large number to indicate invalid answer
      digits = 9999
    # Append the digits value to the digits list
    digits_list.append(digits)

  # Calculate the efficiency for each algorithm using a formula
  efficiency_list = [steps[i] * 10 + digits_list[i] for i in range(len(algorithms))]
  # Find the minimum efficiency value
  efficiency = min(efficiency_list)
  # Find the index of the minimum efficiency value in the efficiency list
  index = efficiency_list.index(efficiency)
  # Access the corresponding search from the search list
  algorithms = algorithms[index]
  # Display the output
  print(f"Among the 4 searches, the efficient search is {algorithms}: "
      f"\nbecause it has fewer steps {steps[index]}, with a greater accuracy of "
      f"{digits_list[index]} digits.") 
  # Return the efficiency value
  return None
