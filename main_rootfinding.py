"""A script to find and compare the roots of different functions.

This script uses various root finding algorithms from the rootfinding module to find and compare the roots of the tan(x) and tanh(x) functions in a given interval. It plots the functions, prints the roots, and evaluates the efficiency of each algorithm.

Attributes:
    x1 (float): The lower bound of the interval.
    x2 (float): The upper bound of the interval.
    y3 (int): A dummy variable that is not used.
    dx (float): The step size for the simple search algorithm.
    xmid (float): The midpoint of the interval.
    x0 (float): The initial guess for the tangent and secant algorithms.
    tol (float): The tolerance for the root finding algorithms.
    algorithms (list): The list of algorithm names.
    functions (list): The list of function objects from the rootfinding module.
    funcs (list): The list of function expressions, derivatives, and names.

"""
from rootfinding import *

# Set the range of x and y values
x1, x2 = -np.pi, np.pi
y3 = 0
dx = 0.01
xmid = (x1+x2)*0.5
x0 = np.pi/100
tol = 1e-5

algorithms = ['Simple Search', 'Bisection Search', 'Tangent Search', \
              'Secant Search']

# Define a list of functions to call
functions = [root_simple, root_bisection, root_tangent, root_secant]

# Define a list of functions and their names
funcs = [(lambda x: np.tan(x), lambda x: 1/(np.cos(x))**2, 'tan(x)'), \
         (lambda x: np.tanh(x), lambda x: 1/(np.cosh(x))**2, 'tanh(x)')]

# Loop over the functions
for f, fp, name in funcs:
    # Plot the function
    plot_function(f, x1, x2, name)
    # Define a list of arguments to pass
    arguments = [(f, xmid, dx, tol, 1000, True), (f, x1, x2, tol, 1000, True), \
     (f, fp, x0, tol, 20, True), (f, x1, x2, tol, 20, True)]
    # Find the roots
    n, root = find_roots(f, x1, x2, tol)
    np.set_printoptions (formatter= {'float': ' {: .2f}'.format})
    root_deg = np.degrees(root)
    # Print the roots
    print_roots(n, root_deg, root, name)
    # Print the results
    answers_in_deg, answers_in_rad, steps = print_results(algorithms, \
                                                          functions, arguments)
    # Print the efficiency
    efficiency(answers_in_rad, root, algorithms, steps, tol)
