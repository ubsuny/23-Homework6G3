def trapezoid_rule(func, a, b, n):
    """
    Approximates the definite integral of a function using the composite trapezoidal rule.

    Parameters:
    - func (function): The integrand function.
    - a (float): The lower limit of integration.
    - b (float): The upper limit of integration.
    - n (int): The number of subintervals for the composite trapezoidal rule.

    Returns:
    float: The approximate integral value.
    """
    h = (b - a) / n
    result = 0.5 * (func(a) + func(b))
    
    # Summation of function values at intermediate points
    for i in range(1, n):
        result += func(a + i * h)
    
    result *= h
    return result


def adaptive_trapezoid_rule(func, a, b, epsilon):
    """
    Uses the adaptive trapezoidal method to compute the definite integral of a function.

    Parameters:
    - func (function): The integrand function.
    - a (float): The lower limit of integration.
    - b (float): The upper limit of integration.
    - epsilon (float): The desired accuracy for stopping the integration.

    Returns:
    float: The approximate integral value.
    """
    n = 1
    old_result = 0
    result = trapezoid_rule(func, a, b, n)
    
    # Iteratively refine the result until desired accuracy is achieved
    while abs(result - old_result) > epsilon:
        n *= 2
        old_result = result
        result = trapezoid_rule(func, a, b, n)
    
    return result


def simpsons_rule(func, a, b, n):
    """
    Approximates the definite integral of a function using composite Simpson's rule.

    Parameters:
    - func (function): The integrand function.
    - a (float): The lower limit of integration.
    - b (float): The upper limit of integration.
    - n (int): The number of subintervals for the composite Simpson's rule.

    Returns:
    float: The approximate integral value.
    """
    h = (b - a) / n
    result = func(a) + func(b)
    
    # Summation of function values at odd and even intermediate points
    for i in range(1, n, 2):
        result += 4 * func(a + i * h)
    
    for i in range(2, n-1, 2):
        result += 2 * func(a + i * h)
    
    result *= h / 3
    return result
