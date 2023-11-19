# calculus.py

def trapezoid_rule(func, a, b, n):
    h = (b - a) / n
    result = 0.5 * (func(a) + func(b))
    for i in range(1, n):
        result += func(a + i * h)
    result *= h
    return result

def adaptive_trapezoid_rule(func, a, b, epsilon):
    n = 1
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
