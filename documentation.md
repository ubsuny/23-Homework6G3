# 1. Simpson's Rule (`simpson` function):

The function uses Simpson's rule to approximate the definite integral of `f` from `a` to `b`. Simpson's rule is a numerical integration technique that uses parabolic arcs to approximate the area under a curve.

```python
def simpson(f, a, b, n):
    h = (b - a) / n
    i = np.arange(0, n)

    s = f(a) + f(b) 
    s += 4 * np.sum(f(a + i[1::2] * h))
    s += 2 * np.sum(f(a + i[2:-1:2] * h))
    
    return s * h / 3`
```

