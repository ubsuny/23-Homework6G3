# Numerical Integration Documentation [1](https://faculty.ksu.edu.sa/sites/default/files/numerical_analysis_9th.pdf), [2](https://github.com/zhufengGNSS/Numerical-Recipes-1/blob/master/Numerical%20Recipes%20-%20The%20Art%20of%20Scientific%20Computing%20-%203rd%20Edition.pdf)

A short summary of the numerical integration methods implemented in the `calculus.py` module, along with examples of their application to specific functions.

## Integration Methods

### Trapezoid Rule

The trapezoid rule is a numerical integration method that approximates the definite integral of a function by dividing the area under the curve into trapezoids.

### Adaptive Trapezoid Rule

The adaptive trapezoid rule is an improvement over the basic trapezoid rule. It dynamically adjusts the number of subintervals to achieve a more accurate result.

### Simpson's Rule

Simpson's rule is a numerical integration method that uses quadratic approximations to estimate the definite integral. It is known for its accuracy, especially with smooth functions.

## Outputs

The following examples demonstrate the application of these integration methods to specific functions:

### Function: exp

- True Integral: None
- Trapezoid Rule: 1.3207062460837202
- Adaptive Trapezoid Rule: 1.3207061872625472
- Simpson's Rule: 1.3207061863727958
  
![exp](https://github.com/Pranjal-Srivastava-2023/23-Homework6G3/assets/143828394/2a691b81-6f4a-4ea3-97d9-9a895a9fc378)

### Function: cos

- True Integral: None
- Trapezoid Rule: 1.7226573096829807
- Adaptive Trapezoid Rule: 1.729503366850286
- Simpson's Rule: 1.7249562895115587

  ![cos](https://github.com/Pranjal-Srivastava-2023/23-Homework6G3/assets/143828394/ae753b1b-0da1-43c3-8307-73f5e6bea1e1)

### Function: poly

- True Integral: None
- Trapezoid Rule: 0.6666666666666667
- Adaptive Trapezoid Rule: 0.6666666666666665
- Simpson's Rule: 0.6666666666666661

  ![pol](https://github.com/Pranjal-Srivastava-2023/23-Homework6G3/assets/143828394/b197feec-71b3-4bb4-89ab-3e9739750d2b)

  References:

[1.Burden, R. L., & Faires, J. D. (2013). Numerical analysis (9th ed.). Brooks/Cole.](https://faculty.ksu.edu.sa/sites/default/files/numerical_analysis_9th.pdf)
[2. Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (2007). Numerical recipes in C (3rd ed.). Cambridge University Press.(https://github.com/zhufengGNSS/Numerical-Recipes-1/blob/master/Numerical%20Recipes%20-%20The%20Art%20of%20Scientific%20Computing%20-%203rd%20Edition.pdf)

