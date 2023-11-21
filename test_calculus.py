"""
This module provides unit tests for the calculus module.
"""

import unittest
import sys
from io import StringIO
from calculus import (
    simpson, trapezoid, adaptive_trapezoid, 
    root_simple, root_bisection, root_secant, 
    root_tangent, root_print_header
)

class TestCalculusFunctions(unittest.TestCase):
    """
    Test various functions from the calculus module.
    """
    
    def test_simpson_known_integral(self):
        """Test the simpson function with a known integral."""
        # Define the integrand function
        def f(x):
            return x**2

        # Integration limits
        a = 0
        b = 1

        # Number of subintervals (even for Simpson's rule)
        n = 100

        # Expected result (known integral of x^2 from 0 to 1)
        expected_integral = 1/3

        # Call the simpson function
        result = simpson(f, a, b, n)

        # Check if the result is close to the expected integral
        self.assertAlmostEqual(result, expected_integral, places=5)

    def test_trapezoid_known_integral(self):
        """Test the trapezoid function with a known integral."""
    
        # Define the integrand function
        def f(x):
            return x**3

        # Integration limits
        a = 0
        b = 1
  
        # Number of subintervals
        n = 100

        # Expected result (known integral of x^3 from 0 to 1)
        expected_integral = 1/4

        # Call the trapezoid function
        result = trapezoid(f, a, b, n)

        # Check if the result is close to the expected integral
        self.assertAlmostEqual(result, expected_integral, places=3)   # Adjust places as needed

    def test_adaptive_trapezoid_known_integral(self):
        """Test the adaptive_trapezoid function with a known integral."""
        # Define the integrand function
        def f(x):
            return x**2

        # Integration limits
        a = 0
        b = 1

        # Expected result (known integral of x^2 from 0 to 1)
        expected_integral = 1/3

        # Call the adaptive_trapezoid function
        result = adaptive_trapezoid(f, a, b, acc=1e-6)

        # Check if the result is close to the expected integral
        self.assertAlmostEqual(result, expected_integral, places=5)

    def test_root_simple_known_root(self):
        """Test the root_simple function with a known root."""
        # Define the function
        def f(x):
            return x**2 - 4

        # Initial guess and step size
        x = 1.0
        dx = 0.1

        # Expected result (known root)
        expected_root = 2.0

        # Call the root_simple function
        result, _ = root_simple(f, x, dx, accuracy=1e-6)

        # Check if the result is close to the expected root
        self.assertAlmostEqual(result, expected_root, places=5)

    def test_root_bisection_known_root(self):
        """Test the root_bisection function with a known root."""
        # Define the function
        def f(x):
            return x**2 - 4

        # Initial interval containing the root
        x1 = 0
        x2 = 3

        # Expected result (known root)
        expected_root = 2.0

        # Call the root_bisection function
        result, _ = root_bisection(f, x1, x2, accuracy=1e-6)

        # Check if the result is close to the expected root
        self.assertAlmostEqual(result, expected_root, places=5)

    def test_root_secant_known_root(self):
        """Test the root_secant function with a known root."""
        # Define the function
        def f(x):
            return x**2 - 4

        # Initial guesses
        x0 = 1.0
        x1 = 3.0

        # Expected result (known root)
        expected_root = 2.0

        # Call the root_secant function
        result, _ = root_secant(f, x0, x1, accuracy=1e-6)

        # Check if the result is close to the expected root
        self.assertAlmostEqual(result, expected_root, places=5)

    def test_root_tangent_known_root(self):
        """Test the root_tangent function with a known root."""
        # Define the function and its derivative
        def f(x):
            return x**2 - 4

        def fp(x):
            return 2*x

        # Initial guess
        x0 = 1.5

        # Expected result (known root)
        expected_root = 2.0

        # Call the root_tangent function
        result = root_tangent(f, fp, x0, accuracy=1e-6)

        # Check if the result is close to the expected root
        self.assertAlmostEqual(result, expected_root, delta=1e-6)

    def test_root_print_header_output(self):
        """Test if root_print_header produces the correct output."""
        # Redirect stdout to capture the printed output
        captured_output = StringIO()
        sys.stdout = captured_output

        # Call the root_print_header function
        algorithm = "Test Algorithm"
        accuracy = 1e-6
        root_print_header(algorithm, accuracy)

        # Get the captured output
        result = captured_output.getvalue()

        # Reset redirect.
        sys.stdout = sys.__stdout__

        # Define the expected output
        expected_output = (
            "\n ROOT FINDING using Test Algorithm"
            "\n Requested accuracy = 1e-06"
            "\n Step     Guess For Root          Step Size      "
            "     Function Value"
            "\n ----  --------------------  --------------------"
            "  --------------------\n"
        )

        # Check if the result matches the expected output
        self.assertEqual(result, expected_output)

if __name__ == "__main__":
    unittest.main()
