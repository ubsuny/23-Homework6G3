"""
test_calculus.py

This module provides unit tests for the calculus module.

"""
import unittest
import math
import numpy as np
from calculus import (
    simpson, trapezoid, adaptive_trapezoid, 
    root_simple, root_bisection, root_secant, 
    root_tangent, root_print_header, root_print_step, root_max_steps
)

class TestCalculusFunctions(unittest.TestCase):

    def test_simpson(self):
        """Test the simpson function with a known integral."""
        f = lambda x: x**2
        result = simpson(f, 0, 1, 1000)
        self.assertAlmostEqual(result, 1/3, places=5)

    def test_trapezoid(self):
        """Test the trapezoid function with a known integral."""
        f = lambda x: x**2
        result = trapezoid(f, 0, 1, 1000)
        self.assertAlmostEqual(result, 1/3, places=5)

    def test_adaptive_trapezoid(self):
        """Test the adaptive_trapezoid function with a known integral."""
        f = lambda x: x**2
        result = adaptive_trapezoid(f, 0, 1, 1e-6)
        self.assertAlmostEqual(result, 1/3, places=5)


