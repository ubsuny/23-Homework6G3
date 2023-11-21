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