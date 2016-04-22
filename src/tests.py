# -*- coding: utf-8 -*-

from classFE import *
from distributionsFE import *
import sympy as sp
import numpy as np

theta = sp.symbols('theta')
phi = sp.symbols('phi')
y = sp.symbols('y')

ns = FENormal(sp.symbols('mu'), sp.symbols('sigma'), symbolic=1)
print ns.cYPhi(2)