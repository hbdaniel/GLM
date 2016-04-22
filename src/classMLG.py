# -*- coding: utf-8 -*-

from math import log, pi, factorial
from numpy import prod
import sympy as sp
from distributionsFE import *

class MLG(object):
	def __init__(self, y, x, family, link, m = None):
		self.y = y
		self.x = x
		self.family = family
		self.link = link
		self.n = y.size
		self.m = m
		self.dist = __defFamily__()

	def __defFamily__(self):
		"""
		If family is kwnon, defines the distribution.
		"""
		if(self.family == 'Normal'):
			return FENormal(0, 1)

		if(self.family == 'Poisson'):
			return FEPoisson(1)

		if(self.family == 'Binomial'):
			if(self.m == None):
				print 'Invalid model! Please define m'
				pass
			return FEBinom(1, 0.5)

		if(self.family == 'InvGaussian'):
			return FEInvGaussian(1, 1)

		if(self.family == 'Gamma'):
			return FEGamma(1,1)

		else
			print('SOON')

	