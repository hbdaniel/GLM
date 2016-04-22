# -*- coding: utf-8 -*-

from math import exp


class FE(object):
	def __init__(self, bTheta, cYPhi, distribution = None, theta = None, phi = None):
		"""
		bTheta and cYPhi are functions of one (theta) and two (y, phi), respectively.
		theta is the localization parameter
		phi is the dispersion parameter
		distribution holds a string with the name of the distribution
		thetaInv and phiInv are functions which gives theta and phi parameter as function of the original parameters
		"""
		self.bTheta = bTheta
		self.cYPhi = cYPhi
		self.theta = theta
		self.phi = phi

	
	def __definedParameters__(self):
		"""
		This method returns True if both theta and phi parameters are definde, else it returns False.
		"""
		if(self.theta == None and self.phi == None):
			return(False)
		else:
			return(True)


	def defineParameters(self, theta, phi):
		"""
		This method is used for define or chance the parameters theta and phi.
		If the parameters are changed.
		"""
		self.theta = theta
		self.phi = phi


	#TODO: Implement the expected value using symPy lib
	def expectedValue(self, theta):
		"""
		Expected Value of distribution as a function of theta.
		We know that: E(Y) = b'(theta).
		"""


	#TODO: Implement the expected value using symPy lib
	def variance(self, theta, phi):
		"""
		Expected Value of distribution as a function of theta and phi.
		We know that: var(Y) = b''(theta)/phi.
		"""


	def f(self, y):
		"""
		Evaluates the density at a point y.
		"""
		if(self.theta == None and self.phi == None):
			print('Please, define theta and phi.')

		else:
			return(exp(self.phi*(y*self.theta - self.bTheta()) + self.cYPhi(y)))

class FESample(FE):
	def likelihood(self, y, logL = 0):
		if(not logL): return(prod([self.f(i) for i in y]))
		else: return(sum([log(self.f(i)) for i in y]))
