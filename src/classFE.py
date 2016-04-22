# -*- coding: utf-8 -*-

from numpy import prod, exp, log
import sympy as sp

class FE(object):
	def __init__(self, bTheta, cYPhi, y = None, distribution = None, theta = None, phi = None):
		"""
		Generic Family Exponential class. This class uses symbolic computation.
		bTheta and cYPhi are functions of one (theta) and two (y, phi), respectively.
		param y: Sample data
		theta is the localization parameter
		phi is the dispersion parameter
		param distribution: holds a string with the name of the distribution
		thetaInv and phiInv are functions which gives theta and phi parameter as function of the original parameters
		"""
		self.bTheta = bTheta
		self.cYPhi = cYPhi
		self.y = y
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
		param theta: new value of theta
		param phi: new value of phi
		"""
		self.theta = theta
		self.phi = phi


	def expectedValue(self, rsymbolic = 0):
		"""
		Expected Value of distribution as a function of theta.
		We know that: E(Y) = b'(theta).
		param rsymbolic: if rsymbolic == True (1), returns the expected value as a function of theta
		"""
		if(self.__definedParameters__ and not rsymbolic):
			return(self.bTheta.diff().subs('theta', self.theta))
		else:
			return(self.bTheta.diff())
		
		
	def variance(self, rsymbolic = 0):
		"""
		Expected Value of distribution as a function of theta and phi.
		We know that: var(Y) = b''(theta)/phi.
		param rsymbolic: if rsymbolic == True (1), returns the expected value as a function of theta
		"""
		if(self.__definedParameters__ and not rsymbolic):
			return((self.bTheta.diff('theta', 2)/sp.symbols('phi')).subs([('theta', self.theta), ('phi', self.phi)]))
		else:
			return(self.bTheta.diff('theta', 2)/sp.symbols('phi'))


	def f(self, y):
		"""
		Evaluates the density at a point y.
		"""
		if(self.theta == None and self.phi == None):
			print('error: Please, define theta and phi.')

		else:
			return(exp(self.phi*(y*self.theta - self.bTheta.subs('theta', self.theta)) + self.cYPhi.subs([('y', y), ('phi', self.phi)])))

	def likelihood(self, logL = 0):
		if(not logL): return(prod([self.f(i) for i in self.y]))
		else: return(sum([log(self.f(i)) for i in self.y]))

# -*- coding: utf-8 -*-


class FE2(object):
	def __init__(self, bTheta, cYPhi, y = None, distribution = None, theta = None, phi = None):
		"""
		Family Exponential class used to the special cases from distributionsFE, which has the most usual distributions from the family.
		bTheta and cYPhi are functions of one (theta) and two (y, phi), respectively.
		param y: Sample data
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

		# FIXME: must return exp
		else:
			return(exp(self.phi*(self.theta*y - self.bTheta()) + self.cYPhi(y)))

	def likelihood(self, y, logL = 0):
		if(not logL): return(prod(self.f(y)))
		else: return(log(prod(self.f(y))))
