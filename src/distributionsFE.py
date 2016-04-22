# -*- coding: utf-8 -*-

from math import log, pi, factorial
from numpy import prod
from classFE import *

"""
Special cases with usual members of exponential family:
	- Normal
	- Poisson
	- Binomial
	- Inverse Gaussian
	- Gamma
"""

class FENormal(FE2):
	"""
	Y ~ Normal(mu, sigma); mu in R, sigma > 0
	E(Y) = mu
	Var(Y) = sigma
	theta = mu
	b(mu) = 0.5*mu**2 => b(theta) = 0.5*theta**2
	phi = 1/sigma**2; constant
	c(y, phi) = 0.5*(log(phi/2pi) + y**2)
	"""
	def __init__(self, mu, sigma):
		self.mu = mu
		self.sigma = sigma
		self.theta = mu
		self.phi = 1.0/sigma**2
		distribution = 'Normal'

	def expectedValue(self):
		return(self.mu)

	def variance(self):
		return(self.sigma)

	def bTheta(self):
		return((self.theta**2)/2.0)

	def cYPhi(self, y):
		return(0.5 * (log(self.phi) - log(2*pi) - (y*self.phi)**2))

	def canonicalLink(self):
		return self.mu

class FEPois(FE2):
	"""
	Y ~ Poisson(mu); mu > 0
	E(Y) = mu
	Var(Y) = mu
	theta = log(mu)
	b(mu) = mu => b(theta) = exY(theta)
	phi = 1; constant
	c(y, phi) = log(y!) = sum(log(i), i = 1,...,n)
	"""
	def __init__(self, mu):
		self.mu = mu
		self.phi = 1
		self.theta = log(mu)
		distribution = 'Poisson'
	
	def expectedValue(self):
		return(mu)

	def variance(self):
		return(mu)

	def bTheta(self):
		return(exp(self.theta))

	def cYPhi(self, y):
		return(log(prod(y)))

	def canonicalLink(self):
		return log(self.mu)

class FEBinom(FE2):
	"""
	Y ~ Binomial(n, mu); mu in [0, 1], n in Z+
	Y* = Y/n
	E(Y) = n*mu
	Var(Y) = n*mu*(1-mu)
	theta = log(mu/(1-mu))
	b(mu) = -log(1-mu) => b(theta) = log(1 + exp(theta))
	phi = n; constant
	c(y, phi) = log(Cn,y)
	"""
	def __init__(self, n, mu):
		self.mu = mu
		self.n = n
		self.theta = log(p/(1-p))
		self.phi = n
		distribution = 'Binomial'
	
	def expectedValue(self):
		return(self.n*self.mu)

	def variance(self):
		return(self.n*self.mu*(1-self.mu))

	def bTheta(self):
		return(log(1+exp(self.theta)))

	def cYPhi(self, y):
		return(log(factorial(self.phi) / factorial(self.phi - y) / factorial(y)))

	def canonicalLink(self):
		return log(self.mu/(1-self.mu))

class FEInvGaussian(FE2):
	"""
	Y ~ NInv(mu, phi); mu > 0, phi > 0
	E(Y) = mu
	Var(Y) = (mu**3)/phi
	theta = -1/2mu**3
	b(mu) = -1/mu**2 => b(theta) = -sqrt(-2*theta)
	c(y, phi) = 0.5*(log(phi) - log(2pi*y**3) - 1/y)
	"""
	def __init__(self, mu, phi):
		self.mu = mu
		self.theta = mu
		self.phi = phi
		distribution = 'Inverse Gaussian'

	def bTheta(self):
		return(-1/(self.mu**2))

	def cYPhi(self, y):
		return(0.5*(log(self.phi) - log(2*pi*(y**3)) - 1/y))

	def canonicalLink(self):
		return 1/self.mu**2

class FEGamma(FE2):
	"""
	Y ~ gamma(a, b); a > 0, b > 0
	E(Y) = ab
	Var(Y) = ab**2
	mu = ab
	phi = a = mu/b
	theta = -1/mu
	b(theta) = -log(-theta)
	c(y, phi) = phi*log(y*phi) - log(y) - log(gamma(phi))
	"""
	def __init__(self, a, b):
		self.mu = a*b
		self.theta = -1/mu
		self.phi = a
		distribution = 'Gamma'

	def expectedValue(self):
		return(self.mu)

	def variance(self):
		return(self.phi*self.mu)

	def bTheta(self):
		return(-log(-self.theta))

	def cYPhi(self, y):
		return(self.phi*log(y*self.phi) - log(y) - log(gamma(self.phi)))

	def canonicalLink(self):
		return 1/self.mu
