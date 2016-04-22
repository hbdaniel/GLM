# -*- coding: utf-8 -*-

from classFE import *

def b(x):
    return(x+1)
    
def c(y, x):
    return(y*x)
    
dist1 = FE(b, c)

class A(object):
    def foo(self, b):
        print "foo", b

class B(A):
    def foo(self, a):
        super(B, self).foo(a)

myB = B()
myB.foo(1)