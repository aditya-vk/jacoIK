# Define functions G, H, TMat

import numpy
from sympy import *


########################
# FUNCTION DEFINITIONS #
########################

# GMatrix Quaternion
def gmat(n,a,d,i):
	"""Takes input axis, angle of rotation, displacement
	and exponent to output quaternion matrix"""
	an = (a+d)/2
	if n is 'z':
		mat = numpy.matrix([[cos(an),-i*sin(an),0,0],[i*sin(an),cos(an),0,0],
							 [0,0,cos(an),i*sin(an)],[0,0,-i*sin(an),cos(an)]])
	elif n is 'y':
		mat = numpy.matrix([[cos(an),0,i*sin(an),0],[0,cos(an),0,i*sin(an)],
							 [-i*sin(an),0,cos(an),0],[0,-i*sin(an),0,cos(an)]])
	elif n is 'x':
		mat = numpy.matrix([[cos(an),0,0,i*sin(an)],[0,cos(an),-i*sin(an),0],
							 [0,i*sin(an),cos(an),0],[-i*sin(an),0,0,cos(an)]])
	for row in range(4):
		for col in range(4):
			mat[row,col] = expand_trig(mat[row,col])
	return mat



# HMatrix Quaternion
def hmat(n,a,d,i):
	"""Takes input axis, angle of rotation, displacement
	and exponent to output quaternion matrix"""
	an = (a-d)/2
	if n is 'z':
		mat = numpy.matrix([[cos(an),i*sin(an),0,0],[-i*sin(an),cos(an),0,0],
							 [0,0,cos(an),i*sin(an)],[0,0,-i*sin(an),cos(an)]])
	elif n is 'y':
		mat = numpy.matrix([[cos(an),0,-i*sin(an),0],[0,cos(an),0,-i*sin(an)],
							 [i*sin(an),0,cos(an),0],[0,i*sin(an),0,cos(an)]])
	elif n is 'x':
		mat = numpy.matrix([[cos(an),0,0,i*sin(an)],[0,cos(an),i*sin(an),0],
							 [0,-i*sin(an),cos(an),0],[-i*sin(an),0,0,cos(an)]])
	for row in range(4):
		for col in range(4):
			mat[row,col] = expand_trig(mat[row,col])
	return mat



# Transformation Matrix
def tmat(za,sz,xa,sx):
	"""Takes DH parameters to output transformation matrix
	for coordinate frame transformation"""
	mat = numpy.identity(4)
	for i in range(6):
		mat = mat*numpy.matrix([[cos(za[i]),-sin(za[i])*cos(xa[i]),-sin(za[i])*sin(xa[i]),sx[i]*cos(za[i])],
							[sin(za[i]), cos(za[i])*cos(xa[i]),-sin(xa[i])*cos(za[i]),sx[i]*sin(za[i])],
							[0,sin(za[i]),cos(za[i]),sz[i]],[0,0,0,1]])
	return mat



# ZXZ Euler Angles Extraction
def zxz(mat):
	"""Take pose and return the
	zxz euler angles for orientation"""
	a = numpy.arctan2(mat[0,2],-mat[1,2])
	b = numpy.arccos(mat[2,2])
	c = numpy.arctan2(mat[2,0], mat[2,1])
	angles = numpy.array([a,b,c])
	return angles

# Substitution Function
def multiSubs(exprNames, varNames, values=[]):

    if ( len(values) == 0 ):                # Get the values from the
        for varName in varNames:            # variables when not defined
            values.append( eval(varName) )  # as argument.
        # End for.
    # End if.

    for exprName in exprNames:                        # Create a temp copy
        expr = eval(exprName)                         # of each expression
        for i in range(len(varNames)):                # and substitute
            expr = expr.subs(varNames[i], values[i])  # each variable.
        # End for.
        yield expr     # Return each expression.
    # End for.