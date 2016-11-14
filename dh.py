# This file contains the DH paramters of the 6R manipulator considered.
from sympy import *
import numpy as np

########
# JACO #
########

dispZ = np.array([0.2755,0,0.0098,0.2393,0.03705,0.1321])
dispX = np.array([0,0.4,0,0,0,0])

angX  = np.array([pi/2,0,pi/2,pi/3,-pi/3,0])