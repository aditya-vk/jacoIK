import numpy as np
from sympy import *

import dh
from functions import gmat, hmat, zxz

##################
# PRE PROCESSING #
##################

# Step 1 : Load DH Parameters
dza = dh.dispZ
dxa = dh.dispX
ax  = dh.angX

# Step 2 : Set delta and scale displacements
delta = 5.0*(sum(dxa) + sum(dza))
dz = dza/delta
dx = dxa/delta

# Step 3 : Find G and H for the given pose
pose   = np.matrix([[1.0,0.0,0.0,0.2],[0.0,0.0,-1.0,-0.15],[0.0,1.0,0.0,0.3]])
xT  = pose[0,3]/delta
yT  = pose[1,3]/delta
zT  = pose[2,3]/delta
(e1,e2,e3)    = zxz(pose[0:3,0:3])

rightG = gmat('x',0,xT,1)*gmat('y',0,yT,1)*gmat('z',0,zT,1)*gmat('z',e1,0,1)*gmat('x',e2,0,1)*gmat('z',e3,0,1)
rightH = hmat('x',0,xT,1)*hmat('y',0,yT,1)*hmat('z',0,zT,1)*hmat('z',e1,0,1)*hmat('x',e2,0,1)*hmat('z',e3,0,1)

# Step 4 : Determine the eight equations
# (a) : Initialize the unknowns in symbolic form
a1,a2,a3,a4,a5,a6 = symbols("a1:7")

# (b) : Generate the two equations
matG1 = gmat('z',a2,dz[1],1)*gmat('x',ax[1],dx[1],1)*gmat('z',a3,dz[2],1)*gmat('x',ax[2],dx[2],1)*gmat('z',a4,dz[3],1)*\
		gmat('x',ax[3],dx[3],1)*gmat('z',a5,dz[4],1)
matG2 = gmat('x',ax[0],dx[0],-1)*gmat('z',a1,dz[0],-1)*rightG*gmat('x',ax[5],dx[5],-1)*gmat('z',a6,dz[5],-1)*gmat('x',ax[4],dx[4],-1)
eqG  = matG1 - matG2

matH1 = hmat('z',a2,dz[1],1)*hmat('x',ax[1],dx[1],1)*hmat('z',a3,dz[2],1)*hmat('x',ax[2],dx[2],1)*hmat('z',a4,dz[3],1)*\
		hmat('x',ax[3],dx[3],1)*hmat('z',a5,dz[4],1)
matH2 = hmat('x',ax[0],dx[0],-1)*hmat('z',a1,dz[0],-1)*rightH*hmat('x',ax[5],dx[5],-1)*hmat('z',a6,dz[5],-1)*hmat('x',ax[4],dx[4],-1)
eqH  = matH1 - matH2

# (c) : Get the 8 equations from the quaternion matrices above (4+4)
eqG1 = expand(expand_trig(eqG[0,0]))
eqG2 = expand(expand_trig(eqG[0,1]))
eqG3 = expand(expand_trig(eqG[0,2]))
eqG4 = expand(expand_trig(eqG[0,3]))

eqH1 = expand(expand_trig(eqH[0,0]))
eqH2 = expand(expand_trig(eqH[0,1]))
eqH3 = expand(expand_trig(eqH[0,2]))
eqH4 = expand(expand_trig(eqH[0,3]))

# #######################
# # ELIMINATION METHODS #
# #######################

# Step 5 : Eliminate angles a1 and a6
# (a) Replace sinusoidals of a1 and a6, s1s6,s1c6,c1s6,c1c6 -> y1,y2,y3,y4
y1,y2,y3,y4 = symbols('y1:5')

eqG1 = (((eqG1.subs(sin(a1/2)*sin(a6/2),y1)).subs(sin(a1/2)*cos(a6/2),y2)).subs(cos(a1/2)*sin(a6/2),y3)).subs(cos(a1/2)*cos(a6/2),y4)
eqG1 = N(eqG1)
eqG2 = (((eqG2.subs(sin(a1/2)*sin(a6/2),y1)).subs(sin(a1/2)*cos(a6/2),y2)).subs(cos(a1/2)*sin(a6/2),y3)).subs(cos(a1/2)*cos(a6/2),y4)
eqG2 = N(eqG2)
eqG3 = (((eqG3.subs(sin(a1/2)*sin(a6/2),y1)).subs(sin(a1/2)*cos(a6/2),y2)).subs(cos(a1/2)*sin(a6/2),y3)).subs(cos(a1/2)*cos(a6/2),y4)
eqG3 = N(eqG3)
eqG4 = (((eqG4.subs(sin(a1/2)*sin(a6/2),y1)).subs(sin(a1/2)*cos(a6/2),y2)).subs(cos(a1/2)*sin(a6/2),y3)).subs(cos(a1/2)*cos(a6/2),y4)
eqG4 = N(eqG4)

eqH1 = (((eqH1.subs(sin(a1/2)*sin(a6/2),y1)).subs(sin(a1/2)*cos(a6/2),y2)).subs(cos(a1/2)*sin(a6/2),y3)).subs(cos(a1/2)*cos(a6/2),y4)
eqH1 = N(eqH1)
eqH2 = (((eqH2.subs(sin(a1/2)*sin(a6/2),y1)).subs(sin(a1/2)*cos(a6/2),y2)).subs(cos(a1/2)*sin(a6/2),y3)).subs(cos(a1/2)*cos(a6/2),y4)
eqH2 = N(eqH2)
eqH3 = (((eqH3.subs(sin(a1/2)*sin(a6/2),y1)).subs(sin(a1/2)*cos(a6/2),y2)).subs(cos(a1/2)*sin(a6/2),y3)).subs(cos(a1/2)*cos(a6/2),y4)
eqH3 = N(eqH3)
eqH4 = (((eqH4.subs(sin(a1/2)*sin(a6/2),y1)).subs(sin(a1/2)*cos(a6/2),y2)).subs(cos(a1/2)*sin(a6/2),y3)).subs(cos(a1/2)*cos(a6/2),y4)
eqH4 = N(eqH4)

print "Solving Linear System"

# (b) Solve the four equation to determine y1,y2,y3,y4
solys = solve([eqG1,eqG2,eqG3,eqG4],(y1,y2,y3,y4),simplify = False,quick = True,rational = False)

# print type(solys)
# print solys.keys()