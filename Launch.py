# __ double spacing refers to division sign

import numpy as np
import matplotlib.pyplot as plt
import sympy

def DotProduct(A, B):
    #Function to take dot product of two vectors A and B of equal length
    assert len(A) == len(B)
    DotProduct = 0
    for i in range(len(A)):
        DotProduct+=A[i]*B[i]
    return DotProduct

f = sympy.symbols('f', positive=True)                  #True Anomaly
t = sympy.symbols('t', positive=True)              #Time
n = sympy.symbols('n', positive=True)                  #Period

G = sympy.symbols('G', positive=True)                  #Gravitation constant
Ms = sympy.symbols('Ms', positive=True)                #Mass of the sun
Mp = sympy.symbols('Mp', positive=True)                #mass of the planet

R = sympy.symbols('R', positive=True)                  #Radius of the planet
r = sympy.symbols('r', positive=True)                  #distance of the planet
a = sympy.symbols('a', positive=True)                  #semi-major axis of the planet
e = sympy.symbols('e', positive=True)                  #eccentricity of the orbit0


sin = sympy.sin                                    #sine of the angle
cos = sympy.cos                                    #cosine function
tan = sympy.tan                                    #tangent function
asin = sympy.asin                                  #arc of sine
acos = sympy.acos                                  #arc of cosine
atan = sympy.atan                                  #arc of tan

theta = sympy.symbols('theta', positive=True)          #colatitude for the surface of the planet
phi = sympy.symbols('phi', positive=True)              #longitude of


#Function need to integrate is following:
Orbit = [r*cos(f), r*sin(f), sympy.Integer(0)]
S_0 = [R*sin(theta)*cos(phi), R*sin(theta)*cos(phi), R*cos(theta)]
S = [R*sin(theta)*cos(phi+n*t), R*sin(theta)*cos(phi+n*t), R*cos(theta)]
S_prime = [R*sin(theta)*cos(phi+n*t), R*sin(theta)*cos(phi+n*t), R*cos(theta)]


#Now describe relation
r_f = a*(1-e**2)/(1+e*cos(f))                            #distance in terms of true anomaly
Orbit_r = [Item.subs(r, r_f) for Item in Orbit]          #Orbit with r subsituted in terms of true anomaly
d_r_f__d_f = sympy.diff(r_f, f)                 #differentiating true anomaly in terms of tye anomaly
d_f__d_t = n*a**2*sympy.sqrt(1-e**2)/r**2       #differentiating true anomaly in terms to time
d_r__d_t = d_r_f__d_f*d_f__d_t                  #differentiating distance with respect to time

f_r = acos((a*(1- e*e)-r)/(e*r))

x = sympy.symbols('x')
P2 = sympy.Rational(1,2)*(3*x**2 - 1)

alpha = sympy.symbols('alpha')
U_T = -G*Ms*R**2/r**3*P2.subs(x,cos(alpha))

#Calculate the angle alpha in terms of r and other function
alpha_r = DotProduct(Orbit_r,S)/(r*R)
alpha_r = alpha_r.subs(f,f_r)
alpha_r = sympy.simplify(sympy.trigsimp(sympy.simplify(alpha_r)))

cos_alpha_r = cos(alpha_r)

U_T_r = -G*Ms*R**2/r**3*P2.subs(x,cos(alpha_r))
U_T_r = sympy.simplify(sympy.trigsimp(U_T_r))

d_U_T_r__d_r = sympy.diff(U_T_r, r)

sympy.pprint(d_U_T_r__d_r)
d_U_T_r__d_r = sympy.simplify(sympy.trigsimp(d_U_T_r__d_r))
d_U_T_r__d_t = d_U_T_r__d_r*d_r__d_t
d_U_T_r__d_t = sympy.simplify(sympy.trigsimp(d_U_T_r__d_t))

print("Printing U_T_r__r")
sympy.pprint(d_U_T_r__d_t)
print("*"*100)
sympy.pprint(d_U_T_r__d_t)

#U_T_r__r = sympy.simplify(U_T_r__r)
#sympy.pprint("The potential function is given by:\n",U_T_r__r)

print("Now printing U_T_r")
sympy.pprint(U_T_r)

print("Now printing alpha")
sympy.pprint(alpha_r)
U_T__r = sympy.diff(U_T, r)

sympy.pprint(U_T)
