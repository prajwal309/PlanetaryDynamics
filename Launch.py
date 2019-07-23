# __ double spacing refers to division sign

import numpy as np
import matplotlib.pyplot as plt
import sympy
import sympy.printing as printing

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
Delta = sympy.symbols('Delta', real=True)

#Function need to integrate is following:
Orbit = [r*cos(f), r*sin(f), sympy.Integer(0)]
S_0 = [R*sin(theta)*cos(phi), R*sin(theta)*cos(phi), R*cos(theta)]
S = [R*sin(theta)*cos(phi+n*t), R*sin(theta)*cos(phi+n*t), R*cos(theta)]
S_prime = [R*sin(theta)*cos(phi+n*t+Delta), R*sin(theta)*cos(phi+n*t+Delta), R*cos(theta)]


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
U_T_r = -G*Ms*R**2/r**3*P2.subs(x,cos(alpha))

#Calculate the angle alpha in terms of r and other function
alpha_r = DotProduct(Orbit_r,S)/(r*R)
alpha_r = alpha_r.subs(f,f_r)

alpha_r_prime = DotProduct(Orbit_r,S_prime)/(r*R)
alpha_r_prime = alpha_r_prime.subs(f,f_r)

print("*"*100)
print("alpha_r \n")
sympy.pprint(alpha_r)
print("In latex \n")
print(sympy.latex(alpha_r))
print("*"*100,"\n\n\n\n")

input("Wait here...")
print("*"*100)
print("Simplified alpha_r \n")
alpha_r = sympy.powsimp(sympy.trigsimp(sympy.simplify(alpha_r)))
sympy.pprint(alpha_r)
print("In latex \n")
print(sympy.latex(alpha_r))
print("*"*100,"\n\n\n\n")

#after simplification
cos_alpha_r = cos(alpha_r)


print("*"*100)
print("alpha_r_prime \n")
sympy.pprint(alpha_r_prime)
print("In latex \n")
print(sympy.latex(alpha_r_prime))
print("*"*100,"\n\n\n\n")


print("*"*100)
print("Simplified alpha_r_prime \n")
alpha_r_prime = sympy.powsimp(sympy.trigsimp(sympy.simplify(alpha_r_prime)))
sympy.pprint(alpha_r)
print("In latex: \n")
print(sympy.latex(alpha_r_prime))
print("*"*100,"\n\n\n\n")

#Taken after simplification
cos_alpha_r_prime = cos(alpha_r_prime)

print("*"*100)
print("U_T_r_prime \n")
U_T_r_prime = -G*Ms*R**2/r**3*P2.subs(x,cos_alpha_r_prime)
sympy.pprint(U_T_r_prime)
print("In latex \n")
print(sympy.latex(U_T_r_prime))
print("*"*100,"\n\n\n\n")

print("*"*100)
print("d_U_T_r_prime/d_r \n")
d_U_T_r_prime__d_r = sympy.diff(U_T_r_prime, r)
sympy.pprint(d_U_T_r_prime__d_r)
print("In latex \n")
print(sympy.latex(d_U_T_r_prime__d_r))
print("*"*100,"\n\n\n\n")


print("*"*100)
print("Simplifying after this")
print("d_U_T_r\'/d_t")
d_U_T_r_prime__d_t = d_U_T_r_prime__d_r*d_r__d_t
sympy.pprint(d_U_T_r_prime__d_t)
print("In latex \n")
print(sympy.latex(d_U_T_r_prime__d_t))
print("*"*100,"\n\n\n\n")

print("*"*100)
print("Simplified d_UT_r\'/d_t")
d_U_T_r_prime__d_t = sympy.powsimp(sympy.trigsimp(sympy.simplify(d_U_T_r_prime__d_t)))
sympy.pprint(d_U_T_r_prime__d_t)
print("In latex \n")
print(sympy.latex(d_U_T_r_prime__d_t))
print("*"*100,"\n\n\n\n")


#New set of variables for defining integrand
rho = sympy.symbols('rho')              #Density of the planet
h2 = sympy.symbols('h2')                #h2 is the vertical displace number
g = sympy.symbols('g')                  #g is the local acceleration due to gravity


print("*"*100)
print("Integrand \n")
Integrand = rho*h2/g*U_T_r*d_U_T_r_prime__d_t
sympy.pprint(Integrand)
print("In latex \n")
print(sympy.latex(Integrand))
print("*"*100,"\n\n\n\n")

print("*"*100)
print("Simplified Integrand")
Integrand = sympy.powsimp(sympy.trigsimp(sympy.simplify(Integrand)))
sympy.pprint(Integrand)
print("In latex \n")
print(sympy.latex(Integrand))
print("*"*100,"\n\n\n\n")

#Performing Integration in theta

print("*"*100)
print("Expression after first integration:")
Integration1 = sympy.integrate(Integrand*r*r*sin(theta), theta)
Result1 = Integration1.subs(theta,sympy.pi) - Integration1.subs(theta,0)
Result1 =  sympy.powsimp(sympy.trigsimp(sympy.simplify(Result1)))
sympy.pprint(Result1)
print("In latex \n")
printing.latex(Result1)
print("*"*100,"\n\n\n\n")

#Performing integration in phi
Integration2 = sympy.integrate(Result1, phi)

print("*"*100)
print("After Second Integration")
FinalValue = Integration.subs(Integration2, 2*sympy.pi) - Integration.subs(Integration2, 0)
sympy.pprint(FinalValue)
print("In latex \n")
printing.latex(FinalValue)
print("*"*100, "\n\n\n\n" )


FinalValue = sympy.powsimp(sympy.trigsimp(sympy.simplify(FinalValue)))
print("*"*100)
print("The simplified final expression is given by:")
sympy.pprint(FinalValue)
print("*"*100, "\n\n\n\n")
