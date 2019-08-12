# __ double spacing refers to division sign
import numpy as np
import matplotlib.pyplot as plt
import sympy

#printing
PPRINT_FLAG = False                #pretty print the  flag in sympy
LATEX_FLAG = True                  #print the variables in latex
SIMPLIFY_FLAG = False               #the simplify flag is given by

def DotProduct(A, B):
    #Function to take dot product of two vectors A and B of equal length
    assert len(A) == len(B)
    DotProduct = 0
    for i in range(len(A)):
        DotProduct+=A[i]*B[i]
    return DotProduct

f = sympy.symbols('f', positive=True)                  #True Anomaly
f_p = sympy.symbols('f_p', positive=True)              #True Anomaly prime due to tidal distortion
t = sympy.symbols('t', positive=True)                  #Time
n = sympy.symbols('n', positive=True)                  #Period

f_p = sympy.symbols('f_p', positive=True)              #True Anomaly prime due to tidal distortion
n_p = sympy.symbols('n_p', positive=True)              #Mean motion prime


G = sympy.symbols('G', positive=True)                  #Gravitation constant
Ms = sympy.symbols('Ms', positive=True)                #Mass of the sun
Mp = sympy.symbols('Mp', positive=True)                #mass of the planet

R = sympy.symbols('R', positive=True)                  #Radius of the planet
r = sympy.symbols('r', positive=True)                  #distance of the planet
a = sympy.symbols('a', positive=True)                  #semi-major axis of the planet
e = sympy.symbols('e', positive=True)                  #eccentricity of the orbit0

r_p = sympy.symbols('r_p', positive=True)                  #distance of the planet

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
S = [R*sin(theta)*cos(phi+n*t), R*sin(theta)*cos(phi+n*t), R*cos(theta)]

#Now describe relation
r_f = a*(1-e**2)/(1+e*cos(f))                            #distance in terms of true anomaly
Orbit_r = [Item.subs(r, r_f) for Item in Orbit]          #Orbit with r subsituted in terms of true anomaly
f_r = acos((a*(1- e*e)-r)/(e*r))


#No need to differentiate for normal orbit i.e. non primed orbit
#d_r_f__d_f = sympy.diff(r_f, f)                 #differentiating true anomaly in terms of tye anomaly
#d_f__d_t = n*a**2*sympy.sqrt(1-e**2)/r**2       #differentiating true anomaly in terms to time
#d_r__d_t = d_r_f__d_f*d_f__d_t                  #differentiating distance with respect to time


#Orbital element in terms of primed factors
Orbit_p = [r_p*cos(f_p), r_p*sin(f_p), sympy.Integer(0)]


#r prime in terms to true anamoly prime
r_p_f = a*(1-e**2)/(1+e*cos(f_p))
d_r_p_f__d_f_p = sympy.diff(r_p_f, f_p)
d_f_p__d_t = n*a**2*sympy.sqrt(1-e**2)/r_p**2        #True Anamoly (prime) relation with respect to the mean anamoly
d_r_p__d_t = d_r_p_f__d_f_p*d_f_p__d_t               #differentiating distance with respect to time

Orbit_r_p = [Item.subs(r_p, r_p_f) for Item in Orbit_p]          #Orbit prime  with r' subsituted in terms of true anomaly


if PPRINT_FLAG:
    sympy.pprint(d_f_p__d_t)


#Surface element for primed factor
S_p = [R*sin(theta)*cos(phi+n*t+Delta), R*sin(theta)*cos(phi+n*t+Delta), R*cos(theta)]

x = sympy.symbols('x')
P2 = sympy.Rational(1,2)*(3*x**2 - 1)

alpha = sympy.symbols('alpha')
alpha_p = sympy.symbols('alpha_p')

U_T_r = -G*Ms*R**2/r**3*P2.subs(x,cos(alpha))
U_T_r_p = -G*Ms*R**2/r_p**3*P2.subs(x,cos(alpha_p))



if SIMPLIFY_FLAG:
    U_T_r_p = sympy.powsimp(sympy.trigsimp(sympy.simplify(U_T_r_p)))

if PPRINT_FLAG:
    print("Differentiation of U(T) prime in terms of r(prime)--- before differentiation")
    sympy.pprint(U_T_r_p)

if LATEX_FLAG:
    print("U(T) prime in terms of r(prime) before differentiation and cos value substitution")
    print(sympy.latex(U_T_r_p))

#Calculate the angle alpha in terms of r and other function
cos_alpha_r = DotProduct(Orbit_r,S)/(r*R)
cos_alpha_r = cos_alpha_r.subs(f,f_r)


cos_alpha_r_p = DotProduct(Orbit_r_p,S_p)/(r_p*R)
cos_alpha_r_p = cos_alpha_r_p.subs(f,f_r)
alpha_r_p = acos(cos_alpha_r_p)



#Need to evaluate
d_alpha_p__d_f_p = sympy.diff(alpha_r_p, f_p)
d_alpha_p__d_t = d_alpha_p__d_f_p*d_f_p__d_t

if SIMPLIFY_FLAG:
    d_alpha_p__d_t = sympy.powsimp(sympy.trigsimp(sympy.simplify(d_alpha_p__d_t)))

if PPRINT_FLAG:
    print("Differentiation of alpha in terms of f_p")
    sympy.pprint(d_alpha_p__d_t)

if LATEX_FLAG:
    print("U(T) in terms of r_p")
    print(sympy.latex(d_alpha_p__d_t))
    print("*"*50+"\n\n")


#Now differentiating the Tidal force
d_U_T_r_p__d_t = sympy.diff(U_T_r_p, r_p)*d_r_p__d_t + sympy.diff(U_T_r_p, alpha_p)*d_alpha_p__d_t


#Replace to sin(alpha_p) with cos(alpha_p)
d_U_T_r_p__d_t = d_U_T_r_p__d_t.subs(sin(alpha_p),sympy.sqrt(1-cos(alpha_p)**2))

#Now substituting the value of cos(alpha'(r'))
d_U_T_r_p__d_t = d_U_T_r_p__d_t.subs(cos(alpha_p),cos_alpha_r_p)

if SIMPLIFY_FLAG:
    d_U_T_r_p__d_t = sympy.powsimp(sympy.trigsimp(sympy.simplify(d_U_T_r_p__d_t)))


if PPRINT_FLAG:
    print("Differentiation of U(T) prime with respect to time before substitution is given by")
    sympy.pprint(d_U_T_r_p__d_t)

if LATEX_FLAG:
    print("Differentiation of U(T) prime with respect to time before substitution is given by")
    print(sympy.latex(d_U_T_r_p__d_t))


#New set of variables for defining integrand
rho = sympy.symbols('rho')              #Density of the planet
h2 = sympy.symbols('h2')                #h2 is the vertical displace number
g = sympy.symbols('g')                  #g is the local acceleration due to gravity


print("*"*100)
print("Integrand \n")
Integrand = rho*h2/g*U_T_r*d_U_T_r_p__d_t

if PPRINT_FLAG:
    sympy.pprint(Integrand)
if LATEX_FLAG:
    print("In latex \n")
    print(sympy.latex(Integrand))
    print("*"*50,"\n\n\n\n")


print("Simplified Integrand")
if SIMPLIFY_FLAG:
    print("Simplifying the integrand.")
    Integrand = sympy.powsimp(sympy.trigsimp(sympy.simplify(Integrand)))
    print(sympy.latex(Integrand))
    print("*"*100,"\n\n\n\n")

#Performing Integration in theta

print("*"*100)
print("Expression after first integration:")
Integration1 = sympy.integrate(Integrand*R*R*sin(theta), theta)
Result1 = Integration1.subs(theta,sympy.pi) - Integration1.subs(theta,0)

if PPRINT_FLAG:
    sympy.pprint(Result1)
if LATEX_FLAG:
    print("Integration 1 in latex \n")
    print(sympy.latex(Result1))
    print("*"*100,"\n\n")

#Performing integration in phi
Integration2 = sympy.integrate(Result1, phi)

print("*"*100)
print("After Second Integration")
FinalValue = Integration2.subs(phi, 2*sympy.pi) - Integration2.subs(phi, 0)
if PPRINT_FLAG:
    sympy.pprint(FinalValue)
if LATEX_FLAG:
    print("After integration in latex: Final Value--- Not simplified yet:: \n")
    print(sympy.latex(FinalValue))
print("*"*100, "\n\n\n\n" )



print("Now simplifying the final value")
FinalValue = sympy.powsimp(sympy.trigsimp(sympy.simplify(FinalValue)))


#simplification in terms of r prime by substitution
#simplification in terms of substituting to r
print("*"*100)
print("The simplified final expression is given by:")
if PPRINT_FLAG:
    sympy.pprint(FinalValue)
if LATEX_FLAG:
    print("Simplified final expression: \n")
    print(sympy.latex(FinalValue))
print("*"*100, "\n\n\n\n")
