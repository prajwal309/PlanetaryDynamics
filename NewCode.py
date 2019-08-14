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

r_p = sympy.symbols('r_p', positive=True)              #distance of the planet

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
Orbit_Norm = [cos(f), sin(f), sympy.Integer(0)]
S_0 = [R*sin(theta)*cos(phi), R*sin(theta)*cos(phi), R*cos(theta)]
S_Norm = [sin(theta)*cos(phi+n*t), sin(theta)*cos(phi+n*t), cos(theta)]
S_p = [R*sin(theta)*cos(phi+n*t+Delta), R*sin(theta)*cos(phi+n*t+Delta), R*cos(theta)]

#Now describe relation
r_f = a*(1-e**2)/(1+e*cos(f))                            #distance in terms of true anomaly


#r prime in terms to true anamoly prime
r_p_f = a*(1-e**2)/(1+e*cos(f_p))
d_r_p_f__d_f_p = sympy.diff(r_p_f, f_p)
d_f_p__d_t = n*a**2*sympy.sqrt(1-e**2)/r_p**2        #True Anamoly (prime) relation with respect to the mean anamoly
d_r_p__d_t = d_r_p_f__d_f_p*d_f_p__d_t               #differentiating distance with respect to time
d_r_p__d_t = d_r_p__d_t.subs(r_p, a*(1-e**2)/(1+e*cos(f_p)))
d_r_p__d_t = sympy.powsimp(sympy.trigsimp(sympy.simplify(d_r_p__d_t)))

print("*"*100)
print("d_r_p__d_t")
sympy.pprint(d_r_p__d_t)
print("+"*100)
print("Latex version")
print(sympy.latex(d_r_p__d_t))
print("*"*100)


#Orbital element in terms of primed factors
Orbit_p_Norm = [cos(f_p), sin(f_p), sympy.Integer(0)]

#Surface element for primed factor
S_Norm = [sin(theta)*cos(phi+n*t), sin(theta)*sin(phi+n*t), cos(theta)]

#Surface element for primed factor
S_p_Norm = [sin(theta)*cos(phi+n*t+Delta), sin(theta)*sin(phi+n*t+Delta), cos(theta)]

x = sympy.symbols('x')
P2 = sympy.Rational(1,2)*(3*x**2 - 1)

alpha = sympy.symbols('alpha')
alpha_p = sympy.symbols('alpha_p')

beta = sympy.symbols('beta')
#U_T_r = -G*Ms*R**2/r**3*P2.subs(x,cos(alpha))
U_T_r_p = -G*Ms*R**2/r_p**3*P2.subs(x,cos(alpha_p))



cos_alpha_r_p = DotProduct(Orbit_p_Norm,S_p_Norm)
d_cos_alpha_r_p__d_t = sympy.diff(cos_alpha_r_p, f_p)*d_f_p__d_t + sympy.diff(cos_alpha_r_p, t)
#Subsituting beta = sqrt(1-e^2)
d_cos_alpha_r_p__d_t = d_cos_alpha_r_p__d_t.subs(sympy.sqrt(1-e*e), beta)
d_cos_alpha_r_p__d_t =  sympy.powsimp(sympy.trigsimp(sympy.simplify(d_cos_alpha_r_p__d_t)))



#Now differentiate the tidal force
FirstTerm = sympy.diff(U_T_r_p,r_p)*d_r_p__d_t

#Just differentiating with respect to
DivideBy = sympy.diff(cos(alpha_p),alpha_p)
SecondTerm = sympy.diff(U_T_r_p,alpha_p)/DivideBy*d_cos_alpha_r_p__d_t

#Now the differentiation term is given by:
d_U_T_r_p__d_t = FirstTerm + SecondTerm


print("*"*100)
print("d_U_T_r_p__d_t")
sympy.pprint(d_U_T_r_p__d_t)
print("+"*100)
print("Latex version")
print(sympy.latex(d_U_T_r_p__d_t))
print("*"*100)

#Calculate the angle alpha in terms of r and other function
cos_alpha_r = DotProduct(Orbit_Norm,S_Norm)
cos_alpha_r = sympy.powsimp(sympy.trigsimp(sympy.simplify(cos_alpha_r)))

U_T_r = -G*Ms*R**2/r**3*P2.subs(x,cos(alpha))
U_T_r = U_T_r.subs(cos(alpha), cos_alpha_r)
U_T_r = sympy.powsimp(sympy.trigsimp(sympy.simplify(U_T_r)))

print("*"*100)
print("U_T_r")
sympy.pprint(U_T_r)
print("+"*100)
print("Latex version")
print(sympy.latex(U_T_r))
print("*"*100)



#New set of variables for defining integrand
rho = sympy.symbols('rho')              #Density of the planet
h2 = sympy.symbols('h2')                #h2 is the vertical displace number
g = sympy.symbols('g')                  #g is the local acceleration due to gravity


Integrand = rho*h2/g*U_T_r*d_U_T_r_p__d_t
#simplifying the integrand
Integrand = sympy.powsimp(sympy.trigsimp(sympy.simplify(Integrand)))
Integrand = Integrand.subs(sympy.sqrt(1-e*e), beta)
print("*"*100)
print("Integrand")
sympy.pprint(Integrand)
print("+"*100)
print("Latex version")
print(sympy.latex(Integrand))
print("*"*100)

#Performing Integration in theta

print("*"*100)
print("Expression after first integration:")
Integration1 = sympy.integrate(Integrand*R*R*sin(theta), theta)
Integration1 = sympy.powsimp(sympy.trigsimp(sympy.simplify(Integration1)))
sympy.pprint(Integration1)

Result1 = Integration1.subs(theta,sympy.pi) - Integration1.subs(theta,0)
print("*"*100)
print("After the first integration")
#sympy.pprint(Result1)
print("Integration 1 in latex \n")
print(sympy.latex(Result1))
print("*"*100,"\n\n")

#Performing integration in phi
Integration2 = sympy.integrate(Result1, phi)
#Integration2 = sympy.powsimp(sympy.trigsimp(sympy.simplify(Integration2)))
FinalValue = Integration2.subs(phi, 2*sympy.pi) - Integration2.subs(phi, 0)
FinalValue = sympy.powsimp(sympy.trigsimp(sympy.simplify(FinalValue)))

print("*"*100)
print("Final expression that is yielded")
#sympy.pprint(FinalValue)
print("*"*100)
print("After integration in latex: Final Value--- Not simplified yet:: \n")
print(sympy.latex(FinalValue))
print("*"*100, "\n\n\n\n" )
