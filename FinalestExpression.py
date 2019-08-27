import sympy
from sympy import *

f = sympy.symbols('f', positive=True)                  #True Anomaly
t = sympy.symbols('t', positive=True)                  #Time
n = sympy.symbols('n', positive=True)                  #Period

f_p = sympy.symbols('f_p', positive=True)              #True Anomaly prime due to tidal distortion
Delta = sympy.symbols('Delta', positive=True)

R = sympy.symbols('R', positive=True)                  #Radius of the planet
r = sympy.symbols('r', positive=True)                  #distance of the planet
a = sympy.symbols('a', positive=True)                  #semi-major axis of the planet
e = sympy.symbols('e', positive=True)                  #eccentricity of the orbit0

G = sympy.symbols('G', positive=True)                  #Gravitation constant
Ms = sympy.symbols('Ms', positive=True)                #Mass of the sun
Mp = sympy.symbols('Mp', positive=True)                #mass of the planet

r_p = sympy.symbols('r_p', positive=True)              #distance of the planet
beta = sympy.symbols('beta', positive=True)

#New set of variables for defining integrand
rho = sympy.symbols('rho', positive=True)              #Density of the planet
h2 = sympy.symbols('h2', positive=True)                #h2 is the vertical displace number
g = sympy.symbols('g', positive=True)                  #g is the local acceleration due to gravity

sin = sympy.sin                                    #sine of the angle
cos = sympy.cos                                    #cosine function

AllTerm = Add(Mul(Rational(3, 5), Pow(Symbol('a', positive=True), Integer(8)), Symbol('beta', positive=True), Pow(Symbol('r', positive=True), Integer(-3)), Pow(Symbol('r_p', positive=True), Integer(-5)), sin(Add(Mul(Integer(2), Symbol('Delta', positive=True)), Mul(Integer(2), Symbol('f', positive=True)), Mul(Integer(-1), Integer(2), Symbol('f_p', positive=True))))), Mul(Integer(-1), Rational(3, 10), Pow(Symbol('a', positive=True), Integer(7)), Pow(Symbol('beta', positive=True), Integer(-1)), Symbol('e', positive=True), Pow(Symbol('r', positive=True), Integer(-3)), Pow(Symbol('r_p', positive=True), Integer(-4)), sin(Symbol('f_p', positive=True))), Mul(Rational(9, 20), Pow(Symbol('a', positive=True), Integer(7)), Pow(Symbol('beta', positive=True), Integer(-1)), Symbol('e', positive=True), Pow(Symbol('r', positive=True), Integer(-3)), Pow(Symbol('r_p', positive=True), Integer(-4)), sin(Add(Mul(Integer(2), Symbol('Delta', positive=True)), Mul(Integer(2), Symbol('f', positive=True)), Mul(Integer(-1), Integer(3), Symbol('f_p', positive=True))))), Mul(Integer(-1), Rational(9, 20), Pow(Symbol('a', positive=True), Integer(7)), Pow(Symbol('beta', positive=True), Integer(-1)), Symbol('e', positive=True), Pow(Symbol('r', positive=True), Integer(-3)), Pow(Symbol('r_p', positive=True), Integer(-4)), sin(Add(Mul(Integer(2), Symbol('Delta', positive=True)), Mul(Integer(2), Symbol('f', positive=True)), Mul(Integer(-1), Symbol('f_p', positive=True))))), Mul(Integer(-1), Rational(3, 5), Pow(Symbol('a', positive=True), Integer(6)), Symbol('beta', positive=True), Pow(Symbol('r', positive=True), Integer(-3)), Pow(Symbol('r_p', positive=True), Integer(-3)), sin(Add(Mul(Integer(2), Symbol('Delta', positive=True)), Mul(Integer(2), Symbol('f', positive=True)), Mul(Integer(-1), Integer(2), Symbol('f_p', positive=True))))))


r_p_f = a*(1-e*e)/(1 + e*cos(f_p))
f_p_nt = beta*(a/r)**2*(n*t+Delta)

print("*"*50)
sympy.pprint(f_p_nt)
print("*"*50)
sympy.pprint(r_p_f)
print("*"*50)

print("*"*50)
print("All of the terms is given by::")
pprint(AllTerm)
print("*"*50)


AllTerm = AllTerm.subs(r_p, r_p_f)
AllTerm = AllTerm.subs(f_p, f_p_nt)


####
AllTermDiff = diff(AllTerm, Delta)
AllTermDiff = AllTermDiff.subs(Delta, 0)
AllTermDiff = AllTermDiff.subs(t, 0)
print("*"*50)
print("Before simplification")
pprint(AllTermDiff)

print("*"*50)
print("After simplification")
AllTermDiff = sympy.powsimp(sympy.simplify(sympy.trigsimp(AllTermDiff)))
#AllTermDiff = sympy.expand(AllTermDiff)
pprint(AllTermDiff)
