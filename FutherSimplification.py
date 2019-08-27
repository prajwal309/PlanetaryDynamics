import sympy

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

Numerator = -3*sympy.pi*G**2*Ms**2*R**6*h2*n*rho
Denominator=20*beta*g*r**3*r_p**5

X = sympy.symbols('x', positive=True)           #Allow for proper substitution

#terms with a^2 beta^2
Factor1 = 16*a**2*beta**2
FirstTerm_1 = sin(f - n*t)**2*sin(Delta - f_p + n*t)*cos(Delta - f_p + n*t)
SecondTerm_1 = sin(f - n*t)*sin(Delta - f_p + n*t)**2*cos(f - n*t)
ThirdTerm_1 = -sin(f - n*t)*cos(f - n*t)*cos(Delta - f_p + n*t)**2
FourthTerm_1 = -sin(Delta - f_p + n*t)*cos(Delta - f_p + n*t)*cos(f - n*t)**2
AllTerm_1 = FirstTerm_1+SecondTerm_1+ThirdTerm_1+FourthTerm_1
AllTerm_1 = sympy.powsimp(sympy.trigsimp(sympy.simplify(AllTerm_1)))
AllTerm_1 = Factor1*AllTerm_1

#terms with beta^2 r_p^2
Factor2 = 16*beta**2*r_p**2
FirstTerm_2 = -sin(f - n*t)**2*sin(Delta - f_p + n*t)*cos(Delta - f_p + n*t)
SecondTerm_2 = -sin(f - n*t)*sin(Delta - f_p + n*t)**2*cos(f - n*t)
ThirdTerm_2 = sin(f - n*t)*cos(Delta - f_p + n*t)**2*cos(f - n*t)
FourthTerm_2 = sin(Delta - f_p + n*t)*cos(Delta - f_p + n*t)*cos(f - n*t)**2
AllTerm_2 = FirstTerm_2+SecondTerm_2+ThirdTerm_2+FourthTerm_2
AllTerm_2 = sympy.powsimp(sympy.trigsimp(sympy.simplify(AllTerm_2)))
AllTerm_2 = Factor2*AllTerm_2

#terms with a e r_p
Factor3 = -3*a*e*r_p
FirstTerm_3 = sin(2*Delta+2*f-3*f_p) - sin(2*Delta + 2*f -f_p)
#SecondTerm_3 =  -sin(2*Delta - 2*f -3*f_p+4*n*t)+ sin(2*Delta - 2*f -f_p+4*n*t)
SecondTerm_3 = 2*cos(2*Delta -2*f+4*n*t - 2*f_p)*sin(f_p) #Using sin(A) - Sin(B) = 2 cos((A+B)/2) sin((A-B)/2)
AllTerm_3 = Factor3*(FirstTerm_3+SecondTerm_3)

#terms with a e r_p sin(f_p)
Factor4 =  a*e*r_p*sin(f_p)
FirstTerm_4 = 36*sin(f-n*t)**2*sin(Delta - f_p+n*t)**2
SecondTerm_4 =  12*sin(f-n*t)**2*cos(Delta - f_p+n*t)**2
ThirdTerm_4 =  -20*sin(f-n*t)**2 #-----
FourthTerm_4 = 12*sin(Delta - f_p+n*t)**2*cos(f-n*t)**2
FifthTerm_4 = -20*sin(Delta - f_p+n*t)**2#-----
SixthTerm_4 = 36*cos(f-n*t)**2*cos(Delta - f_p+n*t)**2
SeventhTerm_4 = -20*cos(f -n*t)**2 #-----
EighthTerm_4 = -20*cos(Delta - f_p+n*t)**2 #-----
NinthTerm_4 = + 20 #-----

Term_12 = FirstTerm_4 + SecondTerm_4
Term_12 = sympy.simplify(sympy.trigsimp(Term_12))
Term_46 = FourthTerm_4 + SixthTerm_4
Term_46 = sympy.simplify(sympy.trigsimp(Term_46))
Term_1246 = Term_12 + Term_46
Term_1246 = sympy.powsimp(sympy.simplify(sympy.trigsimp(Term_1246)))


#Adding term by terms
Term_37 = ThirdTerm_4+SeventhTerm_4
Term_37 = sympy.trigsimp(Term_37)
Term_589 = sympy.trigsimp(FifthTerm_4+EighthTerm_4+NinthTerm_4)  #equal to 0

#Note FifthTerm_4, EighthTerm_4, and NinthTerm_4 cancels out
AllTerm_4 = Term_1246+Term_37+Term_589
AllTerm_4 = sympy.trigsimp(sympy.simplify(sympy.expand(AllTerm_4)))

AllTerm_4 = AllTerm_4.subs(f - n*t, X)
AllTerm_4 = sympy.trigsimp(sympy.simplify(sympy.expand(AllTerm_4)))
AllTerm_4 = AllTerm_4.subs(X, f - n*t)
AllTerm_4 = Factor4*AllTerm_4


AllTerm_34 = AllTerm_3+AllTerm_4
AllTerm_34 = sympy.trigsimp(sympy.simplify(AllTerm_34))



print("First Term is::")
sympy.pprint(Factor1)
sympy.pprint(AllTerm_1)

print("\n\nSecond Term is::")
sympy.pprint(Factor2)
sympy.pprint(AllTerm_2)

print("\n\nThird Term and Fourth Term combined  is::")
sympy.pprint(AllTerm_34)

FactorOut = 2*sympy.pi*rho*h2/g*G**2*Ms**2*R**6*n/a**6

CollectionTerm = Numerator/Denominator*(AllTerm_1+AllTerm_2+AllTerm_34)
CollectionTerm = CollectionTerm/FactorOut
CollectionTerm = sympy.simplify(CollectionTerm)
CollectionTerm = sympy.simplify(CollectionTerm)

print("Factored Out term is:")
sympy.pprint(FactorOut)
print("And rest of the terms are:")
sympy.pprint(CollectionTerm)
Term1, Term2, Term3, Term4, Term5 = sympy.expand(CollectionTerm)

print("\n After expansion")
print("Term1")
sympy.pprint(Term1)
print("*"*50)
print(Term2)
sympy.pprint(Term2)
print(sympy.latex(CollectionTerm))
print("*"*50)

print(sympy.srepr(CollectionTerm))
