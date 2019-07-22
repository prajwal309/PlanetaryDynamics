import sympy

e = sympy.symbols('e', real=True)
E = sympy.symbols('E', real=True)
a = sympy.symbols('a', real=True)

sin = sympy.sin
cos = sympy.cos

x = a*cos(E)
y = a*sin(E)*sympy.sqrt(1-e*e)

Distance = sympy.sqrt((x-a*e)**2+y**2)
sympy.pprint(sympy.powsimp(sympy.simplify(Distance), force=True))


#Test with simple functions
#sympy.pprint(sympy.expand_power_base(sympy.sqrt(x*x*x*x+3*x*x*x*x)))
#sympy.pprint(sympy.powsimp(sympy.sqrt(x*x*x*x+3*x*x*x*x),force=True))
#a = sympy.symbols('a', real='True')
#sympy.pprint(sympy.powsimp(sympy.sqrt(a*a+ 3.0*a*a)))
