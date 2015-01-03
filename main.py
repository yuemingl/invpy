from sympy import *
from sympy.abc import *
from IPython.display import display

init_printing(use_latex='mathjax')

class InvSolve:
	def __init__(self, equation, data_set):
		print equation.lhs.as_coefficient(y)
		print equation.rhs.as_coefficients_dict()	

a,b,c,x,y = symbols('a b c x y')
rlt = InvSolve(Eq(y, a*x**2+b*x+c), {x:[1,2,3,4,5],y:[3,4,2,4,1]})


