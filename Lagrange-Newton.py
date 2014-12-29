# Use Gauss-Newton algorithm to fit a given model
# See http://en.wikipedia.org/wiki/Gauss-Newton_algorithm
from sympy import *
from sympy.abc import *

def Gradient(fun, vars):
	grad = []
	for var in vars:
		grad.append(fun.diff(var))
	return grad

#def Lagrange_Newton(eqs, state_vars, data_dict, init_params_dict, eps=1e-4, max_iter=100):
def Lagrange_Newton(eq, state_var, data_dict, init_params_dict, eps=1e-4, max_iter=100):
	J = [] #Jacobian
	rs = [] #residuals
	dlen = 0
	for key, val in data_dict.items():
		dlen = len(val)
		break
	unknown_params = init_params_dict.keys()
	#
	unknowns = []
	for i in range(dlen):
		unknowns.append(Symbol("y%d"%i))
	for i in range(dlen):
		unknowns.append(Symbol("lambda%d"%i))
	unknowns = unknowns + unknown_params
	print unknowns
	#
	L = 0
	data_y = data_dict[state_var]
	for i in range(dlen):
		L = L + (data_y[i] - unknowns[i])*(data_y[i] - unknowns[i])
	for i in range(dlen):
		lambda_i = unknowns[dlen+i]
		state_eq = eq.lhs.subs(state_var, unknowns[i])
		for key,val in data_dict.items():
			state_eq = state_eq.subs(key, val[i])
		L = L + lambda_i*state_eq
	#
	F = Gradient(L, unknowns)
	JF = []
	for f in F:
		JF.append(Gradient(f, unknowns))
	#print JF
	#
	A = Matrix(JF)
	unknowns_tuple = tuple(unknowns)
	print unknowns_tuple
	NJF = []
	for row in JF:
		NJFrow = []
		for e in row:
			NJFrow.append( lambdify(unknowns_tuple, e) )
		NJF.append(NJFrow)
	print NJF

	NF = []
	for f in F:
		NF.append( lambdify(unknowns_tuple, f) )
	print NF

	init_vals = []
	for idx in range(dlen):
		init_vals.append(data_dict[state_var][idx])
	for idx in range(dlen):
		init_vals.append(0.0)
	for param in unknown_params:
		init_vals.append(init_params_dict[param])
	print init_vals
	
	for it in range(max_iter):
		A = []
		for row in NJF:
			Atmp = []
			for e in row:
				Atmp.append( e(*init_vals) )
			A.append(Atmp)
		bb = []
		for f in NF:
			bb.append( f(*init_vals) )
		#print A
		#print bb
		A = Matrix(A)
		bb = Matrix(bb)
		Ab = A.row_join(bb)
		sol = solve_linear_system(Ab, *unknowns_tuple)
		print sol
		norm = 0
		for idx in range(len(unknowns_tuple)):
			upd = sol[unknowns_tuple[idx]]
			init_vals[idx] = init_vals[idx] - upd
			norm = norm + upd*upd
		if sqrt(norm) < eps:
			print "Number of iteratons: %d" % it
			break
	return init_vals



sol = Lagrange_Newton(
	Eq(y-a*x/(b+x), 0), y,
	{
		x:[0.038,0.194,0.425,0.626,1.253,2.500,3.740],
		y:[0.050,0.127,0.094,0.2122,0.2729,0.2665,0.3317]
	},
	{ a: 0.9, b: 0.2 } 
)
print sol

