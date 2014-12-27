# Use Gauss-Newton algorithm to fit a model:
#		f = a*x/(b+x)
# See http://en.wikipedia.org/wiki/Gauss-Newton_algorithm
from sympy import *
from sympy.abc import *

xs = [0.038,0.194,0.425,0.626,1.253,2.500,3.740]
ys = [0.050,0.127,0.094,0.2122,0.2729,0.2665,0.3317]
J = [] #Jacobian
rs = [] #residuals
for idx in range(len(xs)):
	r = ys[idx] - a*xs[idx]/(b+xs[idx])
	J.append([lambdify((a,b),r.diff(a)), lambdify((a,b), r.diff(b))])
	rs.append(lambdify((a,b),r))

#Initial values for paramters
ai=0.9
bi=0.2
for iter in range(5):
	NJ = []
	for row in J:
		NJ.append([row[0](ai,bi), row[1](ai,bi)])
	A = Matrix(NJ)

	Nrs = []
	for r in rs:
		Nrs.append(r(ai,bi))
	b = Matrix(Nrs)

	Ab2 = (A.T*A).row_join(A.T*b)
	sol = solve_linear_system(Ab2,x,y)

	ai = ai - sol[x]
	bi = bi - sol[y]
	print ai, bi
