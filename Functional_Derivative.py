from sympy import *
from sympy.abc import *
from IPython.display import display

init_printing(use_latex='mathjax')

u = Function('u')(x,y,z)
u0 = Function('u0')(x,y,z)
q = Function('q')(x,y,z)
q0 = Function('q0')(x,y,z)
f = Function('f')(x,y,z)
lamd = Function('lambda')(x,y,z)

reg_term = 0.1*0.5*(q-q0)*(q-q0)

def grad(f):
	g = []
	for i in f.args:
		g.append(f.diff(i))
	return g

def dot(v1, v2):
	r = 0
	for idx in range(len(v1)):
		r += v1[idx]*v2[idx]
	return r

def func_diff(F, f, df):
	alpha = Symbol('alpha')
	F = F.subs(f, f+alpha*df)
	dF = F.diff(alpha)
	return simplify(dF.subs(alpha, 0))

def func_grad(F, xs, dxs):
	G = []
	for i in range(len(xs)):
		G.append(func_diff(F, xs[i], dxs[i]))
	return G

L = 0.5 * (u-u0) * (u-u0) + reg_term + q*dot(grad(u), grad(lamd)) - f*lamd
print 'L='
display(L)

phi = Function('\phi')(x,y,z)
psi = Function('\psi')(x,y,z)
chi = Function('\chi')(x,y,z)
xs = (u, lamd, q)
dxs = (phi, psi, chi)
Lx = func_grad(L, xs, dxs)
print '(Lu, Llamd, Lq)='
display(Matrix(Lx))


du = Function('\delta u')(x,y,z)
dl = Function('\delta \lambda')(x,y,z)
dq = Function('\delta q')(x,y,z)
dxs2 = (du, dl, dq)
Lxx = []
for Lxi in Lx:
	Lxx.append(func_grad(Lxi, xs, dxs2))

Lxx = Matrix(Lxx)
display(Lxx)



