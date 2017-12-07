from sympy import *
x, z, A1, A2, A3, A4, A5, A6 = symbols('x, z, A1, A2, A3, A4, A5, A6')
c_z0 = -(1 + x)**2
c_z1 = A1 + I * A4 * (1 + x)
c_z2 = A2 + I * A5 * (1 + x)
c_z3 = A3 + I * A6 * (1 + x)
eqn = c_z0 + c_z1 * z + c_z2 * z**2 + c_z3 * z**3
z_sol = solve(eqn, z)
z1_sol = z_sol[0]
z2_sol = z_sol[1]
z3_sol = z_sol[2]
# z1_te = series(z1_sol, x)
print(z1_sol)
# print(z1_te)
# print(series(z1_sol))
print('abc')

