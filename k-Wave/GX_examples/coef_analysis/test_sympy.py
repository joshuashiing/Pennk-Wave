from sympy import *

sol_file = 'A_sol.txt'
g = symbols('g')
A1, A2, A3, A4, A5, A6 = symbols('A1, A2, A3, A4, A5, A6')

B = zeros(6, 6)

B[0:3, 0] = [[1], [1 - g], [-1/2*g + 1/2*g**2]]
B[3:6, 0] = B[0:3, 0]
B[0:3, 0] *= cos(-pi * g / 2)
B[3:6, 0] *= sin(-pi * g / 2)

B[0:3, 1] = [[1], [2 - 2*g], [1 - 3*g + 2*g**2]]
B[3:6, 1] = B[0:3, 1]
B[0:3, 1] *= cos(-pi * g)
B[3:6, 1] *= sin(-pi * g)

B[0:3, 2] = [[1], [3 - 3*g], [3 - 15/2*g + 9/2*g**2]]
B[3:6, 2] = B[0:3, 2]
B[0:3, 2] *= cos(-3/2 * pi * g)
B[3:6, 2] *= sin(-3/2 * pi * g)

B[0:3, 3] = [[1], [2 - g], [1 - 3/2*g + 1/2*g**2]]
B[3:6, 3] = B[0:3, 3]
B[0:3, 3] *= cos(pi/2 - pi*g/2)
B[3:6, 3] *= sin(pi/2 - pi*g/2)

B[0:3, 4] = [[1], [3 - 2*g], [3 - 5*g + 2*g**2]]
B[3:6, 4] = B[0:3, 4]
B[0:3, 4] *= cos(pi/2 - pi*g)
B[3:6, 4] *= sin(pi/2 - pi*g)

B[0:3, 5] = [[1], [4 - 3*g], [6 - 21/2*g + 9/2*g**2]]
B[3:6, 5] = B[0:3, 5]
B[0:3, 5] *= cos(pi/2 - 3/2*pi*g)
B[3:6, 5] *= sin(pi/2 - 3/2*pi*g)

c_A = Matrix([A1, A2, A3, A4, A5, A6])
RHS = B * c_A
LHS = Matrix([1, 2, 1, 0, 0, 0])


n = 3
B = B[0:n, 0:n]
c_A = c_A[0:n, 0]
LHS = LHS[0:n, 0]
RHS = B * c_A
c_A_sol = solve(RHS - LHS, c_A)

f = open(sol_file, 'w')
for Ai, Ai_sol in c_A_sol.items():
    f.write(str(Ai) + ' = ')
    f.write(str(Ai_sol) + '\n')

f.write('\n')

for Ai, Ai_sol in c_A_sol.items():
    Ai_sol_taylor = series(Ai_sol)
    f.write(str(Ai) + ' = ')
    f.write(str(Ai_sol_taylor) + '\n')

f.close()