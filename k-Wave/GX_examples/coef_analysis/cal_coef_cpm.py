# w^2 = c1 * k + c2 * k^2 + c3 * k^3 + c4 * (i*w) * k^2

from sympy import *

def gen_B_mat(g):
    B = zeros(6, 4)

    B[0:3, 0] = [[1], [1 - g], [-1 / 2 * g + 1 / 2 * g ** 2]]
    B[3:6, 0] = B[0:3, 0]
    B[0:3, 0] *= cos(-pi * g / 2)
    B[3:6, 0] *= sin(-pi * g / 2)

    B[0:3, 1] = [[1], [2 - 2 * g], [1 - 3 * g + 2 * g ** 2]]
    B[3:6, 1] = B[0:3, 1]
    B[0:3, 1] *= cos(-pi * g)
    B[3:6, 1] *= sin(-pi * g)

    B[0:3, 2] = [[1], [3 - 3 * g], [3 - 15 / 2 * g + 9 / 2 * g ** 2]]
    B[3:6, 2] = B[0:3, 2]
    B[0:3, 2] *= cos(-3 / 2 * pi * g)
    B[3:6, 2] *= sin(-3 / 2 * pi * g)

    B[0:3, 3] = [[1], [3 - 2 * g], [3 - 5 * g + 2 * g ** 2]]
    B[3:6, 3] = B[0:3, 3]
    B[0:3, 3] *= cos(pi / 2 - pi * g)
    B[3:6, 3] *= sin(pi / 2 - pi * g)

    B_cut = zeros(4, 4)
    B_cut[0:4, :] = B[0:4, :].copy()

    # B_cut[0:2, :] = B[0:2, :].copy()
    # B_cut[2:4, :] = B[3:5, :].copy()

    return B_cut

def write_coef(c_A_sol, sol_file):
    # Write the output solution file

    f = open(sol_file, 'w')

    # Output the analytical solution
    for Ai, Ai_sol in c_A_sol.items():
        f.write(str(Ai) + ' = ')
        f.write(str(Ai_sol) + '\n')

    f.write('\n')

    # Output the taylor series of the coefficients
    for Ai, Ai_sol in c_A_sol.items():
        Ai_sol_taylor = series(Ai_sol)
        f.write(str(Ai) + ' = ')
        f.write(str(Ai_sol_taylor) + '\n')

    f.close()

    return

def main():
    sol_file = 'A_sol_comp_2st.txt'
    g = symbols('gamma')
    A1, A2, A3, A4 = symbols('A1, A2, A3, A4')

    B = gen_B_mat(g)
    c_A = Matrix([A1, A2, A3, A4])
    RHS = B * c_A
    LHS = Matrix([1, 2, 1, 0])

    c_A_sol = solve(RHS - LHS, c_A)

    write_coef(c_A_sol, sol_file)
    return

main()
