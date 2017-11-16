"""
Compute the coefficients of the wave equation based on the Kjartansson 
constant-Q model by taking advantage of the Taylor expansion up to O(x^2). 
Conduct coefficient analysis for this approximate visco-acoustic wave
equation and determine the principal component.

Author:         Guangchi Xing
Affiliation:    Penn State University
Date:           11/16/2017 

"""

from sympy import *

sol_file = 'A_sol.txt'  # Solution Output File


def gen_B_mat(g):
    # Generate the B matrix given a specific gamma symbol

    B = zeros(6, 6)

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

    B[0:3, 3] = [[1], [2 - g], [1 - 3 / 2 * g + 1 / 2 * g ** 2]]
    B[3:6, 3] = B[0:3, 3]
    B[0:3, 3] *= cos(pi / 2 - pi * g / 2)
    B[3:6, 3] *= sin(pi / 2 - pi * g / 2)

    B[0:3, 4] = [[1], [3 - 2 * g], [3 - 5 * g + 2 * g ** 2]]
    B[3:6, 4] = B[0:3, 4]
    B[0:3, 4] *= cos(pi / 2 - pi * g)
    B[3:6, 4] *= sin(pi / 2 - pi * g)

    B[0:3, 5] = [[1], [4 - 3 * g], [6 - 21 / 2 * g + 9 / 2 * g ** 2]]
    B[3:6, 5] = B[0:3, 5]
    B[0:3, 5] *= cos(pi / 2 - 3 / 2 * pi * g)
    B[3:6, 5] *= sin(pi / 2 - 3 / 2 * pi * g)

    return B


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
    g = symbols('gamma')    # Gamma
    A1, A2, A3, A4, A5, A6 = symbols('A1, A2, A3, A4, A5, A6')  # Coefficients

    # Assemble the linear system
    B = gen_B_mat(g)
    c_A = Matrix([A1, A2, A3, A4, A5, A6])
    RHS = B * c_A
    LHS = Matrix([1, 2, 1, 0, 0, 0])

    # # Select a portion of the linear system to solve
    # n = 2
    # B = B[0:n, 0:n]
    # c_A = c_A[0:n, 0]
    # RHS = B * c_A
    # LHS = LHS[0:n, 0]

    # Solve the linear system
    c_A_sol = solve(RHS - LHS, c_A)

    # Write the output solution file
    write_coef(c_A_sol, sol_file)
    return


if __name__ == "__main__":
    main()