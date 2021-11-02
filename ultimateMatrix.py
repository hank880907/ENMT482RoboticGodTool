"""
file: ultimateMatrix.py
Author: Hank Wu 30 Oct 2021
"""
from math import pi, floor, ceil
from sympy import Matrix, Number, factor, trigsimp, simplify, sympify, nsimplify, pretty, Eq, nsolve, eye
from typing import List



def round_expr(expr, num_digits):
    return expr.xreplace({n : round(n, num_digits) for n in expr.atoms(Number)})

def matrix_printf(m: Matrix, title = "Matrix", padding = 3, pre_padding = "", print_by_row_thresh = 150, precision = 4) -> str:
    print_by_row = False
    #get the dimention.
    num_rows, num_cols = m.shape
    #get the list of strings
    m_strings = []
    for i in range(num_rows):
        for j in range(num_cols):
            k = simplify(m[i,j])
            k = factor(k)
            k = trigsimp(k)
            k = nsimplify(k, rational = True)
            k = round_expr(k, precision)
            k = pretty(k)
            m_strings.append(k)

    # calculate the width for each column and decide to print by row or not.
    col_widths = [0] * num_cols
    for i in range(num_rows):
        for j in range(num_cols):
            string = m_strings[i*num_cols + j]
            if len(string) > col_widths[j]:
                col_widths[j] = len(string) + padding

    #check if the lenth is exceed print by row thresh
    total_len = sum(col_widths)
    if total_len > print_by_row_thresh:
        print_by_row = True

    # generate the string
    if print_by_row:
        ret_string = "\n"
        ret_string += "_" * floor(max(col_widths)/2 - len(title)) + title + "_" * ceil(max(col_widths)/2 - padding)
        ret_string += "\n"
        for i in range(num_rows):
            ret_string += "\n"
            for j in range(num_cols):
                this_str = m_strings[i*num_cols + j]
                ret_string += f"row{i}, col{j}: {this_str}" + "\n"
        ret_string += "_" * (max(col_widths)-padding)
        ret_string += "\n"
    else:
        ret_string = "\n"
        ret_string += pre_padding + "_" * floor(total_len/2 - len(title)) + title + "_" * ceil(total_len/2 - padding)
        ret_string += "\n"
        for i in range(num_rows):
            ret_string += "\n"
            for j in range(num_cols):
                this_str = m_strings[i*num_cols + j]
                ret_string += pre_padding + f"{this_str}" + " "*(col_widths[j] - len(this_str))
        ret_string += "\n"
        ret_string += pre_padding + "_" * (total_len - padding)
        ret_string += "\n"
    
    return ret_string

def matrix_printf3x1(m: Matrix) -> str:
    mystr = matrix_printf(m, title="vector", pre_padding="\t\t\t")
    return mystr


class tMatrix():
    """transformation matrix generator"""
    def __init__(self, theta, d, alpha, a, print_precision = 3) -> None:
        self.precision = print_precision
        self.my_symbles = []
        theta_mat: Matrix = self.get_theta_mat(theta,d)
        alpha_mat: Matrix = self.get_alpha_mat(alpha, a)
        self.my_T: Matrix = theta_mat * alpha_mat
        pass

    def __str__(self) -> str:
        return matrix_printf(self.my_T, precision = self.precision)

    def __mul__(self, other: 'tMatrix') -> 'tMatrix':
        retval = tMatrix(0,0,0,0)
        retval.my_T = self.my_T * other.my_T
        retval.my_symbles = self.my_symbles + other.my_symbles
        return retval
    
    def get_theta_mat(self, theta, d) -> Matrix:
        if isinstance(theta, str):
            theta_r = theta
        else:
            theta_r = theta * pi / 180

        e00 = sympify(f"cos({theta_r})")
        for i in e00.free_symbols:
            self.my_symbles.append(i)
        e01 = sympify(f"-sin({theta_r})")
        e10 = sympify(f"sin({theta_r})")
        e11 = sympify(f"cos({theta_r})")
        e02 = sympify("0")
        e03 = sympify("0")
        e12 = sympify("0")
        e13 = sympify("0")
        e20 = sympify("0")
        e21 = sympify("0")
        e22 = sympify("1")
        e23 = sympify(f"{d}")
        for i in e23.free_symbols:
            self.my_symbles.append(i)
        e30 = sympify("0")
        e31 = sympify("0")
        e32 = sympify("0")
        e33 = sympify("1")
        my_mat = [[e00,e01,e02,e03],[e10,e11,e12,e13],[e20,e21,e22,e23],[e30,e31,e32,e33]]
        return Matrix(my_mat)

    def get_alpha_mat(self,alpha,a):
        if isinstance(alpha, str):
            alpha_r = alpha
        else:
            alpha_r = alpha * pi / 180
        e11 = sympify(f"cos({alpha_r})")
        for i in e11.free_symbols:
            self.my_symbles.append(i)
        e12 = sympify(f"-sin({alpha_r})")
        e21 = sympify(f"sin({alpha_r})")
        e22 = sympify(f"cos({alpha_r})")


        e00 = sympify("1")
        e01 = sympify("0")
        e02 = sympify("0")
        
        e03 = sympify(f"{a}")
        for i in e11.free_symbols:
            self.my_symbles.append(i)
        e10 = sympify("0")
        
        e13 = sympify("0")
        e20 = sympify("0")

        e23 = sympify("0")
        e30 = sympify("0")
        e31 = sympify("0")
        e32 = sympify("0")
        e33 = sympify("1")
        test = [[e00,e01,e02,e03],[e10,e11,e12,e13],[e20,e21,e22,e23],[e30,e31,e32,e33]]
        return Matrix(test)

    def ik(self, target) -> dict:
        eqs = []
        # extract non-trivial equations.
        for i in range(4):
            for j in range(4):
                if self.my_T[i, j].is_real:
                    pass
                else:
                    eq = Eq(self.my_T[i, j], target[i][j])
                    eqs.append(eq)

        num_syms = len(self.my_symbles)
        init_guess = []
        for i in range(num_syms):
            init_guess.append(1)
        results = nsolve(eqs, self.my_symbles, init_guess, verify = False)
        retval = {}
        for i, result in enumerate(results):
            retval[self.my_symbles[i]] = result
        return retval


def jacobian_god_function(t_matrices_list: List[tMatrix], is_prismatic_list: List[bool], spacing = 30, verbose = True) -> tMatrix:
    T_n_from_0 = eye(4)
    for i in range(len(t_matrices_list)):
        T_n_from_0 = T_n_from_0 * t_matrices_list[i].my_T

    P_n_from_0 = Matrix([T_n_from_0.row(0).col(3), T_n_from_0.row(1).col(3), T_n_from_0.row(2).col(3)])

    cols = []

    for i, t_matrix_n_from_0 in enumerate(t_matrices_list[:-1]):
        T_i_from_0 = eye(4)
        for j in range(i+1):
            T_i_from_0 = T_i_from_0 * t_matrices_list[j].my_T

        Z_i_from_0 = Matrix([T_i_from_0.row(0).col(2), T_i_from_0.row(1).col(2), T_i_from_0.row(2).col(2)])

        P_i_from_0 = Matrix([T_i_from_0.row(0).col(3), T_i_from_0.row(1).col(3), T_i_from_0.row(2).col(3)])
        r_i_from_0 = P_n_from_0 - P_i_from_0

        z_cross_r = Z_i_from_0.cross(r_i_from_0)

        if is_prismatic_list[i]:
            col_i = [Z_i_from_0.row(0).col(0), Z_i_from_0.row(1).col(0), Z_i_from_0.row(2).col(0), 0, 0, 0]

        else:
            col_i = [z_cross_r.row(0), z_cross_r.row(1), z_cross_r.row(2), Z_i_from_0.row(0), Z_i_from_0.row(1), Z_i_from_0.row(2)]

        if verbose:  # for debug
            hank = tMatrix(0,0,0,0)
            hank.my_T = T_i_from_0
            print('i = ', i, '_' * 90)
            print(f'T{i}   =\n{hank}')
            print(f'Z{i}   =', matrix_printf3x1(Z_i_from_0))
            print(f'P{i}   =', matrix_printf3x1(P_i_from_0))
            print(f'r{i}   =', matrix_printf3x1(r_i_from_0))
            print(f'zxr{i} =', matrix_printf3x1(z_cross_r))
            print(f'col{i} =', matrix_printf3x1(Matrix(col_i)))
            print('_' * 100)


        cols.append(col_i)

    jacobian = Matrix(cols).T
    jacobian = nsimplify(jacobian, rational = True)
    jacobian = round_expr(jacobian, 4)
    jacobian = factor(jacobian)
    retval = tMatrix(0,0,0,0)
    retval.my_T = jacobian
    return retval







if __name__ == "__main__":

    # t1 = tMatrix("theta1", 1, 90, 0)
    # t2 = tMatrix("theta2", 0,0,1)
    # t3 = tMatrix("theta3", 0,0,0)
    # t4 = tMatrix(0,0,0,"L3")
    # overall = t1*t2*t3*t4

    # testm = overall.my_T
    # print(overall)
    t1 = tMatrix("theta1", 1, 90 ,0)
    t2 = tMatrix("theta2", 0, 0 , 1)
    t3 = tMatrix("theta3", 0, 0, 0)
    t4 = tMatrix(0, 0, 0, 0.5)
    fk: tMatrix = t1*t2*t3*t4

    target = [
    [0.612,     -0.354,     0.707 ,     0.806],     
    [0.612,     -0.354,     -0.707,     0.806],     
    [0.500,     0.866 ,     0.0   ,     1.957],     
    [0    ,     0     ,     0     ,     1    ]]

    ik = fk.ik(target)
    for key in ik:
        if f"{key}".find("theta") == -1:
            print(f"{key} = {ik[key]:.2f}")
        else:
            print(f"{key} = {ik[key]*180/pi:.2f}")

