import numpy as np
import matplotlib.pyplot as plt

import read_point_x_y as read_p


# 构建系数矩阵函数：用于构建解N+2方程前n项，以及展示矩阵
def para_maker(work_para, work_x, row_n, input_x, columns_n):
    """ make parameter matrix ;
    make it front 4 columns as [1, x_i, (x_i)^2, (x_i)^3,...
    make it last columns as ...((x_i-(x_in)_(k-2))^3)_+]
    :param work_para: work space matrix
    :param work_x: x for put in matrix
    :param row_n: work space row number
    :param input_x: x for compared
    :param columns_n: work space columns number
    :return: matrix just made
    """
    for ii in range(row_n):
        start = 2
        for j in range(columns_n):

            if j < 4:
                work_para[ii, j] = work_x[ii] ** j
            else:
                if work_x[ii] > input_x[(start - 1)]:
                    work_para[ii, j] = (work_x[ii] - input_x[(start - 1)]) ** 3
                start = start + 1
    pass

    print('make a parameter matrix successfully')



# 自由端条件添加；在构建n+2方程组的前n个运行之后运行
def free_end(work_para_x, work_para_y, input_x, input_n):

    """add free_end condition to matrix
    using the second derivative of two end is 0;
    add two equation to the parameter matrix
    :param work_para_x: work space matrix of x
    :param work_para_y: work space vector of y
    :param input_x: coordinate value of x known
    :param input_n: number of coordinate value of x known
    :return: matrix added free_end condition
    """

    # 左端点处
    work_para_x[-2, 2] = 2.
    work_para_x[-2, 3] = 6. * x_in[0]
    work_para_y[-2, 0] = 0.
    # 右端点处
    work_para_x[(-1), 2] = 2.
    work_para_x[(-1), 3] = 6. * x_in[(-1)]
    work_para_y[(-1), 0] = 0.
    start = 2
    for j in range(4, input_n + 2):
        if input_x[(-1)] > input_x[(start - 1)]:
            work_para_x[(input_n + 1), j] = 6. * (input_x[(-1)] - input_x[(start - 1)])
        start += 1
    print('add free end condition successfully')


def fix_end(work_para_x, work_para_y, input_x, input_n, deg_1=None, deg_n=None, dy_1=None, dy_n=None):
    """add fixed_end condition to matrix
    using the first derivative of two end is fixed value known;
    add two equation to the parameter matrix
    :param work_para_x: work space matrix of x
    :param work_para_y: work space vector of y
    :param input_x: coordinate value of x known
    :param input_n: number of coordinate value of x known
    :param deg_1: the angle at point 1
    :param deg_n: the angle at point n
    :param dy_1: the first derivative of point 1
    :param dy_n: the first derivative of point n
    :return: matrix added fixed_end condition
    """
    if deg_1:
        dy_1 = np.tan(deg_1 * np.pi / 180)
        dy_n = np.tan(-deg_n * np.pi / 180)
    else:
        dy_1 = dy_1
        dy_n = dy_n
    # 左端点处
    work_para_x[-2, 1] = 1.
    work_para_x[-2, 2] = 2. * x_in[0]
    work_para_x[-2, 3] = 3. * x_in[0] ** 2
    work_para_y[-2, 0] = dy_1
    # 右端点处
    work_para_x[-2, 1] = 1.
    work_para_x[-2, 2] = 2. * x_in[-1]
    work_para_x[-2, 3] = 3. * x_in[-1] ** 2
    work_para_y[(-1), 0] = dy_n
    start = 2
    for j in range(4, input_n + 2):
        if input_x[(-1)] > input_x[(start - 1)]:
            work_para_x[(input_n + 1), j] = 3. * (input_x[(-1)] - input_x[(start - 1)]) ** 2
        start += 1
    print('add fixed end condition successfully')


def draw_curve(x_input, sol_a, columns_n):
    """ give show (x,y) using solved f(x)
    :param x_input: the array of x given
    :param sol_a: the array of a solved
    :param columns_n: the number of columns of show matrix
    :return: show(x,y)
    """
    show_x = np.arange(x_input[0], x_input[(-1)] + 1, 0.5)
    show_n = len(show_x)
    show_para = np.zeros((show_n, columns_n))
    # 根据坐标值对应得到坐标系数矩阵--调用para_maker
    para_maker(show_para, show_x, show_n, x_input, columns_n)
    show_y = show_para.dot(sol_a)
    return show_x, show_y


# --------------------解方程部分-----------------------------------------
def curve_fairing(x_input, y_input, deg_1=None, deg_n=None, dy_1=None, dy_n=None):
    """curve fairing
    :param x_input: coordinate value of x known
    :param y_input: coordinate value of y known
    :param deg_1: the angle at point 1
    :param deg_n: the angle at point n
    :param dy_1: the first derivative of point 1
    :param dy_n: the first derivative of point n
    :return:
    """
    # 创建方程系数矩阵
    N = len(x_input)
    para = np.zeros((N + 2, N + 2))  # x矩阵
    para_y = np.zeros((N + 2, 1))  # y向量

    # 构建x次幂系数矩阵前n行
    para_maker(para, x_input, N, x_input, N+2)

    # 构建y列向量前n行
    for i in range(N):
        para_y[i, 0] = y_input[i]

    # 根据端点条件再得到两个方程
    # 自由端条件y''=0
    # free_end(para, para_y, x_in, N)
    # 固定端条件
    fix_end(para, para_y, x_input, N, deg_1, deg_n, dy_1, dy_n)

    # -------求逆解方程
    para_inv = np.linalg.inv(para)
    sol_a = para_inv.dot(para_y)
    print('solve a successfully')

    plt.plot(x_input, y_input, 'ro')
    show_x, show_y = draw_curve(x_input, sol_a,  N + 2)
    plt.plot(show_x, show_y)
    plt.show()
    print('Finish fairing successfully')
    return 0


# ------------------------main------------------------------
if __name__ == '__main__':
    deg_1 = 60
    deg_n = 60

    x_in, y_in = read_p.read_txt('x_y.txt')
    print("x=", x_in)
    print("y=", y_in)

    curve_fairing(x_in, y_in, deg_1, deg_n)
