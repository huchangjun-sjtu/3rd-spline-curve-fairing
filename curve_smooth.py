import numpy as np
import matplotlib.pyplot as plt
import read_point_x_y as read_p
import spline_curve as sp


def transfer(input_x, input_y):
    """coordinate transformation
    :param input_x: coordinate value of x known
    :param input_y: coordinate value of y known
    :return: coordinate value of x,y for smooth method
    """
    theta = np.arctan2((input_y[-1] - input_y[0]), (input_x[-1] - input_x[0]))
    out_x = (input_x - input_x[0]) * np.cos(theta) + (input_y - input_y[0]) * np.sin(theta)
    out_y = (input_y - input_y[0]) * np.cos(theta) - (input_x - input_x[0]) * np.sin(theta)
    print('transfer successfully')
    return out_x, out_y


def para_dx_maker(work_x, row_n, input_x, columns_n):
    """ make parameter matrix dx;
    make it front 4 columns as [1, 2*x_i, 3*(x_i)^2,...
    make it last columns as ...(3*(x_i-(x_in)_(k-2))^2)_+]
    :param work_x: x for put in matrix
    :param row_n: work space row number
    :param input_x: x for compared
    :param columns_n: work space columns number
    :return: matrix just made
    """
    work_para = np.zeros((row_n, columns_n))
    for ii in range(row_n):
        start = 2
        for j in range(columns_n):
            if j < 3:
                work_para[ii, j] = work_x[ii] ** j * (j+1)
            else:
                if work_x[ii] > input_x[(start - 1)]:
                    work_para[ii, j] = (work_x[ii] - input_x[(start - 1)]) ** 2 * 3
                start = start + 1
    pass
    print('make a parameter matrix of dx successfully')
    return work_para


def para_ddx_maker(work_x, row_n, input_x, columns_n):
    """ make parameter matrix ddx;
    make it front 4 columns as [0, 2, 6*(x_i),...
    make it last columns as ...(6*(x_i-(x_in)_(k-2)))_+]
    :param work_x: x for put in matrix
    :param row_n: work space row number
    :param input_x: x for compared
    :param columns_n: work space columns number
    :return: matrix just made
    """
    work_para = np.zeros((row_n, columns_n))
    for ii in range(row_n):
        start = 2
        for j in range(columns_n):
            if j == 0:
                work_para[ii, j] = 0
            elif j == 1:
                work_para[ii, j] = 2
            elif j == 2:
                work_para[ii, j] = work_x[ii] * 6
            else:
                if work_x[ii] > input_x[(start - 1)]:
                    work_para[ii, j] = (work_x[ii] - input_x[(start - 1)]) * 6
                start = start + 1
    pass
    print('make a parameter matrix of ddx successfully')
    return work_para


def a_matrix_maker(matrix_x, input_y, row_n, columns_n, w):
    """make A parameter matrix for value difference at each point
    :param matrix_x: the matrix of x has been make in spline curve
    :param input_y: the array of y known
    :param row_n: number of row of A
    :param columns_n: number of column of A
    :param w: wight of each point
    :return:A parameter matrix ; an array of y about a
    """
    work_para = np.zeros((row_n, columns_n))
    work_y_a = np.zeros(row_n)
    matrix_wx = np.zeros(matrix_x.shape)
    array_wy = np.zeros(input_y.shape)
    for ii in range(len(matrix_x)):
        matrix_wx[ii, :] = w[ii] * matrix_x[ii, :]  # 对 X 矩阵 加权 w
        array_wy[ii] = w[ii] * input_y[ii]  # 对 Y 向量 加权 w
    matrix_sum = np.zeros(matrix_wx.shape)
    array_sum = np.zeros(input_y.shape)
    for ii in range(row_n):
        for i_sum in range(len(matrix_wx)):
            matrix_sum[i_sum, :] = matrix_wx[i_sum, :] * matrix_x[i_sum, ii + 1]  # 乘 f 对 a 的一阶导
            array_sum[i_sum] = array_wy[i_sum] * matrix_x[i_sum, ii + 1]  # 乘 f 对 a 的一阶导
        for j in range(columns_n):
            work_para[ii, j] = sum(matrix_sum[:, j + 1])
            work_y_a[ii] = sum(array_sum[:])
    pass
    print(work_y_a)
    print('make a parameter matrix of A and y successfully')
    return work_para, work_y_a


def b_matrix_maker(input_x, row_n, columns_n, q, dy):
    """make B parameter matrix for first derivative at each point
    :param input_x: the array of x known
    :param row_n: number of row of B parameter matrix
    :param columns_n: number of columns of B parameter matrix
    :param q: wight of first derivative of first and last point known
    :param dy: first derivative of first and last point known
    :return: B parameter matrix
    """
    work_para = np.zeros((row_n, columns_n))
    work_y_b = np.zeros(row_n)
    para_dx = para_dx_maker(input_x, row_n-1, input_x, columns_n)  # 调用一阶导系数矩阵
    para_dxn = para_dx[-1, :]  # 使用 Xn 的一阶导系数
    for i in range(row_n):
        for j in range(columns_n):
            work_para[i, j] = para_dxn[i] * para_dxn[j] * q[-1]
        work_y_b[i] = q[-1] * para_dxn[i] * dy[-1]
    work_para[0, 0] += q[0]
    work_y_b[0] += q[0] * dy[0]
    print(work_y_b)
    print('make a parameter matrix of B successfully')
    return work_para, work_y_b


def s_matrix_maker(row_n, columns_n, s):
    """make S parameter matrix for smooth at each point
    :param row_n: number of row of S parameter matrix
    :param columns_n: number of columns of S parameter matrix
    :param s: the parameter of smooth level given and adjust until smooth enough
    :return: S parameter matrix
    """
    s_array = np.zeros(row_n)
    for i in range(3, columns_n):
        s_array[i] = s
    work_para = np.diag(s_array)
    print('make a parameter matrix of S successfully')
    return work_para


def check_smooth(input_x, sol_a, row_n, columns_n):
    dx = para_dx_maker(input_x, row_n, input_x, columns_n)
    ddx = para_ddx_maker(input_x, row_n, input_x, columns_n)
    dy = dx.dot(sol_a)
    print('dy =', dy)
    ddy = ddx.dot(sol_a)
    print('ddy =', ddy)
    k = ddy / ((1 + dy ** 2) ** (3 / 2))
    print('k = ', k)
    if (ddy[0] * ddy[1] < 0) and (ddy[1] * ddy[2] < 0):
        print('f1')
        return False
    for i in range(2, len(input_x)-1):
        if (ddy[i-1] * ddy[i] < 0) and (ddy[i] * ddy[i+1] < 0):
            print('f2-', i)
            return False
        elif (ddy[i-2] * ddy[i-1] < 0) and (ddy[i-1] * ddy[i] > 0) and(ddy[i] * ddy[i+1] < 0):
            print('f3', i)
            return False
        elif (k[i]**2 < min([k[i-1] ** 2, k[i+1] ** 2])) or (k[i]**2 > max([k[i-1] ** 2, k[i+1] ** 2])):
            print('f4', i)
            return False
    return True


def solve_curve_smooth(work_x, work_y, input_w, input_q, dy):
    s_in = 0
    sol_a = np.zeros(len(work_x))
    check = False
    n = len(work_x)
    para_x = np.zeros((n, n + 2))  # x矩阵
    sp.para_maker(para_x, work_x, n, work_x, n + 2)
    print(para_x)
    # A型值偏离矩阵
    para_a, array_a = a_matrix_maker(para_x, work_y, n + 1, n + 1, input_w)
    # B斜率偏离矩阵
    para_b, array_b = b_matrix_maker(work_x, n + 1, n + 1, input_q, dy)
    # 判断光顺进入循环
    while check is False:
        # S光顺系数矩阵
        para_s = s_matrix_maker(n + 1, n + 1, s_in)
        # 解方程部分：
        para_all = para_a + para_b + para_s
        array_all = array_a + array_b
        para_inv = np.linalg.inv(para_all)
        sol_a = para_inv.dot(array_all)
        # 判别光顺：
        check = check_smooth(work_x, sol_a, n, n + 1)
        # check = True
        s_in += 1
        if s_in > 50:
            break
    print('Finish spline curve sueccessfully')
    return sol_a, s_in


# ------------------main--------------------
# 由文件读入坐标数据
# x_in, y_in = read_p.read_txt('x_y.txt')
'''
x_in = np.array([1., 2., 4., 7., 8., 9.], dtype=float)
y_in = np.array([2., 2.5, 3., 2.5, 2., 2.], dtype=float)
'''
q_in = np.array([1., 1.], dtype=float)
dy_in = np.array([0.95, 0], dtype=float)
# dy_in = np.array([0.95, -0.4], dtype=float)

# 对坐标进行坐标转换得到工作坐标
# [tran_x, tran_y] = transfer(x_in, y_in)
tran_x = np.array([0., 2., 4., 7., 8., 9.], dtype=float)
tran_y = np.array([0., 1.5, 2., 2., 2., 2.], dtype=float)
w_in = np.ones(len(tran_x))
print(tran_x)
print(tran_y)
n_in = len(tran_x)
# 解方程
solve_a, s_now = solve_curve_smooth(tran_x, tran_y, w_in, q_in, dy_in)
#  输出图像
show_x = np.arange(tran_x[0], tran_x[(-1)] + 1, 0.1)
show_n = len(show_x)
show_para = np.zeros((show_n, n_in + 2))
# 根据坐标值对应得到坐标系数矩阵--调用para_maker
sp.para_maker(show_para, show_x, show_n, tran_x, n_in + 2)
show_a = np.hstack(([0], solve_a))  # 此处表达式与spline不同，需要在前面加0元素
show_y = show_para.dot(show_a)
print('s_in = ', s_now)
print('sol_a ', solve_a)
plt.plot(tran_x, tran_y, 'ro')
plt.plot(show_x, show_y)
plt.show()
