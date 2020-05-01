import numpy as np
import matplotlib.pyplot as plt


# 构建系数矩阵函数：用于构建解N+2方程前n项，以及展示矩阵
def para_maker(work_para, work_x, work_n, input_x, input_n):
    for ii in range(work_n):
        start = 2
        for j in range(input_n + 2):
            if j < 4:
                work_para[ii, j] = work_x[ii] ** j
            else:
                if work_x[ii] > input_x[(start - 1)]:
                    work_para[ii, j] = (work_x[ii] - input_x[(start - 1)]) ** 3
                start = start + 1
    pass


# 自由端条件添加；在构建n+2方程组的前n个运行之后运行
def free_end(work_para_x, work_para_y, input_x, input_n):
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


# -------------------输入坐标部分，之后改进为读入文件或输入数据
#'''
x_in = np.array([-6919, 0, 7312.5,	9750, 14625, 19500,	24375,	29250,	39000,	48750,
              58500,	68250,	78000, 87750, 97500, 107250,	117000,	126750,	136500,
              146250,	156000,	165750,	175500,	185250,	190125,	195000], dtype=float)
y_in = np.array([602,	7220,	10438,	11282,	12687,	13837,	14663,	15344,	16192,	16400,
              16400, 16400, 16400, 16400, 16400, 16392, 16389, 16392, 16383,
              15784, 14322, 11436, 7652, 3715, 2035, 14], dtype=float)
'''
x_in = np.array([1, 3, 5, 7, 9], dtype=float)
y_in = np.array([1, 3, 0, -2, 6], dtype=float)
 '''

# --------------------解方程部分-----------------------------------------
# 创建方程系数矩阵
N = len(x_in)
para = np.zeros((N + 2, N + 2))  # x矩阵
para_y = np.zeros((N + 2, 1))  # y向量

# 构建x次幂系数矩阵前n行
para_maker(para, x_in, N, x_in, N)

# 构建y列向量前n行
for i in range(N):
    para_y[i, 0] = y_in[i]

# 根据端点条件再得到两个方程
# 自由端条件y''=0
free_end(para, para_y, x_in, N)


# -------求逆解方程
para_inv = np.linalg.inv(para)
sol_a = para_inv.dot(para_y)

# -------画出拟合图像
show_x = np.arange(x_in[0], x_in[(-1)] + 1, 0.5)
show_N = len(show_x)
show_para = np.zeros((show_N, N + 2))
# 根据坐标值对应得到坐标系数矩阵--调用para_maker
para_maker(show_para, show_x, show_N, x_in, N)

'''
for i in range(show_N):
    start = 2
    for j in range(N + 2):
        if j < 4:
            show_para[i, j] = show_x[i] ** j
        else:
            if show_x[i] > x_in[(start - 1)]:
                show_para[i, j] = (show_x[i] - x_in[(start - 1)]) ** 3
            start = start + 1
'''
show_y = show_para.dot(sol_a)

print(para)
print(sol_a)
print(show_para)
print(para_y)
print(show_y)

x_range = max(x_in) - min(x_in)
y_range = max(y_in) - min(y_in)
# plt.figure(figsize=(10, 10*y_range/x_range), dpi=80)
plt.plot(x_in, y_in, 'ro')
plt.plot(show_x, show_y)

plt.show()
