import numpy as np


# 判断字符是否为数字
def is_number(s):
    try:
        if (s is '-') or (s is '.'):
            return True
    except ValueError:
        pass

    try:
        float(s)
        return True
    except ValueError:
        pass

    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass

    return False


# 从txt文件读入坐标
def read_txt(work_file):
    with open(work_file) as f:
        read_data = f.read()
    x_y = []
    string = ''
    for str_ in read_data:
        if is_number(str_):
            string += str_
        else:
            if string:
                x_y.append(string)
                string = ''
    input_x = np.zeros(int(len(x_y) / 2))
    input_y = np.zeros(int(len(x_y) / 2))
    for i in range(0, len(x_y), 2):
        input_x[int(i / 2)] = float(x_y[i])
        input_y[int(i / 2)] = float(x_y[i + 1])
    print('read points txt successfully')
    return input_x, input_y

