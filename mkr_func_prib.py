from math import pow
##описываем приближенные значения производных с первого по 4 порядок.
##ниже описываемые функции используются в цикле или отдельно

#получение приближенного значения производной первого порядка
def PORYADOK_1(step, func_left_1, func_right_1): 
    return (func_right_1 - func_left_1) / (2 * step)

#получение приближенного значения производной второго порядка
def PORYADOK_2(step, func_left_1, func_right_1, func_centr): 
    return (func_right_1 - 2 * func_centr - func_left_1) / pow(step, 2)

#получение приближенного значения производной третьего порядка
def PORYADOK_3(step, func_left_1, func_right_1, func_right_2, func_left_2):
    return (func_left_2 - 2 * func_left_1 + 2 * func_right_1 - func_right_2) / (2 * pow(step, 3))

#получение приближенного значения производной четвертого порядка
def PORYADOK_4(step, func_left_1, func_right_1, func_right_2, func_left_2, func_centr):
    return (func_left_2 - 4 * func_left_1 + 6 * func_centr - 4 * func_right_1 + func_right_2) / pow(step, 4)


