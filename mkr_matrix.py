import numpy as np
from math import pow
import matplotlib.pyplot as plt

#матричная форма конечно-разностных соотношений;

#функция матрицы типа D1 размером n*(n+2);
def MATRIX_D1(n, step):
    
    ## коэфициенты матрицы и нужные для нее значения;
    h = 1 / (2 * step)
    
    ## создание матрицы D1;
    matrix = []
    for i in range(n):
        line = []
        for j in range(n+2):
            ## запись чисел в матрицу
            if i == j:
                x = -1 * h
            elif i+2 == j:
                x = 1 * h
            else:
                x = 0
            line.append(x)
        matrix.append(line)
    return matrix

#функция матрицы типа D2 размером n*(n+2);
def MATRIX_D2(n, step):
    ## коэфициенты матрицы и нужные для нее значения;
    h = 1 / (step ** 2)
    
    ## создание матрицы D2;
    matrix = []
    for i in range(n):
        line = []
        for j in range(n+2):
            ## запись чисел в матрицу
            if i == j or i + 2 == j:
                x = 1 * h
            elif i + 1 == j:
                x = -2 * h
            else:
                x = 0
            line.append(x)
        matrix.append(line)
    return matrix    

#функция матрицы типа D3 размером n*(n+2);
def MATRIX_D3(n, step):
    ## коэфициенты матрицы и нужные для нее значения;
    h = 1 / (2 * step)
    
    ## создание матрицы D3;
    matrix = []
    for i in range(n):
        line = []
        for j in range(n+2):
            ## запись чисел в матрицу
            if i == j:
                x = 1 * h
            elif i + 1 == j:
                x = -2 * h
            elif i + 3 == j:
                x = 2 * h
            elif i + 4 == j:
                x = -1 * h            
            else:
                x = 0
            line.append(x)
        matrix.append(line)
    return matrix  

#функция матрицы типа D4 размером n*(n+2);
def MATRIX_D4(n, step):
    ## коэфициенты матрицы и нужные для нее значения;
    h = 1 / (2 * step)
    
    ## создание матрицы D4;
    matrix = []
    for i in range(n):
        line = []
        for j in range(n+2):
            ## запись чисел в матрицу
            if i == j or i + 4 == j:
                x = 1 * h
            elif i + 1 == j or i + 3 == j:
                x = -4 * h
            elif i + 2 == j:
                x = 6 * h                     
            else:
                x = 0
            line.append(x)
        matrix.append(line)
    return matrix


'''# тут сосредоточенная продольная сила S в матричном виде;
def MATRIX_S(buffer_S, n):
    matrix = []
    for i in range(n):
        line = []
        for j in range(n):
            if i == j:
                x = buffer_S[j]
            else:
                x = 0
            line.append(x)
        matrix.append(line)
    return matrix

# тут расширенная S - S~. Запись такая -> [{0}[S]{0}]
def MATRIX_S_(buffer_S, n):
    matrix = []
    for i in range(n):
        line = []
        for j in range(n + 2):
            if i + 1 == j:
                x = buffer_S[i]
            else:
                x = 0
            line.append(x)
        matrix.append(line)
    return matrix
'''
# тут, исходя из записи EJ*D2*Y~ - S~*Y~ = M0, имеем матрицу A, которая состоит из (EJ*D2 - S~)
def MATRIX_A(mat_EJ, n, step):
    matrix_d2 = np.array(MATRIX_D2(n, step))
    matrix_res_d2_ej = np.dot(mat_EJ, matrix_d2) ##тут надо реализовать метод подгонки для возможной экономи
    matrix_A = matrix_res_d2_ej
    return matrix_A

#расширенная A - A~. запись такая -> [[010 ... 000], [A], [000 ... 010]]
def MATRIX_A_(mat_EJ, n, step):
    matrix_a = MATRIX_A(mat_EJ, n, step)
    ##создаем линии для расширения
    line1, line2 = np.zeros(n+2), np.zeros(n+2)
    line1[1], line2[-2] = 1, 1
    matrix_a = np.insert(matrix_a, 0, [line1], axis=0)
    matrix_a = np.append(matrix_a, [line2], axis=0)
    return matrix_a

#finally, деление A~ на M~
def MATRIX_Y_(mat_EJ, mat_m0_, n, step):
    mat_a_ = MATRIX_A_(mat_EJ, n, step)
    may_y_ = np.linalg.solve(mat_a_, mat_m0_)
    return may_y_

#проверка на сходимость
def Proverka(p, leng, quantity_0, quantity_end):
    cheker_list = []
    ej = 2 * (10 ** 11) * (0.005 * (0.01 ** 3) / 12)
    while(quantity_0 < quantity_end):
        step = leng / quantity_0
        x = np.linspace(0, leng, num=(quantity_0))
        
        m0 = np.zeros((quantity_0+2,1))
        buf_ej = np.eye(n, k=0) * ej
        
        m0[1:quantity_0//2 + 1, 0] = (p / 2) * x[0:quantity_0//2]
        m0[quantity_0//2 + 1:quantity_0+1, 0] = (p / 2) * (leng - x[quantity_0//2:quantity_0]) #M0
        res = MATRIX_Y_(buf_ej, m0, quantity_0, step)
        res_num = res[n//2,0]
        
        cheker_list = np.append(cheker_list, res_num)
        n += 2
        
    plt.plot(cheker_list, 'b--') #вывод на экран график
    plt.show()    
    return 0

'''w(z) = -Pl^3 / 48EJ * (3 z/l - 4 z^3/l^3) '''
#реализация графика функции W(z) (визи)
def WIZI(p, leng, quantity_0):
    ej = 2 * (10 ** 11) * (0.005 * (0.01 ** 3) / 12)
    x = np.linspace(0, leng, num=(quantity_0))
    const_1 = (-1)*(p / leng**3) / 48 * ej
    w_z = np.zeros((quantity_0, 1))
    w_z[0:quantity_0//2 + 1, 0] = const_1 * ( (3 * x[0:quantity_0//2 + 1] / leng) - (4 * ( x[0:quantity_0//2 + 1]**3 / leng**3 ) ) )
    w_z[quantity_0//2 + 1, 0] = const_1
    w_z[quantity_0//2 + 1:quantity_0, 0] = np.flip(w_z[0:quantity_0//2, 0], axis=0)
    
    plt.plot(w_z)
    plt.show()
    return 0
#================================================================================================================
## тут мы создаем информацию (набор данных)
p, l, n = 1000, 1, 11
h = l / (n - 1)
#================================================================================================================
#получение значения (EJy``)`` - (Ny`)`;
def _Y_(ej_buf, q_buf, N_buf, h, n):
    K_buf = np.zeros((n))
    k2, k1, k0 = [], [], []
    
    #поулучение частей для матрицы K
    for i in range(0,n+2,1):
        k2.append( ( ej_buf[i+1] + 2*ej_buf[i] + ej_buf[i-1] ) / (2 * pow(h,4)) )
        k1.append( ( (-6 * ej_buf[i]) / pow(h,4) ) + ( N_buf[i+1] + 4*N_buf[i] + N_buf[i-1] ) / (4 * pow(h,2)) )
        k0.append( ( 10*ej_buf[i] - 2*ej_buf[i+1] - 2*ej_buf[i-1] ) / pow(h,4) - ( 2*N_buf[i] / pow(h,2) ) )
    
    for i in range(n):
        line = []
        for j in range(n+2):
            if i == j or i+4 == j:
                x = k2[j]
            elif i+1 == j or i+3 == j:
                x = k1[j]
            elif i+2 == j:
                x = k0[j]
            else:
                x = 0
            line.append(x)
        K_buf[i] = line
    
    line1, line2 = np.zeros(n+2), np.zeros(n+2)
    line0, line3 = line1, line2
    
    line0[2] = 1
    line1[1:4,0] = [-1,0-1]
    line2[-2:-4] = [1,-2,1]
    line3[-1:-5] = [-1,2,0,-2,1]
    #получение расширенной матрицы K~
    K_buf = np.insert(K_buf, 0, [line0], axis=0)
    K_buf = np.insert(K_buf, 1, [line1], axis=0)
    K_buf = np.append(K_buf, [line2], axis=0)
    K_buf = np.append(K_buf, [line3], axis=0) 
    
    return 0


#================================================================================================================
