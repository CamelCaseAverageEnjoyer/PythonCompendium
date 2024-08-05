import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------- ТРИВИАЛЬНЫЕ ФУНКЦИИ
def clip(a: float, bot: float, top: float) -> float:
    if a < bot:
        return bot
    if a > top:
        return top
    return a


# ----------------------------------------- ФУНКЦИИ МАТРИЦ И ВЕКТОРОВ
def quart2dcm(L):
    """Функция ищет матрицу поворота из кватерниона поворота; \n
    Кватернион L передаётся вектором длины 4; \n
    Возвращает матрицу 3х3."""
    w, x, y, z = L
    A = np.eye(3)
    A[0][0] = 1 - 2 * y ** 2 - 2 * z ** 2
    A[0][1] = 2 * x * y + 2 * z * w
    A[0][2] = 2 * x * z - 2 * y * w
    A[1][0] = 2 * x * y - 2 * z * w
    A[1][1] = 1 - 2 * x ** 2 - 2 * z ** 2
    A[1][2] = 2 * y * z + 2 * x * w
    A[2][0] = 2 * x * z + 2 * y * w
    A[2][1] = 2 * y * z - 2 * x * w
    A[2][2] = 1 - 2 * x ** 2 - 2 * y ** 2
    return A

def q_dot(L1, L2):
    """Функция является кватернионным умножением; \n
    Кватернион L1,L2 передаются векторами длины 4; \n
    Возвращает кватернион L[0]..L[3]."""
    return np.array([L1[0] * L2[0] - L1[1] * L2[1] - L1[2] * L2[2] - L1[3] * L2[3],
                     L1[0] * L2[1] + L1[1] * L2[0] + L1[2] * L2[3] - L1[3] * L2[2],
                     L1[0] * L2[2] + L1[2] * L2[0] + L1[3] * L2[1] - L1[1] * L2[3],
                     L1[0] * L2[3] + L1[3] * L2[0] + L1[1] * L2[2] - L1[2] * L2[1]])

def get_U_S(A, w_0, t, ECCENTRICITY = None, INCLINATION = None):
    from sympy import sin, cos, Matrix, atan, sqrt, tan
    if ECCENTRICITY is None and INCLINATION is None:
        U = Matrix([[0, 1, 0],
                    [0, 0, 1],
                    [1, 0, 0]]) @ \
            Matrix([[cos(t * w_0), sin(t * w_0), 0],
                    [-sin(t * w_0), cos(t * w_0), 0],
                    [0, 0, 1]])
    else:
        e = 0 if ECCENTRICITY is None else ECCENTRICITY
        i = 0 if INCLINATION is None else INCLINATION
        E = t * w_0  # Эксцентрическая аномалия
        f = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2))  # Истинная аномалия
        U = Matrix([[0, 1, 0],  # Поворот к экваториальной плоскости
                    [0, 0, 1],
                    [1, 0, 0]]) @ \
            Matrix([[cos(f), sin(f), 0],  # Разница между истинной аномалией и местной
                    [-sin(f), cos(f), 0],
                    [0, 0, 1]]) @ \
            Matrix([[1, 0, 0],
                    [0, cos(i), sin(i)],  # Поворот к плоскости орбиты
                    [0, -sin(i), cos(i)]])
    S = A @ U.T
    return U, S

# ----------------------------------------- ПЕРЕХОДЫ МЕЖДУ СИСТЕМАМИ КООРДИНАТ
def i_o(self, a, U=None):
    """Инерциальная -> Орбитальная"""
    a_np = np.array(a)
    U = self.U if (U is None) else U
    return U @ a_np - np.array([0, 0, self.Re])

def o_i(self, a, U=None):
    """Орбитальная -> Инерциальная"""
    a_np = np.array(a)
    U = self.U if (U is None) else U
    return U.T @ (a_np + np.array([0, 0, self.Re]))
