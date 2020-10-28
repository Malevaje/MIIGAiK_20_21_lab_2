# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 09:07:04 2020

@author: Oleg_Chekin
"""
from Linking_to_two_wall_markers.func_with_angles import *

position_of_the_plumb = {
    'I': {
        'a, m': 5.1170,
        'b, m': 9.4250,
        'c, m': 4.3090,
        'α': [1, 40, 23],
        'ω': [127, 47, 55],
        'a_1, m': 5.1200,
        'b_1, m': 5.3930,
        'c_1, m': 10.5110,
        'α_1': [0, 52, 21],
        'ω_1': [173, 7, 7],
    },
    'II': {
        'a, m': 5.1210,
        'b, m': 9.4250,
        'c, m': 4.3080,
        'α': [1, 52, 20],
        'ω': [127, 35, 57],
        'a_1, m': 5.1210,
        'b_1, m': 5.3930,
        'c_1, m': 10.5110,
        'α_1': [0, 47, 25],
        'ω_1': [173, 2, 14],
    },
    'III': {
        'a, m': 5.1170,
        'b, m': 9.4250,
        'c, m': 4.3090,
        'α': [1, 28, 25],
        'ω': [127, 59, 52],
        'a_1, m': 5.1180,
        'b_1, m': 5.3930,
        'c_1, m': 10.5090,
        'α_1': [0, 57, 13],
        'ω_1': [173, 11, 59],
    },
}
source_data = {
    'α_T-PZ20': [125, 57, 1],
    'X_PZ20': 669.962,
    'Y_PZ20': 1044.082,
}
designations_position = [
    'a, m',
    'b, m',
    'c, m',
    'α',
    'ω',
    'a_1, m',
    'b_1, m',
    'c_1, m',
    'α_1',
    'ω_1',
]
designations = [
    'α',
    'γ',
    'β',
    'Невязка',
    'α_1',
    'β_1',
    'γ_1',
    'Невязка_1',
]
designations_t5 = [
    'd_1',
    'd_2',
    'С_выч',
    'f_s',
    'd_1_1',
    'd_1_2',
    'С_выч_1',
    'f_s_1',
]
designations_t6 = [
    '(a)',
    '(b)',
    '(c)',
    '(a_1)',
    '(b_1)',
    '(c_1)',
]
designations_t7 = [
    'a',
    'b',
    'c',
    'a_1',
    'b_1',
    'c_1',
]
designations_t8 = [
    'α',
    'γ',
    'β',
    'Сумма',
    'α_1',
    'γ_1',
    'β_1',
    'Сумма',
]


def measure_log():
    print('  Таблица 1. Журнал измерений')
    print('| {:^10} | {:^51} |'.format('', 'Положение отвеса'))
    print('| {:^10} | {:_^51} |'.format('Обозна-', ''))
    print('| {:^10} | {:^15} | {:^15} | {:^15} |'.format('чение',
                                                         'I',
                                                         'II',
                                                         'III'))
    print('| {:-^10} | {:-^15} | {:-^15} | {:-^15} |'.format('', '', '', ''))
    a = list(position_of_the_plumb.values())
    log = []
    for i in a:
        log.append(list(i.values()))

    for i in range(len(log)):
        for k in range(len(log[i])):
            if type(log[i][k]) == list:
                log[i][k] = sek_r_grad_str(grad_r_sek_list(log[i][k]))

    for i in range(len(designations_position)):
        print('| {:^10} | {:^15} | {:^15} | {:^15} |'.format(designations_position[i],
                                                             log[0][i],
                                                             log[1][i],
                                                             log[2][i],
                                                             )
              )
    print()


def omega_check():
    """
    Выполняет проверку условия: теоретическое значение угла ω не должно превышать
    фактическое более чем на 10" - для надземных измерений, 15" - для подземных.
    :return: True - если измерения не прошли проверку, иначе - False
    """
    print('  Таблица 2. Проверка ω')
    print('| {:-^73} |'.format(''))
    print('| {:^15} | {:^15} | {:^16} | {:^18} |'.format('|ω_I - ω_II|',
                                                         '|ω_I - ω_III|',
                                                         '(15 * ro) / c',
                                                         'measurement'))
    print('| {:^15} | {:^15} | {:^16} | {:^18} |'.format('',
                                                         '',
                                                         '',
                                                         'tolerance'))
    print('| {:-^15} | {:-^15} | {:-^16} | {:-^18} |'.format('', '', '', ''))
    # on the surface
    a = abs(grad_r_sek_list(position_of_the_plumb['II']['ω'])
            - grad_r_sek_list(position_of_the_plumb['I']['ω']))
    b = abs(grad_r_sek_list(position_of_the_plumb['I']['ω'])
            - grad_r_sek_list(position_of_the_plumb['III']['ω']))
    c = (15 * constants['ro']) / (position_of_the_plumb['I']['c, m'] * 10 ** 3)
    d = (0.1 * constants['ro']) / (position_of_the_plumb['I']['c, m'] * 10 ** 3)
    print('| {:^15} | {:^15} | {:^16} | {:^18} |'.format(sek_r_grad_str(a),
                                                         sek_r_grad_str(b),
                                                         sek_r_grad_str(abs(c)),
                                                         round(d, 3)))

    # underground
    a_1 = abs(grad_r_sek_list(position_of_the_plumb['II']['ω_1'])
              - grad_r_sek_list(position_of_the_plumb['I']['ω_1']))
    b_1 = abs(grad_r_sek_list(position_of_the_plumb['I']['ω_1'])
              - grad_r_sek_list(position_of_the_plumb['III']['ω_1']))
    c_1 = (15 * constants['ro']) / (position_of_the_plumb['I']['c_1, m'] * 10 ** 3)
    d_1 = (0.1 * constants['ro']) / (position_of_the_plumb['I']['c_1, m'] * 10 ** 3)
    print('| {:^15} | {:^15} | {:^16} | {:^18} |'.format(sek_r_grad_str(a_1),
                                                         sek_r_grad_str(b_1),
                                                         sek_r_grad_str(c_1),
                                                         round(d_1, 3)))
    print()
    if abs(a - c) > 10 or abs(a_1 - c_1) > 15:
        print('Неудовлетворительные результаты измерений!')
        return True
    return False


def alpha_check():
    """
    Выполняет проверку условия: теоретическое значение угла α не должно превышать
    фактическое более чем на 10" - для надземных измерений, 15" - для подземных.
    :return: True - если измерения не прошли проверку, иначе - False
    """
    print('  Таблица 3. Проверка α')
    print('| {:-^73} |'.format(''))
    print('| {:^15} | {:^15} | {:^16} | {:^18} |'.format('|α_I - α_II|',
                                                         '|α_I - α_III|',
                                                         '(15 * ro) / c - ',
                                                         'measurement'))
    print('| {:^15} | {:^15} | {:^16} | {:^18} |'.format('',
                                                         '',
                                                         '(15 * ro) / b',
                                                         'tolerance'))
    print('| {:-^15} | {:-^15} | {:-^16} | {:-^18} |'.format('', '', '', ''))
    # on the surface
    a = abs(grad_r_sek_list(position_of_the_plumb['II']['α'])
            - grad_r_sek_list(position_of_the_plumb['I']['α']))
    b = abs(grad_r_sek_list(position_of_the_plumb['I']['α'])
            - grad_r_sek_list(position_of_the_plumb['III']['α']))
    c = abs(((15 * constants['ro']) / (position_of_the_plumb['I']['c, m'] * 10 ** 3) -
             (15 * constants['ro']) / (position_of_the_plumb['I']['b, m'] * 10 ** 3)))
    d = 0.14 * np.sqrt((constants['ro'] / (position_of_the_plumb['I']['c, m'] * 10 ** 3)) ** 2 +
                       (constants['ro'] / (position_of_the_plumb['I']['b, m'] * 10 ** 3)) ** 2)
    print('| {:^15} | {:^15} | {:^16} | {:^18} |'.format(sek_r_grad_str(a),
                                                         sek_r_grad_str(b),
                                                         sek_r_grad_str(c),
                                                         round(d, 3)))

    # underground
    a_1 = abs(grad_r_sek_list(position_of_the_plumb['II']['α_1'])
              - grad_r_sek_list(position_of_the_plumb['I']['α_1']))
    b_1 = abs(grad_r_sek_list(position_of_the_plumb['I']['α_1'])
              - grad_r_sek_list(position_of_the_plumb['III']['α_1']))
    c_1 = abs(((15 * constants['ro']) / (position_of_the_plumb['I']['c_1, m'] * 10 ** 3) -
               (15 * constants['ro']) / (position_of_the_plumb['I']['b_1, m'] * 10 ** 3)))
    d_1 = 0.14 * np.sqrt((constants['ro'] / (position_of_the_plumb['I']['c_1, m'] * 10 ** 3)) ** 2 +
                         (constants['ro'] / (position_of_the_plumb['I']['b_1, m'] * 10 ** 3)) ** 2)
    print('| {:^15} | {:^15} | {:^16} | {:^18} |'.format(sek_r_grad_str(a_1),
                                                         sek_r_grad_str(b_1),
                                                         sek_r_grad_str(c_1),
                                                         round(d_1, 3)))
    if abs(a - c) > 10 or abs(a_1 - c_1) > 15:
        print('Неудовлетворительные результаты измерений!')
        return True
    return False


def adjustment_a():
    """
    Вычисляет углы треугольника бета и гамма, невязки.
    Формирует таблицу 4. Вычисления углов и невязки.
    :return: f_s
    """
    f = position_of_the_plumb
    a = []
    c = []
    # Добавляем для печати углы алфа
    for i in f:
        c.append(grad_r_sek_list(f[i]['α']))
    a.append(c)

    c = []
    # добавляем для печати углы бета
    for i in f:
        cos_beta = ((f[i]['a, m'] ** 2 + f[i]['c, m'] ** 2 - f[i]['b, m'] ** 2) /
                    (2 * f[i]['a, m'] * f[i]['c, m']))
        if cos_beta < 0:
            beta = (gradVrad(constants['180']) - np.arcsin(
                    (f[i]['b, m'] * np.sin(gradVrad(grad_r_sek_list(f[i]['α'])))) /
                    f[i]['a, m'])
                    )
        else:
            beta = np.arcsin(f[i]['b, m'] * np.sin(gradVrad(grad_r_sek_list(f[i]['α']))) /
                             f[i]['a, m'])
        c.append(radVgrad(beta))
    a.append(c)

    c = []
    # Добавляем для печати углы гамма
    for i in f:
        cos_gamma = ((f[i]['a, m'] ** 2 + f[i]['b, m'] ** 2 - f[i]['c, m'] ** 2) /
                     (2 * f[i]['a, m'] * f[i]['b, m']))
        if cos_gamma < 0:
            gamma = (gradVrad(constants['180']) - np.arcsin(
                (f[i]['c, m'] * np.sin(gradVrad(grad_r_sek_list(f[i]['α'])))) /
                f[i]['a, m'])
                     )
        else:
            gamma = np.arcsin(f[i]['c, m'] * np.sin(gradVrad(grad_r_sek_list(f[i]['α']))) /
                              f[i]['a, m'])
        c.append(radVgrad(gamma))
    a.append(c)

    # Считаем невязку в треугольнике на дневной поверхности
    c = []
    j = 0
    while j < 3:
        n = a[0][j] + a[1][j] + a[2][j] - constants['180']
        c.append(n)
        j += 1
    a.append(c)

    # Добавляем для печати углы алфа_1
    c = []
    for i in f:
        c.append(grad_r_sek_list(f[i]['α_1']))
    a.append(c)

    c = []
    # добавляем для печати углы бета_1
    for i in f:
        cos_beta_1 = ((f[i]['a_1, m'] ** 2 + f[i]['c_1, m'] ** 2 - f[i]['b_1, m'] ** 2) /
                      (2 * f[i]['a_1, m'] * f[i]['c_1, m']))
        if cos_beta_1 < 0:
            beta_1 = (gradVrad(constants['180']) - np.arcsin(
                    (f[i]['b_1, m'] * np.sin(gradVrad(grad_r_sek_list(f[i]['α_1'])))) /
                    f[i]['a_1, m'])
                    )
        else:
            beta_1 = np.arcsin(f[i]['b_1, m'] * np.sin(gradVrad(grad_r_sek_list(f[i]['α_1']))) /
                               f[i]['a_1, m'])
        c.append(radVgrad(beta_1))
    a.append(c)

    c = []
    # Добавляем для печати углы гамма_1
    for i in f:
        cos_gamma_1 = ((f[i]['a_1, m'] ** 2 + f[i]['b_1, m'] ** 2 - f[i]['c_1, m'] ** 2) /
                       (2 * f[i]['a_1, m'] * f[i]['b_1, m']))
        if cos_gamma_1 < 0:
            gamma_1 = (gradVrad(constants['180']) - np.arcsin(
                (f[i]['c_1, m'] * np.sin(gradVrad(grad_r_sek_list(f[i]['α_1'])))) /
                f[i]['a_1, m'])
                     )
        else:
            gamma_1 = np.arcsin(f[i]['c_1, m'] * np.sin(gradVrad(grad_r_sek_list(f[i]['α_1']))) /
                                f[i]['a_1, m'])
        c.append(radVgrad(gamma_1))
    a.append(c)

    # Считаем невязку в треугольнике в шахте
    c = []
    j = 0
    while j < 3:
        n = a[4][j] + a[5][j] + a[6][j] - constants['180']
        c.append(n)
        j += 1
    a.append(c)

    # Готовим список который вернёт функция для дальнейшей обработки
    angels = [
        a[0],
        a[2],
        a[4],
        a[5],
    ]

    # Формируем табличку для печати
    q = []
    j = 0
    while j < 8:
        c = []
        if j == 3 or j == 7:
            c.append(str(round(a[j][0], 4)) + '"')
            c.append(str(round(a[j][1], 4)) + '"')
            c.append(str(round(a[j][2], 4)) + '"')
        else:
            c.append(sek_r_grad_str(a[j][0]))
            c.append(sek_r_grad_str(a[j][1]))
            c.append(sek_r_grad_str(a[j][2]))
        j += 1
        q.append(c)

    print()
    print('  Таблица 4. Вычисление угловых невязок')
    print('| {:^10} | {:^60} |'.format('', 'Положение отвеса'))
    print('| {:^10} | {:_^60} |'.format('Обозна-', ''))
    print('| {:^10} | {:^18} | {:^18} | {:^18} |'.format('чение',
                                                         'I',
                                                         'II',
                                                         'III'))
    print('| {:-^10} | {:-^18} | {:-^18} | {:-^18} |'.format('', '', '', ''))
    for i in range(len(designations)):
        print('| {:^10} | {:^18} | {:^18} | {:^18} |'.format(designations[i],
                                                             q[i][0],
                                                             q[i][1],
                                                             q[i][2],
                                                             )
              )
    print()
    return angels


def adjustment_l(angels):
    f = position_of_the_plumb
    a = []

    # Вычисляем d_1
    k = 0
    c = []
    for i in f:
        c.append(f[i]['c, m'] * np.cos(gradVrad(angels[0][k])))
        k += 1
    a.append(c)

    # Вычисляем d_2
    k = 0
    c = []
    for i in f:
        c.append(f[i]['a, m'] * np.cos(gradVrad(angels[1][k])))
        k += 1
    a.append(c)

    # Вычисляем с_выч
    o = [a[0][0] + a[1][0],
         a[0][1] + a[1][1],
         a[0][2] + a[1][2],
         ]
    a.append(o)

    # Вычисляем f_s
    k = 0
    c = []
    for i in f:
        c.append(a[2][k] - f[i]['b, m'])
        k += 1
    a.append(c)

    # Вычисляем d_1_1
    k = 0
    c = []
    for i in f:
        c.append(f[i]['b_1, m'] * np.cos(gradVrad(angels[2][k])))
        k += 1
    a.append(c)

    # Вычисляем d_1_2
    k = 0
    c = []
    for i in f:
        c.append(f[i]['a_1, m'] * np.cos(gradVrad(angels[3][k])))
        k += 1
    a.append(c)

    # Вычисляем с_выч_1
    m = [a[4][0] + a[5][0],
         a[4][1] + a[5][1],
         a[4][2] + a[5][2],
         ]
    a.append(m)

    # Вычисляем f_s_1
    k = 0
    c = []
    for i in f:
        c.append(a[6][k] - f[i]['c_1, m'])
        k += 1
    a.append(c)

    f_s = [
        a[3],
        a[7],
    ]

    print('  Таблица 5. Вычисление линейных невязок')
    print('| {:^10} | {:^60} |'.format('', 'Положение отвеса'))
    print('| {:^10} | {:_^60} |'.format('Обозна-', ''))
    print('| {:^10} | {:^18} | {:^18} | {:^18} |'.format('чение',
                                                         'I',
                                                         'II',
                                                         'III'))
    print('| {:-^10} | {:-^18} | {:-^18} | {:-^18} |'.format('', '', '', ''))
    for i in range(len(designations_t5)):
        print('| {:^10} | {:^18} | {:^18} | {:^18} |'.format(designations_t5[i],
                                                             round(a[i][0], 4),
                                                             round(a[i][1], 4),
                                                             round(a[i][2], 4),
                                                             )
              )
    return f_s


def calculating_corrections(f_s):
    c = []
    a = []
    for i in range(len(f_s[0])):
        c.append(- (f_s[0][i] / 3))
    a.append(c)
    a.append(c)
    c = []
    for i in range(len(f_s[0])):
        c.append((f_s[0][i] / 3))
    a.append(c)

    c =[]
    for i in range(len(f_s[1])):
        c.append(- (f_s[1][i] / 3))
    a.append(c)
    a.append(c)
    c = []
    for i in range(len(f_s[1])):
        c.append((f_s[1][i] / 3))
    a.append(c)

    print()
    print('  Таблица 6. Вычисление поправок в стороны')
    print('| {:^10} | {:^60} |'.format('', 'Положение отвеса'))
    print('| {:^10} | {:_^60} |'.format('Обозна-', ''))
    print('| {:^10} | {:^18} | {:^18} | {:^18} |'.format('чение',
                                                         'I',
                                                         'II',
                                                         'III'))
    print('| {:-^10} | {:-^18} | {:-^18} | {:-^18} |'.format('', '', '', ''))
    for i in range(len(designations_t6)):
        print('| {:^10} | {:^18} | {:^18} | {:^18} |'.format(designations_t6[i],
                                                             round(a[i][0], 4),
                                                             round(a[i][1], 4),
                                                             round(a[i][2], 4),
                                                             )
              )

    # Вносим поправки в измеренные стороны
    f = position_of_the_plumb
    b = []
    c = []
    k = 0
    for i in f:
        c.append(f[i]['a, m'] + a[0][k])
        k += 1
    b.append(c)

    c = []
    k = 0
    for i in f:
        c.append(f[i]['b, m'] + a[1][k])
        k += 1
    b.append(c)

    c = []
    k = 0
    for i in f:
        c.append(f[i]['c, m'] + a[2][k])
        k += 1
    b.append(c)

    c = []
    k = 0
    for i in f:
        c.append(f[i]['a_1, m'] + a[3][k])
        k += 1
    b.append(c)
    c = []
    k = 0

    for i in f:
        c.append(f[i]['b_1, m'] + a[4][k])
        k += 1
    b.append(c)

    c = []
    k = 0
    for i in f:
        c.append(f[i]['c_1, m'] + a[5][k])
        k += 1
    b.append(c)

    print()
    print('  Таблица 7. Исправленные стороны')
    print('| {:^10} | {:^60} |'.format('', 'Положение отвеса'))
    print('| {:^10} | {:_^60} |'.format('Обозна-', ''))
    print('| {:^10} | {:^18} | {:^18} | {:^18} |'.format('чение',
                                                         'I',
                                                         'II',
                                                         'III'))
    print('| {:-^10} | {:-^18} | {:-^18} | {:-^18} |'.format('', '', '', ''))
    for i in range(len(designations_t7)):
        print('| {:^10} | {:^18} | {:^18} | {:^18} |'.format(designations_t7[i],
                                                             round(b[i][0], 4),
                                                             round(b[i][1], 4),
                                                             round(b[i][2], 4),
                                                             )
              )

    return b


def calculating_log(parties, angels):
    angels_alpha = angels[::2]

    a = []
    # Вычисляем c/a b/a для надземного треугольника
    k = 2
    while k >= 1:
        c = []
        for i in range(3):
            c.append(parties[k][i] / parties[0][i])
        a.append(c)
        k -= 1

    # Вычисляем c/a b/a для подземного треугольника
    k = 5
    while k >= 4:
        c = []
        for i in range(3):
            c.append(parties[k][i] / parties[3][i])
        a.append(c)
        k -= 1

    # Вычисляем синусы alpha
    sin_alpha = []
    for i in angels_alpha:
        c = []
        for k in i:
            c.append(np.sin(gradVrad(k)))
        sin_alpha.append(c)

    # Вычисляем косинусы углов beta и gamma ЗДЕСЬ БАГ (считает нижний косинус по верхней стороне a)
    cos_angels = []
    j = 1
    while j < 6:
        c = []
        for i in range(3):
            c.append(
                (parties[0][i] ** 2 + parties[j][i] ** 2 - parties[j + 1][i] ** 2) /
                (2 * parties[0][i] * parties[j][i])
            )
        cos_angels.append(c)
        c = []
        for i in range(3):
            c.append(
                (parties[0][i] ** 2 + parties[j+1][i] ** 2 - parties[j][i] ** 2) /
                (2 * parties[0][i] * parties[j+1][i])
            )
        cos_angels.append(c)
        j += 3

    # Вычисляем синусы углов beta и gamma
    sin_angels = []
    j = 0
    while j < 4:
        c = []
        if j < 2:
            for i in range(3):
                c.append(a[j][i] * sin_alpha[0][i])
        else:
            for i in range(3):
                c.append(a[j][i] * sin_alpha[0][i])
        sin_angels.append(c)
        j += 1

    sin_angels = np.array(sin_angels)
    cos_angels = np.array(cos_angels)
    angels_alpha = np.array(angels_alpha)
    t = np.arcsin(sin_angels)
    for i in range(len(t)):
        for j in range(len(t[i])):
            t[i][j] = radVgrad(t[i][j])

    for i in range(len(cos_angels)):
        for j in range(len(cos_angels[i])):
            if cos_angels[i][j] < 0:
                t[i][j] = (constants['180'] - t[i][j])

    result = []
    result.append(angels_alpha[0])
    result.append(t[0])
    result.append(t[1])
    c = []
    for i in range(3):
        c.append(result[0][i] + result[1][i] + result[2][i])
    result.append(c)

    result.append(angels_alpha[1])
    result.append(t[2])
    result.append(t[3])
    c = []
    for i in range(3):
        c.append(result[4][i] + result[5][i] + result[6][i])
    result.append(c)

    print()
    print('  Таблица 8. Результаты уравнивания')
    print('| {:^10} | {:^60} |'.format('', 'Положение отвеса'))
    print('| {:^10} | {:_^60} |'.format('Обозна-', ''))
    print('| {:^10} | {:^18} | {:^18} | {:^18} |'.format('чение',
                                                         'I',
                                                         'II',
                                                         'III'))
    print('| {:-^10} | {:-^18} | {:-^18} | {:-^18} |'.format('', '', '', ''))
    for i in range(len(designations_t8)):
        print('| {:^10} | {:^18} | {:^18} | {:^18} |'.format(designations_t8[i],
                                                             sek_r_grad_str(result[i][0]),
                                                             sek_r_grad_str(result[i][1]),
                                                             sek_r_grad_str(result[i][2]),
                                                             )
              )
    print()

    # print(t)
    # print(angels_alpha)
    # print(sek_r_grad_str(t[1][2] + t[0][2] + angels_alpha[0][2]))


    # print()
    # for i in a:
    #     print('{} {} {}'.format(round(i[0], 4), round(i[1], 4), round(i[2], 4)))
    # print()
    # for i in cos_angels:
    #     print('{} {} {}'.format(round(i[0], 4), round(i[1], 4), round(i[2], 4)))
    # print()
    # for i in sin_angels:
    #     print('{} {} {}'.format(round(i[0], 4), round(i[1], 4), round(i[2], 4)))



def the_connecting_triangle():
    measure_log()
    omega_check()
    alpha_check()
    angels = adjustment_a()
    f_s = adjustment_l(angels)
    parties = calculating_corrections(f_s)
    calculating_log(parties, angels)

if __name__ == "__main__":
    result = the_connecting_triangle()
