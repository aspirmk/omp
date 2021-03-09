# -*- coding: utf-8 -*-
"""
Модуль определения места повреждения (КЗ) на высоковольтных ВЛ
Базируется на упрощенных способах одностороннего замера, описанных в
книге Аржанникова Е.А. Лукоянова В.Ю. Мисриханова М.Ш. Определение места
короткого замыкания на высоковольтных линиях электропередачи /
Под ред. В.А. Шуина. - М: Энергоатомиздат, 2003.
"""
import numpy as np

r2d = 180/np.pi
a =  -0.5 + 0.5j*np.sqrt(3)
a2 = -0.5 - 0.5j*np.sqrt(3)
a0 = 1.0 + 0.0j
v1 = np.array([a0,a,a2])/3
v2 = np.array([a0,a2,a])/3
v0 = np.array([a0,a0,a0])/3
vA = np.array([a0,a0,a0])
vB = np.array([a2,a ,a0])
vC = np.array([a, a2,a0])
vAB = vA - vB
vBC = vB - vC
vCA = vC - vA
Mf2s = np.array([v1,v2,v0])
Ms2f = np.array([vA,vB,vC])
Ms2ff = np.array([vAB,vBC,vCA])


class uch:
    """Класс реализующий участок линии с однородными параметрами """
    def __init__(self, line, name, L, Z1, Z0, Z0otv=0, M0=0, Z0otvpl=0, I1n=0):
        """Конструктор участка линии с однородными параметрами
        line - ссылка на объект линии к которой будет добавлен данный участок
        name - наименование участка ВЛ
        L - длина участка линии, км
        Z1 - удельное сопротивление прямой последовательности участка ВЛ
        Z0 - удельное сопротивление нулевой последовательности участка ВЛ
        Z0otv - сопротивление нулевой последовательности на землю отпайки с
           трансформатором с глухозаземленной нейтралью
        M0 - удельное сопротивление взаимоиндукции нулевой последовательности
           с параллельной ВЛ, предполагается, что данная параллельная ВЛ
           подключена к той же секции шин и обладает соответствующим
           удельным сопротивлением нулевой последовательности
         Z0otvpl - сопротивление нулевой последовательности на землю отпайки с
           трансформатором с глухозаземленной нейтралью от параллельной ВЛ
         I1n - вектор  тока нагрузки подстанции фазы А
        """
        line.bu.append(self)
        self.name = name
        self.L = L
        self.Z1 = Z1
        self.Z0 = Z0
        self.Z0otv = Z0otv
        self.M0 = M0
        self.Z0otvpl = Z0otvpl
        self.I1n = I1n
        self.K = (Z0 - Z1) / Z1
        self.Km = M0 / Z1


class lin:
    """Класс реализующий линию и предназначен для хранения ее мат.модели и для
    расчета ОМП"""
    def __init__(self, name):
        """Конструктор ВЛ
        name - Наименование ВЛ
        """
        self.name = name
        self.bu = []

    def calc(self, vidKZ, Uabc, Iabc, I1prav=0, I03pl=0):
        """Метод осуществляющий расчет ОМП
        vidKZ - Вид КЗ:
            'ABC' - Для расчета ОМП при трехфазных КЗ
            'AB', 'BC', 'CA' - Для расчета ОМП при междуфазных КЗ, в качестве
                поляризующего тока используется составляющая тока обратной
                последовательности неповрежденной фазы
            'A0', 'B0', 'C0' - Для расчета ОМП при КЗ на землю, в качестве
                поляризующего тока используется составляющая тока нулевой
                последовательности
            'A0I2', 'B0I2', 'C0I2' - Для расчета ОМП при КЗ на землю, в качестве
                поляризующего тока используется составляющая тока обратной
                последовательности
            'A0I1av', 'B0I1av', 'C0I1av' - Для расчета ОМП при КЗ на землю, в качестве
                поляризующего тока используется аварийная составляющая тока
                прямой последовательности
        Uabc - вектор фазных напряжений в начале ВЛ, месте установки РАС
        Iabc - вектор фазных токов в начале ВЛ, месте установки РАС
        I1prav - вектор предаварийного тока фазы А прямой последовательност,
           используется для расчета ОМП при КЗ на землю с полярицией аварийной
           составляющей тока прямой последовательности
        I03pl - вектор утроенного тока нулевой последовательности параллельной ВЛ
        """
        U120 = Mf2s @ Uabc
        I120 = Mf2s @ Iabc
        U0pl = U120[2]
        I0pl = I03pl/3
        print('РАСЧЕТ МЕСТА ПОВРЕЖДЕНИЯ ВЛ')
        print('Наименование ВЛ:', self.name)
        print('Параметры предаварийного и аварийного режимов:')
        print(Str(U120, I120, I1prav, I03pl))
        if vidKZ in ('A0I1av','B0I1av','C0I1av'):
            I120[0] -= I1prav
        Lres = 0
        for ku in self.bu:
            Lu = fun_omp[vidKZ](U120,I120,I0pl,ku)
            if Lu <= 0:
                break
            elif Lu <= ku.L:
                Lres += Lu
                break
            else:# Lu > ku.L
                if not ku is self.bu[-1]:
                    Lres += ku.L
                    U120 -= ku.L * np.array([ku.Z1,ku.Z1,ku.Z0]) * I120
                    U120[2] -= ku.L * ku.M0 * I0pl
                    U0pl -= ku.L * (ku.M0 * I120[2] + ku.Z0 * I0pl)
                    if not vidKZ in ('A0I1av','B0I1av','C0I1av'):
                        I120[0] -= ku.I1n
                    if np.abs(ku.Z0otv) > 0:
                        I120[2] -= U120[2] / ku.Z0otv
                    if np.abs(ku.Z0otvpl) > 0:
                        I0pl -= U0pl / ku.Z0otvpl
                else:
                    Lres += Lu
        print('Результат расчета ОМП для вида КЗ вида {0} - {1:6.2f} км'.format(vidKZ, Lres))

def fomp(u,ud,iop,val):
    if val=='imag':
        return np.imag(u/iop)/np.imag(ud/iop)
    elif val=='real':
        return np.real(u/iop)/np.real(ud/iop)
    else:
        return 0

fun_omp=dict({'ABC' : lambda U120,I120,I0pl,ku : fomp(vAB @ U120, ku.Z1 * (vAB @ I120), vC @ I120, 'real'),
               'AB' : lambda U120,I120,I0pl,ku : fomp(vAB @ U120, ku.Z1 * (vAB @ I120), a2 * I120[1], 'real'),
               'BC' : lambda U120,I120,I0pl,ku : fomp(vBC @ U120, ku.Z1 * (vBC @ I120), a0 * I120[1], 'real'),
               'CA' : lambda U120,I120,I0pl,ku : fomp(vCA @ U120, ku.Z1 * (vCA @ I120), a * I120[1], 'real'),
               'A0' : lambda U120,I120,I0pl,ku : fomp(vA @ U120, ku.Z1 * (vA @ I120 + ku.K * I120[2] + ku.Km * I0pl), I120[2],'imag'),
               'B0' : lambda U120,I120,I0pl,ku : fomp(vB @ U120, ku.Z1 * (vB @ I120 + ku.K * I120[2] + ku.Km * I0pl), I120[2],'imag'),
               'C0' : lambda U120,I120,I0pl,ku : fomp(vC @ U120, ku.Z1 * (vC @ I120 + ku.K * I120[2] + ku.Km * I0pl), I120[2],'imag'),
               'A0I2' : lambda U120,I120,I0pl,ku : fomp(vA @ U120, ku.Z1 * (vA @ I120 + ku.K * I120[2] + ku.Km * I0pl), I120[1],'imag'),
               'B0I2' : lambda U120,I120,I0pl,ku : fomp(vB @ U120, ku.Z1 * (vB @ I120 + ku.K * I120[2] + ku.Km * I0pl), a * I120[1],'imag'),
               'C0I2' : lambda U120,I120,I0pl,ku : fomp(vC @ U120, ku.Z1 * (vC @ I120 + ku.K * I120[2] + ku.Km * I0pl), a2 * I120[1],'imag'),
               'A0I1av' : lambda U120,I120,I0pl,ku : fomp(vA @ U120, ku.Z1 * (vA @ I120 + ku.K * I120[2] + ku.Km * I0pl), I120[0],'imag'),
               'B0I1av' : lambda U120,I120,I0pl,ku : fomp(vB @ U120, ku.Z1 * (vB @ I120 + ku.K * I120[2] + ku.Km * I0pl), a2 * I120[0],'imag'),
               'C0I1av' : lambda U120,I120,I0pl,ku : fomp(vC @ U120, ku.Z1 * (vC @ I120 + ku.K * I120[2] + ku.Km * I0pl), a * I120[0],'imag')
              })

def Str(u120,i120,I1prav,I03pl):
    strI1prav = "I1пр.ав = {0:>7.0f} ∠ {1:>6.1f}\n"
    strUABC = "| UA  = {0:>7.0f} ∠ {1:>6.1f} | UB  = {2:>7.0f} ∠ {3:>6.1f} | UC  = {4:>7.0f} ∠ {5:>6.1f} |\n"
    strU120 = "| U1  = {0:>7.0f} ∠ {1:>6.1f} | U2  = {2:>7.0f} ∠ {3:>6.1f} | 3U0 = {4:>7.0f} ∠ {5:>6.1f} |\n"
    strUAB_BC_CA = "| UAB = {0:>7.0f} ∠ {1:>6.1f} | UBC = {2:>7.0f} ∠ {3:>6.1f} | UCA = {4:>7.0f} ∠ {5:>6.1f} |\n"
    strIABC = "| IA  = {0:>7.0f} ∠ {1:>6.1f} | IB  = {2:>7.0f} ∠ {3:>6.1f} | IC  = {4:>7.0f} ∠ {5:>6.1f} |\n"
    strI120 = "| I1  = {0:>7.0f} ∠ {1:>6.1f} | I2  = {2:>7.0f} ∠ {3:>6.1f} | 3I0 = {4:>7.0f} ∠ {5:>6.1f} |\n"
    strIAB_BC_CA = "| IAB = {0:>7.0f} ∠ {1:>6.1f} | IBC = {2:>7.0f} ∠ {3:>6.1f} | ICA = {4:>7.0f} ∠ {5:>6.1f} |\n"
    str3I0parvl = "3I0пар.вл = {0:>7.0f} ∠ {1:>6.1f}"
    u1,u2,u0 = u120
    uA,uB,uC = Ms2f @ u120
    uAB,uBC,uCA = Ms2ff @ u120
    i1,i2,i0 = i120
    iA,iB,iC = Ms2f @ i120
    iAB,iBC,iCA = Ms2ff @ i120
    resstr = []
    resstr.append(strI1prav.format(np.abs(I1prav),r2d*np.angle(I1prav)))
    resstr.append(strUABC.format(np.abs(uA),r2d*np.angle(uA),np.abs(uB),r2d*np.angle(uB),np.abs(uC),r2d*np.angle(uC)))
    resstr.append(strU120.format(np.abs(u1),r2d*np.angle(u1),np.abs(u2),r2d*np.angle(u2),np.abs(3*u0),r2d*np.angle(u0)))
    resstr.append(strUAB_BC_CA.format(np.abs(uAB),r2d*np.angle(uAB),np.abs(uBC),r2d*np.angle(uBC),np.abs(uCA),r2d*np.angle(uCA)))
    resstr.append(strIABC.format(np.abs(iA),r2d*np.angle(iA),np.abs(iB),r2d*np.angle(iB),np.abs(iC),r2d*np.angle(iC)))
    resstr.append(strI120.format(np.abs(i1),r2d*np.angle(i1),np.abs(i2),r2d*np.angle(i2),np.abs(3*i0),r2d*np.angle(i0)))
    resstr.append(strIAB_BC_CA.format(np.abs(iAB),r2d*np.angle(iAB),np.abs(iBC),r2d*np.angle(iBC),np.abs(iCA),r2d*np.angle(iCA)))
    resstr.append(str3I0parvl.format(np.abs(I03pl),r2d*np.angle(I03pl)))
    return ''.join(resstr)