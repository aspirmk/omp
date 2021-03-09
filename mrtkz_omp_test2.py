#Импортирование модуля расчета ТКЗ (mrtkz3.py) и
#модуля расчета параметров воздушных линий  PVL,
#которые должны находиться в той же папке, где и настоящий файл
# import PVL5 as PVL
import mrtkz3 as mrtkz
import momp as omp

L1 = omp.lin('Тестовая ВЛ')
U1 = omp.uch(L1, 'Участок 1', 10, 0.2+0.4j, 0.4+1.2j, M0=0.15+0.9j, Z0otvpl=60j, I1n=0)
U2 = omp.uch(L1, 'Участок 2', 15, 0.3+0.4j, 0.5+1.3j, Z0otv=60j, M0=0.15+0.85j, Z0otvpl=60j)
U3 = omp.uch(L1, 'Участок 3', 15, 0.3+0.4j, 0.5+1.3j, Z0otv=60j, I1n=0)
# U4 = omp.uch('Участок 4', L1, 15, 0.3+0.4j, 0.5+1.3j)

#Создание расчетной модели
mdl=mrtkz.Model()

#Создание узлов
q1 = mrtkz.Q(mdl,'Sys1')
q2 = mrtkz.Q(mdl,'u1')
q3 = mrtkz.Q(mdl,'u2')
q4 = mrtkz.Q(mdl,'u3')

q2p = mrtkz.Q(mdl,'u1')
q3p = mrtkz.Q(mdl,'u2')
q4p = mrtkz.Q(mdl,'Sys2')

#Создание ветвей энергосистем
Sys1 = mrtkz.P(mdl,'Sys1',0,q1,(2j,2j,3j),E=(65000,0,0))
LU1 = mrtkz.P(mdl,'U1',q1,q2,(2+4j,2+4j,4+12j))

LU2 = mrtkz.P(mdl,'U2',q2,q3,(4.5+6j,4.5+6j,7.5+19.5j))
TU2 = mrtkz.P(mdl,'T2',q3,0,(9999,9999,60j))

LU3 = mrtkz.P(mdl,'U3',q3,q4,(4.5+6j,4.5+6j,7.5+19.5j))
TU3 = mrtkz.P(mdl,'T3',q4,0,(9999,9999,60j))


LU1p = mrtkz.P(mdl,'U1',q1,q2p,(2+4j,2+4j,4+12j))
TU1p = mrtkz.P(mdl,'T1p',q2p,0,(9999,9999,60j))
LU2p = mrtkz.P(mdl,'U2',q2p,q3p,(4.5+6j,4.5+6j,7.5+19.5j))
TU2p = mrtkz.P(mdl,'T2p',q3p,0,(9999,9999,60j))
LU3p = mrtkz.P(mdl,'U3',q3p,q4p,(4.5+6j,4.5+6j,7.5+19.5j))
Sys2 = mrtkz.P(mdl,'Sys2',0,q4p,(2j,2j,3j),E=(65000,0,0))

mrtkz.M(mdl,'LU1-LU1p',LU1,LU1p,10*(0.15+0.9j),10*(0.15+0.9j))
mrtkz.M(mdl,'LU2-LU2p',LU2,LU2p,15*(0.15+0.85j),15*(0.15+0.85j))

print('Однофазные КЗ')
#Создание однофазного КЗ
KZ1 = mrtkz.N(mdl,'Тестовое КЗ',q2,'A0')
#  Проверка на вырожденность
mdl.Test4Singularity()
#Формирование разреженной СЛАУ и расчет электрических параметров
mdl.Calc()
#Вывод результатов расчета для короткого замыкания
# KZ1.res()
# LU1.res1()
L1.calc('A0I1av', LU1.q1UABC, LU1.q1IABC,I03pl=3*LU1p.q1I0)

print()
#Очистка модели от КЗ и обрывов за исключением типа 'N0'
mdl.ClearN()
KZ1 = mrtkz.N(mdl,'Тестовое КЗ',q3,'B0')
mdl.Calc()
L1.calc('B0I1av', LU1.q1UABC, LU1.q1IABC,I03pl=3*LU1p.q1I0)

print()
#Очистка модели от КЗ и обрывов за исключением типа 'N0'
mdl.ClearN()
KZ1 = mrtkz.N(mdl,'Тестовое КЗ',q4,'C0')
mdl.Calc()
L1.calc('C0I1av', LU1.q1UABC, LU1.q1IABC,I03pl=3*LU1p.q1I0)

print()
print('Двухфазные КЗ на землю')
mdl.ClearN()
KZ1 = mrtkz.N(mdl,'Тестовое КЗ',q2,'AB0')
mdl.Calc()
L1.calc('A0I1av', LU1.q1UABC, LU1.q1IABC,I03pl=3*LU1p.q1I0)
L1.calc('B0I1av', LU1.q1UABC, LU1.q1IABC,I03pl=3*LU1p.q1I0)

print()
mdl.ClearN()
KZ1 = mrtkz.N(mdl,'Тестовое КЗ',q3,'BC0')
mdl.Calc()
# LU1.res1()
L1.calc('B0I1av', LU1.q1UABC, LU1.q1IABC,I03pl=3*LU1p.q1I0)
L1.calc('C0I1av', LU1.q1UABC, LU1.q1IABC,I03pl=3*LU1p.q1I0)

print()
mdl.ClearN()
KZ1 = mrtkz.N(mdl,'Тестовое КЗ',q4,'CA0')
mdl.Calc()
L1.calc('C0I1av', LU1.q1UABC, LU1.q1IABC,I03pl=3*LU1p.q1I0)
L1.calc('A0I1av', LU1.q1UABC, LU1.q1IABC,I03pl=3*LU1p.q1I0)

print()
print('Двухфазные КЗ')
mdl.ClearN()
KZ1 = mrtkz.N(mdl,'Тестовое КЗ',q2,'AB')
mdl.Calc()
L1.calc('AB', LU1.q1UABC, LU1.q1IABC,I03pl=3*LU1p.q1I0)

print()
mdl.ClearN()
KZ1 = mrtkz.N(mdl,'Тестовое КЗ',q3,'BC')
mdl.Calc()
L1.calc('BC', LU1.q1UABC, LU1.q1IABC,I03pl=3*LU1p.q1I0)

print()
mdl.ClearN()
KZ1 = mrtkz.N(mdl,'Тестовое КЗ',q4,'CA')
mdl.Calc()
L1.calc('CA', LU1.q1UABC, LU1.q1IABC,I03pl=3*LU1p.q1I0)

print('Однофазное КЗ через R')
#Создание однофазного КЗ
KZ1 = mrtkz.N(mdl,'Тестовое КЗ',q3,'A0r',r=1)
#  Проверка на вырожденность
mdl.Test4Singularity()
#Формирование разреженной СЛАУ и расчет электрических параметров
mdl.Calc()
#Вывод результатов расчета для короткого замыкания
# KZ1.res()
# LU1.res1()
L1.calc('A0I1av', LU1.q1UABC, LU1.q1IABC,I03pl=3*LU1p.q1I0)


print('Двухфазное КЗ через R')
print()
mdl.ClearN()
KZ1 = mrtkz.N(mdl,'Тестовое КЗ',q3,'BCr',r=5)
mdl.Calc()
L1.calc('BC', LU1.q1UABC, LU1.q1IABC,I03pl=3*LU1p.q1I0)
