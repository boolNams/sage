#
# ПРОГРАММА РАСЧЕТА БАЛКИ С ПОМОЩЬЮ СИСТЕМЫ КОМПЬЮТЕРНОЙ АЛГЕБРЫ | SAGE
# ДЛИНА БАЛКИ | 25
#
#================================================== ПРОГИБ
M, P, k, qA, qB, E, J = var('M P k qA qB E J')							# переменные заданные по условию
M = 1
P = -2
k = 3
qA = -1
qB = 3
E = 100
J = 100
w0, wk, teta0, RA, RB, RC, RD = var('w0 wk teta0 R1 R2 R3 R4')			# переменные подлежащие определению

x = var('x')															# координата (0,25)     
prd = E*J
div = 1/(E*J)  

a = var('a')															# параметр многочленов 
y2(x, a) = (1/2) * (x - a)^2
y3(x, a) = (1/6) * (x - a)^3
y4(x, a) = (1/24) * (x - a)^4
y5(x, a) = (1/120) * (x - a)^5                                                   

AR, BR, CR, DR = 1, 8, 16, 25											# координаты упругих опор
hAR(x,RA) = RA*heaviside(x - AR)										# функции хевисайда для реакций
hBR(x,RB) = RB*heaviside(x - BR)
hCR(x,RC) = RC*heaviside(x - CR)
hDR(x,RD) = RD*heaviside(x - DR)
wR = y3(x,AR)*hAR + y3(x,BR)*hBR + y3(x,CR)*hCR + y3(x,DR)*hDR			# прогиб (опоры)

AM = 13																	# координата приложенного момента M
hAM(x) = M*heaviside(x - AM)											# функция хевисайда для момента M
wM = y2(x,AM)*hAM														# прогиб (момент)

AP = 20																	# координата приложенной силы P
hAP(x) = P*heaviside(x - AP)											# функция хевисайда для силы P
wP = y3(x,AP)*hAP             											# прогиб (сила)

Aq, Bq = 10, 22 														# координаты распределенной нагрузки
tg = (qB - qA) / (Bq - Aq)												# тангенс угла наклона
hAq(x) = heaviside(x - Aq)												# функция хевисайда основной нагрузки (+)
hBq(x) = heaviside(x - Bq)												# функция хевисайда компенсирующей нагрузки (-)
wAq = (qA*y4(x,Aq) + k*y5(x,Aq))*hAq									# прогиб (основная нагрузка +) 
wBq = -(qB*y4(x,Bq) + k*y5(x,Bq))*hBq									# прогиб (компенсирующая нагрузка -)

Ak = 23																	# координата пружины жетскости k
hAk(x, wk) = k*wk*heaviside(x - Ak)										# функция хевисайда пружины k
wAk = y3(x,Ak)*hAk														# прогиб (пружина)

wConst(x, w0, teta0) = w0 + teta0*x										# прогиб (константы интегрирвания)

#==================================================
w = wConst + div * (wR + wM + wP + wAq + wBq + wAk)						# прогиб
#================================================== СИСТЕМА

eq1 = RA+RB+RC+RD + P + (1/2)*(Bq-Aq)*(qA+qB) + k*wk == 0   			# уравнение равновесия (1)

MR = RA*AR+RB*BR+RC*CR+RD*DR                                      		# момент создаваемый опорами
MQ = (1/2)*(qA*(Bq^2-Aq^2) + (qB-qA)*(Bq - Aq))							# момент создаваемый распределенной нагрузкой
eq2 = MR + P*AP + M + k*wk*Ak + MQ == 0									# уравнение моментов (2)

eq3 = w(x = AR) == 0              										# уравнения на опорах (3-6)
eq4 = w(x = BR) == 0
eq5 = w(x = CR) == 0
eq6 = w(x = DR) == 0

eq7 = w(x = Ak) - wk == 0      											# уравнение на пружину (7)

eq = [eq1,eq2,eq3,eq4,eq5,eq6,eq7] 										# система уравнений
sol = solve(eq, w0, wk, teta0, RA, RB, RC, RD, solution_dict = True)[0]
w0, wk, teta0 = sol[w0], sol[wk], sol[teta0]
RA, RB, RC, RD = sol[RA], sol[RB], sol[RC], sol[RD]

w = w(RA=RA, RB=RB, RC=RC, RD=RD, teta0=teta0, w0=w0, wk=wk)			# функция прогиба от x

#================================================== ГРАФИКА

import numpy as np
import matplotlib.pyplot as plt

N = 100                                                                 # количество точек на графике

arg = np.linspace(0,25,N)												# построение прогиба
val = np.zeros(N)
for i in range(0,N):
	val[i] = w(x = arg[i])
plt.plot(arg, val)

arg = np.array([AR, BR, CR, DR], dtype = float)							# построение опор
val = np.zeros(4)
plt.plot(arg, val, 'ro', label = "R")

arg = np.array([AP], dtype = float)										# построение силы P
val = np.zeros(1)
plt.plot(arg,val, 'bo', label = "P")

arg = np.array([AM], dtype = float)										# построение момента M
val = np.zeros(1)
plt.plot(arg,val, 'mo', label = "M")

arg = np.array([Ak], dtype = float)										# построение пружины k
val = np.zeros(1)
plt.plot(arg,val, 'yo', label = "k")

arg = np.array([Aq, Bq], dtype = float)									# построение распределенной нагрузки
val = np.array([qA, qB], dtype = float)
plt.plot(arg,val, color = 'g', alpha = 0.3, label = "q")
plt.fill_between(arg, val,color = 'g', alpha = 0.3)



plt.grid()
plt.legend()
plt.xlim(0,25)
plt.show()

#==================================================
print("complete")