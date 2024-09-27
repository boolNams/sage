#
# ПРОГРАММА РАСЧЕТА БАЛКИ С ПОМОЩЬЮ СИСТЕМЫ КОМПЬЮТЕРНОЙ АЛГЕБРЫ | SAGE
# ДЛИНА БАЛКИ | 25
#

M, P, k, qA, qB, E, J = var('M P k qA qB E J')							# переменные заданные по условию
M = 2
P = 2
k = 1/3
qA = 2/3
qB = 1
E = 5
J = 1
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
#==================================================

eq1 = RA+RB+RC+RD + P + (1/2)*(Bq-Aq)*(qA+qB) + k*wk == 0
eq2 = RA*AR+RB*BR+RC*CR+RD*DR + P*AP + M + k*wk*Ak + 

print("complete")