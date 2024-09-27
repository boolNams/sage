M, P, k, q1, q2, E, J = var('M P k q1 q2 E J')		# переменные заданные по условию
omega0, teta0, R1, R2, R3, R4 = var('omega0 teta0 R1 R2 R3 R4')			# переменные подлежащие определению

x = var('x')										# координата
omega = function('omega')(x)						# прогиб балки
teta  = omega.diff(x)								# углы поворота сечений балки
MM    = E * J * teta.diff(x)						# момент балки

print("complete")