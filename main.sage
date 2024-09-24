M, P, k, q1, q2, E, J = var('M P k q1 q2 E J')		# переменные заданные по условию
omega0, teta0, R = var('omega0 teta0 R')			# переменные подлежащие определению

x = var('x')										# координата
omega = function('omega')(x)						# прогиб балки
teta  = omega.diff(x)								# углы поворота сечений балки
MM    = E * J * teta.diff(x)						# момент балки

print("complete")