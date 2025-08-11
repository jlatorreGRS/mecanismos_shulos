import numpy as np

# Trabajo Teoría de Mecanismos
# Problemas de velocidad y aceleración

# Datos del mecanismo
L2 = 25  # mm
L4 = 130
L5 = 175
L6 = 110
x1 = 17.67
y1 = 17.67
x2 = -16.15
y2 = 74.91
x3 = 156.433
y3 = 45.95
xA = 0
yA = 0
xB = 50
yB = -37
xC = 186
yC = -60
w2 = 4  # rad/s
alpha = np.deg2rad(45)
beta = np.deg2rad(49.8863)
sigma = np.deg2rad(120.5891)  # sigma = 59.4109
ro = np.deg2rad(115.117)
tita = np.deg2rad(180 - 74.4079)
delta = np.deg2rad(9.5245)

# Vector de coordenadas generalizadas
q = np.array([x1, y1, x2, y2, alpha, sigma, x3, y3, tita])

# Evaluación de ecuaciones de restricción "phi"
phi = np.zeros(8)
phi[0] = x1**2 + y1**2 - L2**2
phi[1] = (x2 - xB) ** 2 + (y2 - yB) ** 2 - L4 ** 2
phi[2] = (x2 - xB) * (y1 - yB) - (x1 - xB) * (y2 - yB)
if abs(np.sin(alpha)) > 0.7:
    phi[3] = x1 - L2 * np.cos(alpha)
else:
    phi[3] = y1 - L2 * np.sin(alpha)
phi[4] = x2 - xB - L4 * np.cos(sigma)
phi[5] = (x3 - x2) ** 2 + (y3 - y2) ** 2 - L5 ** 2
phi[6] = (x3 - xC) ** 2 + (y3 - yC) ** 2 - L6 ** 2
phi[7] = x3 - xC - L6 * np.cos(tita)

# Evaluación del Jacobiano
phiq = np.zeros((9, 9))
phiq[0, :] = [2 * x1, 2 * y1, 0, 0, 0, 0, 0, 0, 0]
phiq[1, :] = [0, 0, 2 * (x2 - xB), 2 * (y2 - yB), 0, 0, 0, 0, 0]
phiq[2, :] = [-(y2 - yB), (x2 - xB), (y1 - yB), -(x1 - xB), 0, 0, 0, 0, 0]
if abs(np.sin(alpha)) > 0.7:
    phiq[3, :] = [1, 0, 0, 0, L2 * np.sin(alpha), 0, 0, 0, 0]
else:
    phiq[3, :] = [0, 1, 0, 0, -L2 * np.cos(alpha), 0, 0, 0, 0]
phiq[4, :] = [0, 0, 1, 0, 0, L4 * np.sin(sigma), 0, 0, 0]
phiq[5, :] = [0, 0, -2 * (x3 - x2), -2 * (y3 - y2), 0, 0, 2 * (x3 - x2), 2 * (y3 - y2), 0]
phiq[6, :] = [0, 0, 0, 0, 0, 0, 2 * (x3 - xC), 2 * (y3 - yC), 0]
phiq[7, :] = [0, 0, 0, 0, 0, 0, 1, 0, L6 * np.sin(tita)]
phiq[8, :] = [0, 0, 0, 0, 1, 0, 0, 0, 0]

# Vector derecho de la ecuación de velocidades
rhsv = np.array([0, 0, 0, 0, 0, 0, 0, 0, w2])

# Resolución de las velocidades
qp = np.linalg.solve(phiq, rhsv)

# Extracción de valores del vector qp
x1p, y1p, x2p, y2p, w2, w4, x3p, y3p, w6 = qp

# Evaluación de la derivada temporal del Jacobiano
phiqp = np.zeros((8, 9))
phiqp[0, :] = [2 * x1p, 2 * y1p, 0, 0, 0, 0, 0, 0, 0]
phiqp[1, :] = [0, 0, 2 * x2p, 2 * y2p, 0, 0, 0, 0, 0]
phiqp[2, :] = [-y2p, x2p, y1p, -x1p, 0, 0, 0, 0, 0]
if abs(np.sin(alpha)) > 0.7:
    phiqp[3, :] = [0, 0, 0, 0, L2 * np.sin(alpha) * w2, 0, 0, 0, 0]
else:
    phiqp[3, :] = [0, 0, 0, 0, L2 * np.cos(alpha) * w2, 0, 0, 0, 0]
phiqp[4, :] = [0, 0, 0, 0, 0, L4 * np.cos(sigma) * w4, 0, 0, 0]
phiqp[5, :] = [0, 0, -2 * (x3p - x2p), -2 * (y3p - y2p), 0, 0, 2 * (x3p - x2p), 2 * (y3p - y2p), 0]
phiqp[6, :] = [0, 0, 0, 0, 0, 0, 2 * x3p, 2 * y3p, 0]
phiqp[7, :] = [0, 0, 0, 0, 0, 0, 0, 0, L6 * np.cos(tita) * w6]

# Producto de phiqp*qp
phiqpqp = phiqp @ qp

# Solución del sistema phiqqpp=rhsa
rhsa = np.append(-phiqpqp, 0)
qpp = np.linalg.solve(phiq, rhsa)

# Resultados
print(
    "\n-----------\nVelocidades\n-----------\n"
    f"v1x = {qp[0]:.3f} mm/s\n"
    f"v1y = {qp[1]:.3f} mm/s\n"
    f"v2x = {qp[2]:.3f} mm/s\n"
    f"v2y = {qp[3]:.3f} mm/s\n"
    f"w2 = {qp[4]:.3f} rad/s\n"
    f"w4 = {qp[5]:.3f} rad/s\n"
    f"v3x = {qp[6]:.3f} mm/s\n"
    f"v3y = {qp[7]:.3f} mm/s\n"
    f"w6 = {qp[8]:.3f} rad/s\n"
)

print(
    "\n-------------\nAceleraciones\n-------------\n"
    f"a1x = {qpp[0]:.3f} mm/s^2\n"
    f"a1y = {qpp[1]:.3f} mm/s^2\n"
    f"a2x = {qpp[2]:.3f} mm/s^2\n"
    f"a2y = {qpp[3]:.3f} mm/s^2\n"
    f"aang2 = {qpp[4]:.3f} rad/s^2\n"
    f"aang4 = {qpp[5]:.3f} rad/s^2\n"
    f"a3x = {qpp[6]:.3f} mm/s^2\n"
    f"a3y = {qpp[7]:.3f} mm/s^2\n"
    f"aang6 = {qpp[8]:.3f} rad/s^2\n"
)
