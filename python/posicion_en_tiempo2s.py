import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend("Agg")

# Trabajo Teoría de Mecanismos
# Posición en 2s

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
sigma = np.deg2rad(120.5891)
ro = np.deg2rad(115.117)
tita = np.deg2rad(180 - 74.4079)
delta = np.deg2rad(9.5245)

# Datos de la simulación
duracion = 2
At = 0.01
npasos = int(duracion / At)
t = np.zeros(npasos)
avance = np.zeros(npasos)
avancez2 = np.zeros(npasos)

# Bucle principal
for i in range(npasos):
    # Problema de posición
    alpha += w2 * At
    error = 10000.0
    niter = 0

    while error > 1e-6 and niter < 50:
        # Restricciones
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

        # Jacobiano
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

        # Aq
        Aq = np.linalg.solve(phiq, np.append(-phi, 0))

        # Actualizar coordenadas
        x1 += Aq[0]
        y1 += Aq[1]
        x2 += Aq[2]
        y2 += Aq[3]
        sigma += Aq[5]
        x3 += Aq[6]
        y3 += Aq[7]
        tita += Aq[8]

        # Calcular error
        error = np.linalg.norm(Aq)
        niter += 1

    plt.clf()
    plt.plot([xA, x1], [yA, y1], color="r", linewidth=4, marker="o", markersize=6)
    plt.plot([x2, xB], [y2, yB], color="k", linewidth=4, marker="o", markersize=6)
    plt.plot([x2, x3], [y2, y3], color="k", linewidth=4, marker="o", markersize=6)
    plt.plot([x3, xC], [y3, yC], color="k", linewidth=4, marker="o", markersize=6)
    plt.xlim([-200, 300])
    plt.ylim([-150, 200])
    plt.pause(1e-10)
    t[i] = i * At
    avance[i] = y2
    avancez2[i] = x2

plt.figure()
plt.plot(t, avance, t, avancez2, "--")
plt.title("Avance/retroceso de herramienta")
plt.xlabel("tiempo (s)")
plt.ylabel("Coordenada X [- -](mm), Y [-](mm)")

print(
    "\n----------\nPosiciones\n----------\n"
    f"x1 = {x1:.3f} mm\n"
    f"y1 = {y1:.3f} mm\n"
    f"x2 = {x2:.3f} mm\n"
    f"y2 = {y2:.3f} mm\n"
    f"x3 = {x3:.3f} mm\n"
    f"y3 = {y3:.3f} mm\n"
    f"alpha = {alpha:.3f} rad\n"
    f"sigma = {sigma:.3f} rad\n"
    f"tita = {tita:.3f} rad\n"
)
