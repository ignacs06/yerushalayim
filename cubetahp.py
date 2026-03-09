# =============================================================================
# Formulas Practica 1: Cubeta de Ondas
# =============================================================================
#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import t


def min_cuadrados_avanzado(x, y, t_student, x_0):
    n = len(x)
    x_prom = np.mean(x)
    y_prom = np.mean(y)

    D = np.sum((x - x_prom)**2)
    E = np.sum((x - x_prom) * (y - y_prom))

    m = E / D
    c = y_prom - m * x_prom

    y_predicha = m * x + c
    s_res_2 = np.sum((y - y_predicha)**2) / (n - 2)

    y_0 = m * x_0 + c
    error_y0 = t_student * np.sqrt(s_res_2 * (1 / n + (x_0 - x_prom)**2 / D))

    return y_0, error_y0, m, c


# --- Ejemplo de uso con tus datos ---
# Supongamos que T y S son tus arrays de temperatura y tension
# T = np.array([...])
# S = np.array([...])
# ts1 = valor de la t de Student

# =============================================================================
# Calculo de sigma mediante regresion (NO USAR)
# =============================================================================
# ts1 = t.ppf(1 - 0.01 / 2, 11)
# T = np.array([16.20,16.90,17.18,17.45,19.61,19.63,20.00,24.00,24.80,24.90,25.00,25.40])
# S = np.array([73.41,73.06,73.17,73.06,72.78,72.75,72.75,71.50,72.00,72.00,72.00,72.00])
# te = np.linspace(16.00,26.00,1000)
# temp_objetivo = 20.5
# tension, error, pendiente, ordenada = min_cuadrados_avanzado(T, S, ts1, temp_objetivo)
# =============================================================================

# Introduce aqui sigma y su incertidumbre medidas en laboratorio (mN/m)
tension = 72.6
error = 4.2

h = 8.19e-3
ih = 0.01e-3
g = 9.81
rho = 1000

# Datos de Excel (foto)
f = np.array([47, 52, 57, 61, 67, 75])
delt_fr = np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
lambdap = np.array([5.60, 5.20, 4.63, 4.61, 4.22, 3.95]) * 1e-3
delt_lambda = np.array([1, 1, 1, 1, 1, 1], dtype=float) * 1e-4
lambdac = np.array([5.88, 5.48, 5.27, 4.82, 4.42, 4.12]) * 1e-3

v1p = lambdap * f
iv1p = np.sqrt((f * delt_lambda)**2 + (lambdap * delt_fr)**2)
v1c = lambdac * f
iv1c = np.sqrt((f * delt_lambda)**2 + (lambdac * delt_fr)**2)

v3 = np.sqrt(g * h)
iv3 = np.sqrt(g * ih / (2 * v3))


def v4(l):
    # Ecuacion (4): v = sqrt((g*lambda)/(2*pi) + (2*pi*sigma)/(rho*lambda))
    sigma = tension / 1000.0
    dsigma = error / 1000.0

    termino1 = (g * l) / (2 * np.pi)
    termino2 = (2 * np.pi * sigma) / (rho * l)
    valor_v4 = np.sqrt(termino1 + termino2)

    # Propagacion: (dv)^2 = (dv/dlambda dlambda)^2 + (dv/dsigma dsigma)^2
    dv_dlambda = ((g / (2 * np.pi)) - (2 * np.pi * sigma) / (rho * l**2)) / (2 * valor_v4)
    dv_dsigma = (2 * np.pi / (rho * l)) / (2 * valor_v4)
    incertidumbre = np.sqrt((dv_dlambda * delt_lambda)**2 + (dv_dsigma * dsigma)**2)

    return valor_v4, incertidumbre


v4p, iv4p = v4(lambdap)
v4c, iv4c = v4(lambdac)


plt.figure()
plt.tick_params(axis='both', labelsize=20)

plt.errorbar(lambdap*1e3, v1p*1e2, yerr=iv1p*1e2, xerr=delt_lambda*1e3, fmt='o', color='red',
             ecolor='red', capsize=5, label='_nolegend_')
plt.scatter(lambdap*1e3, v1p*1e2, color='red', label=r'$v_{\mathrm{ec.(1)}}$')

plt.errorbar(lambdap*1e3, v4p*1e2, yerr=iv4p*1e2, xerr=delt_lambda*1e3, fmt='o', color='blue',
             ecolor='blue', capsize=5, label='_nolegend_')
plt.scatter(lambdap*1e3, v4p*1e2, color='blue', label=r'$v_{\mathrm{ec.\,(4)}}$')
plt.xlabel(r'$\lambda\;(\mathrm{mm})$', fontsize=25)
plt.ylabel(r'$v\;(\mathrm{cm\,s^{-1}})$', fontsize=25)
plt.grid('True')
plt.legend(fontsize=20)
plt.gca().set_facecolor('aliceblue')
plt.show()

# Comparacion con valores leidos del Excel (foto)
v1p_excel = np.array([0.26320, 0.27040, 0.263910, 0.281210, 0.28274, 0.29625])
v1c_excel = np.array([0.27636, 0.28496, 0.30039, 0.29402, 0.29614, 0.30900])

print("v1p calculada:", np.round(v1p, 5))
print("v1p excel   :", np.round(v1p_excel, 5))
print("dif v1p     :", np.round(v1p - v1p_excel, 5))
print("v1c calculada:", np.round(v1c, 5))
print("v1c excel   :", np.round(v1c_excel, 5))
print("dif v1c     :", np.round(v1c - v1c_excel, 5))


def v5(l):
    termino1 = (g * l) / (2 * np.pi)
    termino2 = (2 * np.pi * tension) / (1000 * rho * l)
    valor_v4 = np.sqrt(termino1 + termino2)
    return valor_v4


def v_grupo(l):
    # v_g = (1/(2*w)) * (g + 3*k^2*sigma/rho), con k=2*pi/lambda
    sigma = tension / 1000.0
    k = 2 * np.pi / l
    omega = np.sqrt(g * k + (sigma / rho) * k**3)
    return (g + 3 * k**2 * sigma / rho) / (2 * omega)


def punto_corte_numerico(x, y1, y2):
    diff = y1 - y2
    idx = np.where(np.diff(np.sign(diff)) != 0)[0]

    if len(idx) > 0:
        i = idx[0]
        x1, x2 = x[i], x[i + 1]
        d1, d2 = diff[i], diff[i + 1]
        xc = x1 - d1 * (x2 - x1) / (d2 - d1)
        yc = y1[i] + (y1[i + 1] - y1[i]) * (xc - x1) / (x2 - x1)
        return xc, yc

    i = np.argmin(np.abs(diff))
    return x[i], y1[i]


lambdal = np.linspace(3e-3, 20e-3, 1000)
v_fase_l = v5(lambdal)
v_grupo_l = v_grupo(lambdal)
lambda_corte, v_corte = punto_corte_numerico(lambdal, v_fase_l, v_grupo_l)

plt.figure()
plt.tick_params(axis='both', labelsize=20)
plt.plot(lambdal*1e3, v_fase_l*1e2, color='red', label=r'$v_{\mathrm{ec.\,(4)}}$',zorder=2)
plt.plot(lambdal*1e3, v_grupo_l*1e2, color='blue', label=r'$v_{\mathrm{grupo}}$',zorder=2)
plt.scatter(lambda_corte*1e3, v_corte*1e2, color='purple', s=70, label='Punto de corte',zorder=10)
plt.xlabel(r'$\lambda\;(\mathrm{mm})$', fontsize=25)
plt.ylabel(r'$v\;(\mathrm{cm\,s^{-1}})$', fontsize=25)
plt.grid('True')
plt.legend(fontsize=20)
plt.gca().set_facecolor('aliceblue')

plt.show()

def punto_corte_analitico(sigma=tension/1000.0, rho=1000, g=9.81):
    lambda_corte = 2 * np.pi * np.sqrt(sigma / (rho * g))
    v_corte,delt_v_corte = v4(lambda_corte)
    delt_lambda_corte = (np.pi / (rho*g*np.sqrt(sigma / (rho * g)))) * (error / 1000.0)
    return lambda_corte, v_corte, delt_lambda_corte, delt_v_corte  
print("Punto de corte analitico:", punto_corte_analitico()[0]*1e3,'$\pm$',punto_corte_analitico()[2]*1e3, punto_corte_analitico()[1]*1e2, '$\pm$', punto_corte_analitico()[3][0]*1e2)
print("Punto de corte numerico:", lambda_corte*1e3, v_corte*1e2)

def angulos_interferencia(lamda, D, N):
    # Minimos de interferencia de dos fuentes coherentes:
    # D * sin(theta_N) = (2N+1) * lambda / 2
    seno = (2 * N + 1) * lamda / (2 * D)
    theta_deg = np.full_like(seno, np.nan, dtype=float)
    delt_theta_deg = np.full_like(seno, np.nan, dtype=float)
    validos = np.abs(seno) <= 1
    theta_deg[validos] = np.degrees(np.arcsin(seno[validos]))

    # Propagacion de incertidumbre en radianes y conversion final a grados
    # Delta lambda = 0.1 mm, Delta D = 0.2 cm
    delt_lambda = 0.1e-3  # m
    delt_D = 0.1e-2       # m
    den = np.sqrt(1 - seno[validos]**2)
    dtheta_dlambda = (2 * N[validos] + 1) / (2 * D * den)
    dtheta_dD = -((2 * N[validos] + 1) * lamda) / (2 * D**2 * den)
    delt_theta_rad = np.sqrt((dtheta_dlambda * delt_lambda)**2 + (dtheta_dD * delt_D)**2)
    delt_theta_deg[validos] = np.degrees(delt_theta_rad)

    return theta_deg, delt_theta_deg


lambda_67 = 4.42e-3
lambda_47 = 5.9e-3
D_1 = 2.3e-2
D_2 = 3.0e-2
N = np.array([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])

theta_67_1, delt_theta_67_1 = angulos_interferencia(lambda_67, D_1, N)
theta_47, delt_theta_47 = angulos_interferencia(lambda_47, D_1, N)
theta_67_2, delt_theta_67_2 = angulos_interferencia(lambda_67, D_2, N)

print("N:", N)
print("theta teoricos 67Hz, D1 (deg):", np.round(theta_67_1, 2))
print("theta teoricos 47Hz, D1 (deg):", np.round(theta_47, 2))
print("theta teoricos 67Hz, D2 (deg):", np.round(theta_67_2, 2))
print("delt_theta 67Hz, D1 (deg):", np.round(delt_theta_67_1, 2))
print("delt_theta 47Hz, D1 (deg):", np.round(delt_theta_47, 2))
print("delt_theta 67Hz, D2 (deg):", np.round(delt_theta_67_2, 2))

v4p=v4p*1e2
iv4p=iv4p*1e2

# %%
