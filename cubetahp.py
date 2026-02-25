# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal.
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import t

def min_cuadrados_avanzado(x, y, t, x_0):
    n = len(x)
    x_prom = np.mean(x)
    y_prom = np.mean(y)
    
    # Cálculos intermedios (Vectorizados)
    D = np.sum((x - x_prom)**2)
    E = np.sum((x - x_prom) * (y - y_prom))
    
    # Pendiente y ordenada
    m = E / D
    c = y_prom - m * x_prom
    
    # Cálculo de residuos y varianza residual
    y_predicha = m * x + c
    s_res_2 = np.sum((y - y_predicha)**2) / (n - 2)
    
    # 1. Valor buscado a la temperatura x_0
    y_0 = m * x_0 + c
    
    # 2. Error de la predicción en x_0 (Fórmula de laboratorio)
    # Esta fórmula tiene en cuenta que el error es mínimo en el promedio de x
    # y aumenta a medida que te alejas de los datos medidos.
    error_y0 = t * np.sqrt(s_res_2 * (1/n + (x_0 - x_prom)**2 / D))
    
    # Errores de m y c (por si los sigues necesitando)

    return y_0, error_y0, m, c

# --- Ejemplo de uso con tus datos ---
# Supongamos que T y S son tus arrays de temperatura y tensión
# T = np.array([...])
# S = np.array([...])
# ts1 = valor de la t de Student

ts1=t.ppf(1-0.01/2,11)
T=np.array([16.20,16.90,17.18,17.45,19.61,19.63,20.00,24.00,24.80,24.90,25.00,25.40])
S=np.array([73.41,73.06,73.17,73.06,72.78,72.75,72.75,71.50,72.00,72.00,72.00,72.00])
te=np.linspace(16.00,26.00,1000)
h=10.87e-3
ih=1e-5
g=9.80
rho=1000

temp_objetivo = 20.5
tension, error, pendiente, ordenada = min_cuadrados_avanzado(T, S, ts1, temp_objetivo)
lambdac=np.array([6,5.4,4.9,4.6,4.2,3.9])*1e-3
lambdap=np.array([6.1,5.4,4.8,4.6,4.1,3.8])*1e-3
ilambdas=1e-5
f=np.array([40,47,54,61,68,75])
ifr=1

v1p=lambdap*f
iv1p=np.sqrt((f*ilambdas)**2+(lambdap*ifr)**2)
v1c=lambdac*f
iv1c=np.sqrt((f*ilambdas)**2+(lambdac*ifr)**2)

v3=np.sqrt(g*h)
iv3=np.sqrt(g*ih/(2*v3))



def v4(l):
    # Cálculo de la velocidad v4
    # Fórmula: sqrt( (g*l / 2pi) + (2pi * tension / (rho * l)) )
    termino1 = (g * l) / (2 * np.pi)
    termino2 = (2 * np.pi * tension) / (1000*rho * l)
    valor_v4 = np.sqrt(termino1 + termino2)
    
    # Cálculo de la incertidumbre iv4
    # He separado los términos para que sea legible y evitar errores de sintaxis
    parte_a = (g / (2 * np.pi)) - (2 * np.pi * (tension/1000)) / (rho * l**2)
    incertidumbre = np.sqrt( (parte_a**2 * ilambdas**2) + (2 * np.pi * error/1000 / (rho * l))**2 ) / (2 * valor_v4)
    
    return valor_v4, incertidumbre
v4p,iv4p=v4(lambdap)
v4c,iv4c=v4(lambdac)


plt.figure()
plt.grid(True)
plt.plot(lambdap,v1p,'o',color='red',label='v=$\lambda*f$ para ondas planas')
plt.errorbar(lambdap,v1p,yerr=iv1p,xerr=ilambdas,fmt='o',color='red',ecolor='red',capsize=5)

plt.plot(lambdap,v1c,'o',color='blue',label='v=$\lambda*f$ para ondas circulares')
plt.errorbar(lambdap,v1c,yerr=iv1c,xerr=ilambdas,fmt='o',color='blue',ecolor='blue',capsize=5)

plt.plot(lambdap,v4p,'o',color='lightgreen',label='$v_{prof}$ para ondas planas')
plt.errorbar(lambdap,v4p,yerr=iv4p,xerr=ilambdas,fmt='o',color='lightgreen',ecolor='lightgreen',capsize=5)

plt.plot(lambdap,v4c,'o',color='green',label='$v_{prof}$ para ondas circulares')
plt.errorbar(lambdap,v4c,yerr=iv4c,xerr=ilambdas,fmt='o',color='green',ecolor='green',capsize=5)

plt.plot(lambdap,np.array([v3,v3,v3,v3,v3,v3]),color='orange',label='$v_{sup}$')
plt.errorbar(lambdap,np.array([v3,v3,v3,v3,v3,v3]),yerr=iv3,fmt='o',color='orange',ecolor='orange',capsize=5)



def vg(l):
    k=2*np.pi/l
    vg1=g/(2*np.sqrt(g*k+tension*k**3/(1000*rho)))
    vg2=3*k**2*tension/(1000*rho*(2*np.sqrt(g*k+tension*k**3/(1000*rho))))
    return vg1+vg2
vgp=vg(lambdap)
vgc=vg(lambdac)

def v5(l):
    # Cálculo de la velocidad v4
    # Fórmula: sqrt( (g*l / 2pi) + (2pi * tension / (rho * l)) )
    termino1 = (g * l) / (2 * np.pi)
    termino2 = (2 * np.pi * tension) / (1000*rho * l)
    valor_v4 = np.sqrt(termino1 + termino2)
    
    return valor_v4



lambdal=np.linspace(3e-3,20e-3,1000)
plt.figure()
plt.grid(True)
plt.plot(lambdal,vg(lambdal),label='$v_g$')
plt.plot(lambdal,v5(lambdal),label='$v_{fase}$')
plt.legend()









