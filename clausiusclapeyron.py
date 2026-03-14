import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from scipy import stats
from scipy.optimize import curve_fit
import pandas as pd
plt.close('all')
#Cte de los gases ideales
R=8.314 #bar*L/mol*K
#Variación de entalpía molar para la vaporización del agua (valor tabulado)
DHv=-42000

'''Datos para obtener la entalpía tabulada'''
''' a partir de una regresión lineal'''

Th=np.linspace(0,100,32)+273.16
ph=np.array([
    0.61165, 0.70599, 0.81355, 0.93356, 1.0730,
    1.2282, 1.5990, 2.0647, 2.3393, 2.9858,
    3.7831, 4.2470, 5.3251, 6.6328, 7.3849,
    9.1124, 11.177, 12.352, 15.022, 18.171,
    19.946, 23.943, 28.599, 31.201, 37.009,
    43.703, 55.635, 65.017, 70.182, 81.541,
    94.390, 101.42])

resh = stats.linregress(1/Th,np.log(ph))
pendienteh = resh.slope
#Calculo la entalpía a partir del ajuste lineal
Delta_Htab=pendienteh*R
interseccion = resh.intercept
r_pearsonh = resh.rvalue
error_pendienteh = resh.stderr
error_htab=error_pendienteh*R
error_interseccion = resh.intercept_stderr

'''datos = np.loadtxt('archivo.txt')'''

# Puedes separar las columnas fácilmente
'''x = datos[:, 0]
y = datos[:, 1]'''

#Defino una función de desviación relativa para el apartado 2

def desvrel(x,y):
    a=np.abs(x-y)
    return np.abs(a/y)

pv = np.array([0.19,0.25,0.31,0.36,0.42,0.48,0.54,0.61,
               0.67,0.73,0.79,0.83,0.99])
T = np.array([59.47,66.3,70.57,73.32,77.41,81.75,85.82,88.60,
              91.31,92.09,94.36,95.67,100.77])+273.16

#Para LaTeX:
data1={
      r"$p_v (bar)$":pv,
      r"$T (K)$": T
      }

df=pd.DataFrame(data1)

print(df.to_latex(
    index=False,
    caption="Datos experimentales de presión y temperatura según el Método Experimental 1",
    label="tab:cuestion1",
    column_format='|c|c|', 
    position='H',
    bold_rows=True,
    escape=False,
    float_format="%.2f" #Los dos decimales
    ))

#Defino variables para el ajuste lineal
lnpv=np.log(pv)
invT=1/T
#He importado una función más compacta que min_cuadrados, lo siento Samu.
#Se ve que devuelve como argumentos los parámetros que pide el guion

res = stats.linregress(invT,lnpv)
pendiente = res.slope
#Calculo la entalpía a partir del ajuste lineal
Delta_Hlin=pendiente*R
interseccion = res.intercept
r_pearson = res.rvalue
error_pendiente = res.stderr
error_h=error_pendiente*R
error_interseccion = res.intercept_stderr
dnl=desvrel(Delta_Hlin, DHv)
print(f"{'='*30}")
print("RESULTADOS AJUSTE LINEAL")
print(f"{'='*30}")
print(f"Pendiente (m): {pendiente:.2f} ± {error_pendiente:.2f}")
print(f"ΔHv:       {Delta_Hlin:.2f} ± {error_h:.2f} J/mol")
print('Desviación relativa:', dnl)
print(f"Ordenada(b): {interseccion:.2f} ± {error_interseccion:.2f}")
print(f"Coef. Pearson(r): {r_pearson:.6f}")
print(f"Coef. de Determinación(\u03c7\u00b2): {r_pearson**2:.6f}")

#Las seis cifras significativas son de motu proprio,
#he puesto lo que me ha parecido adecuado para al menos
#ver que difiere ligeramente de uno y que no son iguales

#\u03c7\ub002 es para que aparezca chi^2 con letras griegas 
#en el terminal, sin más.

#N.B.: Para escribir ± se puede usar \u00b1

plt.figure()
plt.scatter(pv,T,label='Datos experimentales', color='red')
plt.errorbar(pv, T, 
             xerr=0.01, 
             yerr=0.01, 
             fmt='none',           # 'o' mantiene el estilo de puntos (scatter)
             color='red', 
             ecolor='red',    # Color de las barras de error
             elinewidth=1,      # Grosor de las barras
             capsize=3)
plt.xlabel('Presión de vapor ($bar$)',fontsize=14)
plt.ylabel('Temperatura ($K$)',fontsize=14)
plt.legend()
plt.show()

plt.figure()
plt.scatter(invT,lnpv,label='Datos experimentales', color='red')
plt.errorbar(invT, lnpv, 
             xerr=0.00001, 
             yerr=0.01/T**2, 
             fmt='none',           # 'o' mantiene el estilo de puntos (scatter)
             color='red', 
             ecolor='red',    # Color de las barras de error
             elinewidth=1,      # Grosor de las barras
             capsize=3         # Pestañas en los extremos de las barras
)
plt.ylabel('ln(pv)\n$[p_v]=bar$ ',fontsize=14)
linspace=np.linspace(invT[0],invT[-1],1000)
plt.xlabel('1/T\n$[T]=K$',fontsize=14)
plt.plot(linspace,interseccion+pendiente*linspace,color='blue',
         label='Ajuste lineal')
plt.legend()
plt.show()

'''AJUSTE NO LINEAL'''




#La función que he importado de scipy para el ajuste NO lineal
#funciona por iteraciones, y necesita un valor donde empezar
#más o menos cercano al que queremos obtener para la entalpía, y
#de ahí comenzar a iterar para obtener la gráfica y los datos.

def ecuacion_3(T, p0, delta_H, T0):
    return p0 * np.exp(-(delta_H / R) * (1/T - 1/T0))

# 3. Ajuste no lineal
# p0_guess: valores iniciales aproximados para que el algoritmo no se pierda
p0_guess = [np.mean(pv), np.abs(Delta_Hlin), np.mean(T)] 
popt, pcov = curve_fit(ecuacion_3, T, pv, p0=p0_guess)

# Extraemos los resultados
p0_fit, delta_H_fit, T0_fit = popt
#LA FUNCIÓN EXTRAE CON popt LOS VALORES CON MENOR ERROR SEGÚN pcov
errores = np.sqrt(np.diag(pcov)) # Errores estándar

# 4. Cálculo de Chi (R²) manual para ajustes no lineales
residuos = pv - ecuacion_3(T, *popt)
ss_res = np.sum(residuos**2)
ss_tot = np.sum((pv - np.mean(pv))**2)
chi_r2 = 1 - (ss_res / ss_tot)

dl=desvrel(-delta_H_fit,DHv)
print(f"{'='*30}")
print("RESULTADOS AJUSTE NO LINEAL")
print(f"{'='*30}")
print(f"p0:        {p0_fit:.2f} ± {errores[0]:.2f}")
print(f"ΔHv:       {delta_H_fit:.2f} ± {errores[1]:.2f} J/mol")
print("Desviación relativa:", dl)
print(f"T0:        {T0_fit:.2f} ± {errores[2]:.2f} K")
print(f"Coef. de Determinación(\u03c7\u00b2):   {chi_r2:.6f}")
print(f"{'='*30}")


# 6. Gráfica
plt.figure()
plt.scatter(T, pv, label='Datos Experimentales', color='red')
T_fino = np.linspace(min(T), max(T), 100)
plt.plot(T_fino, ecuacion_3(T_fino, *popt), label='Ajuste según la ec.(3)')
plt.errorbar(T, pv, 
             xerr=0.01, 
             yerr=0.01, 
             fmt='none',           # 'o' mantiene el estilo de puntos (scatter)
             color='red', 
             ecolor='red',    # Color de las barras de error
             elinewidth=1,      # Grosor de las barras
             capsize=3,         # Pestañas en los extremos de las barras
             label='Barras de error')
plt.xlabel('Temperatura ($K$)',fontsize=14)
plt.ylabel('Presión de Vapor ($p_v$)',fontsize=14)
plt.legend()
plt.show()



'''Método 2'''

W=np.array([4.5,3.9,3.47,2.99,2.50])*np.array([62,52,45,40,33])
dmdt=(np.array([145.4,131.7,122.1,113.8,106.8])-98.8)/(388)

res2 = stats.linregress(dmdt,W)
dhv2 = res2.slope*18.2
#Calculo la entalpía a partir del ajuste lineal
qppunto = res2.intercept
r_pearson2 = res2.rvalue
error_pendiente2 = res2.stderr
error_interseccion2 = res2.intercept_stderr
dnl=desvrel(dhv2, DHv)
s=np.sqrt(1*np.array([4.5,3.9,3.47,2.99,2.5])**2+(0.01*np.array([62,52,45,4,33]))**2)

linspace2=np.linspace(dmdt[0],dmdt[-1],1000)
plt.figure()
plt.scatter(dmdt,W,label='Datos experimentales', color='red')
plt.plot(linspace2,qppunto+res2.slope*linspace2,color='blue',label='Ajuste lineal')
plt.errorbar(dmdt, W, 
             xerr=0.005, 
             yerr=0, 
             fmt='none',           # 'o' mantiene el estilo de puntos (scatter)
             color='red', 
             ecolor='red',    # Color de las barras de error
             elinewidth=1,      # Grosor de las barras
             capsize=3)
plt.ylabel('Potencia ($W$)',fontsize=14)
plt.xlabel('Pérdida de masa por u. de tiempo ($kg s^{-1}$)',fontsize=14)
plt.legend()
plt.show()

dr2=desvrel(dhv2,DHv)
print("RESULTADOS METODO 2")
print(f"{'='*30}")
print(f"Pendiente (m): {dhv2:.2f} ± {error_pendiente2:.2f}")
print('Desviación relativa:', dnl)
print(f"Ordenada(b): {qppunto:.2f} ± {error_interseccion2:.2f}")
print(f"Coef. Pearson(r): {r_pearson2:.6f}")
print(f"Coef. de Determinación(\u03c7\u00b2): {r_pearson2**2:.6f}")

data2={
      r"$Potencia (W)$":W,
      r"$\frac{dm}{dt} (kg s^{-1})$": dmdt
      }

dg=pd.DataFrame(data2)

print(dg.to_latex(
    index=False,
    caption="Datos experimentales de potencia eléctrica y pérdida de masa por unidad de tiempo",
    label="tab:cuestion2",
    column_format='|c|c|', 
    position='H',
    bold_rows=True,
    escape=False,
    float_format="%.2f"
    ))
print(dhv2)





