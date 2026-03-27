## RESULTADOS PARA LAS ISOTERMAS DE UN GAS REAL:

#%%

## ASPECTOS A TENER EN CUENTA:
# 1) Definir todas las temperaturas y presiones de coexistencia así como los errores instrumentales.
# 2) Comprobar que las raíces de las ecuaciones VdW y RK están en el orden esperado de pequeña a grande.
# 3) Comprobar que las unidades y los decimales son los correctos en los resultados.
# 4) Falta añadir punto crítico experimental


import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
plt.close('all')

## TEMPERATURAS EN GRADOS CELSIUS:

T1=21
T2=25.5
T3=29.6
T4=34
T5=38.4
T6=46

## CARGA Y VISUALIZACIÓN DE LOS DATOS:

datos_original=pd.read_excel('Isotermas.xlsx',sheet_name=None)
print(datos_original.keys())

isoterma_1=datos_original['Isoterma 1']
isoterma_2=datos_original['Isoterma 2']
isoterma_3=datos_original['Isoterma 3']
isoterma_4=datos_original['Isoterma 4']
isoterma_5=datos_original['Isoterma 5']
isoterma_6=datos_original['Isoterma 6']

print(isoterma_1.head(3))
print(isoterma_2.head(3))
print(isoterma_3.head(3))
print(isoterma_4.head(3))
print(isoterma_5.head(3))
print(isoterma_6.head(3))

datos_ajuste=pd.read_excel('Isotermas_ajuste.xlsx',sheet_name=None)
print(datos_ajuste.keys())

isoterma_1_ajuste=datos_ajuste['Isoterma 1']
isoterma_2_ajuste=datos_ajuste['Isoterma 2']
isoterma_3_ajuste=datos_ajuste['Isoterma 3']
isoterma_4_ajuste=datos_ajuste['Isoterma 4']
isoterma_5_ajuste=datos_ajuste['Isoterma 5']
isoterma_6_ajuste=datos_ajuste['Isoterma 6']

print(f"Número de filas en isoterma_1_ajuste: {len(isoterma_1_ajuste)}")
print(f"Número de filas en isoterma_2_ajuste: {len(isoterma_2_ajuste)}")
print(f"Número de filas en isoterma_3_ajuste: {len(isoterma_3_ajuste)}")
print(f"Número de filas en isoterma_4_ajuste: {len(isoterma_4_ajuste)}")
print(f"Número de filas en isoterma_5_ajuste: {len(isoterma_5_ajuste)}")
print(f"Número de filas en isoterma_6_ajuste: {len(isoterma_6_ajuste)}")

## CONVERSIÓN DE TEMPERATURAS A KELVIN:

T1=T1+273.15
T2=T2+273.15
T3=T3+273.15
T4=T4+273.15
T5=T5+273.15
T6=T6+273.15

print(f"Temperatura T1 en Kelvin: {T1:.2f} K",)
print(f"Temperatura T2 en Kelvin: {T2:.2f} K",)
print(f"Temperatura T3 en Kelvin: {T3:.2f} K",)
print(f"Temperatura T4 en Kelvin: {T4:.2f} K",)
print(f"Temperatura T5 en Kelvin: {T5:.2f} K",)
print(f"Temperatura T6 en Kelvin: {T6:.2f} K",)

## GRÁFICOS DE LAS ISOTERMAS:

delt_P=0.5
delt_V=0.05

V1=isoterma_1.iloc[:,0].values
P1=isoterma_1.iloc[:,1].values

V2=isoterma_2.iloc[:,0].values
P2=isoterma_2.iloc[:,1].values

V3=isoterma_3.iloc[:,0].values
P3=isoterma_3.iloc[:,1].values

V4=isoterma_4.iloc[:,0].values
P4=isoterma_4.iloc[:,1].values

V5=isoterma_5.iloc[:,0].values
P5=isoterma_5.iloc[:,1].values

V6=isoterma_6.iloc[:,0].values
P6=isoterma_6.iloc[:,1].values

print(f"V1: {V1}")
print(f"P1: {P1}")

plt.figure()
plt.tick_params(axis='both', labelsize=20)

plt.errorbar(V1,P1,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='blue',ecolor='blue',capsize=3,label='_nolegend_')
plt.scatter(V1,P1,color='blue',label=f"Isoterma $T1={T1:.1f}~K$")

plt.errorbar(V2,P2,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='green',ecolor='green',capsize=3,label='_nolegend_')
plt.scatter(V2,P2,color='green',label=f"Isoterma $T2={T2:.1f}~K$")

plt.errorbar(V3,P3,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='red',ecolor='red',capsize=3,label='_nolegend_')
plt.scatter(V3,P3,color='red',label=f"Isoterma $T3={T3:.1f}~K$")

plt.errorbar(V4,P4,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='cyan',ecolor='cyan',capsize=3,label='_nolegend_')
plt.scatter(V4,P4,color='cyan',label=f"Isoterma $T4={T4:.1f}~K$")

plt.errorbar(V5,P5,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='magenta',ecolor='magenta',capsize=3,label='_nolegend_')
plt.scatter(V5,P5,color='magenta',label=f"Isoterma $T5={T5:.1f}~K$")

plt.errorbar(V6,P6,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='orange',ecolor='orange',capsize=3,label='_nolegend_')
plt.scatter(V6,P6,color='orange',label=f"Isoterma $T6={T6:.1f}~K$")

plt.title('Medidas experimentales', fontsize=25)
plt.xlabel(r'$V\;(\mathrm{mL})$', fontsize=25)
plt.ylabel(r'$P\;(\mathrm{bar})$', fontsize=25)
plt.legend(fontsize=12, ncol=2, loc='upper right')
plt.gca().set_facecolor('aliceblue')
plt.show()

## DEFINICIÓN DE FUNCIONES:

T=None

def isoterma_Van_der_Waals(V,A,B,C):
    if T is None:
        raise ValueError("La temperatura T no ha sido definida.")
    return (C*T)/(V-B)-A/(V**2)

def Redlich_Kwong(V,A,B,C):
    if T is None:
        raise ValueError("La temperatura T no ha sido definida.")
    return (C*T)/(V-B)-A/(V*(V+B)*np.sqrt(T))

## AJUSTE DE LAS ISOTERMAS:

# Datos del ajuste:

V1_ajuste=isoterma_1_ajuste.iloc[:,0].values
P1_ajuste=isoterma_1_ajuste.iloc[:,1].values    

V2_ajuste=isoterma_2_ajuste.iloc[:,0].values
P2_ajuste=isoterma_2_ajuste.iloc[:,1].values

V3_ajuste=isoterma_3_ajuste.iloc[:,0].values
P3_ajuste=isoterma_3_ajuste.iloc[:,1].values

V4_ajuste=isoterma_4_ajuste.iloc[:,0].values
P4_ajuste=isoterma_4_ajuste.iloc[:,1].values

V5_ajuste=isoterma_5_ajuste.iloc[:,0].values
P5_ajuste=isoterma_5_ajuste.iloc[:,1].values

V6_ajuste=isoterma_6_ajuste.iloc[:,0].values
P6_ajuste=isoterma_6_ajuste.iloc[:,1].values

print(f"V1_ajuste: {V1_ajuste}")
print(f"P1_ajuste: {P1_ajuste}")

V0 = np.median(V1_ajuste)
P0 = np.median(P1_ajuste)
C0 = (P0*V0)/T1
A0 = P0*V0**2
B0 = 0.05*V0
p0 = [A0, B0, C0]


# Ajuste para la isoterma 1 (T1):

T=T1 # Ajuste
params,cov=curve_fit(isoterma_Van_der_Waals,V1_ajuste,P1_ajuste,p0=p0,bounds = ([0, 0, 0], [np.inf, 0.9*np.min(V1_ajuste), np.inf]),maxfev=100000)

A1,B1,C1=params #Parámetros
delt_A1=np.sqrt(np.diag(cov))[0]
delt_B1=np.sqrt(np.diag(cov))[1]
delt_C1=np.sqrt(np.diag(cov))[2]
print(f"Parámetros ajustados para T1: A={A1:.3e} ± {delt_A1:.3e} bar mL², B={B1:.3e} ± {delt_B1:.3e} mL, C={C1:.3e} ± {delt_C1:.3e} bar mL/K")

print(f"Covarianza para T1:\n{cov}") #Matriz de covarianza

print(f"La incertidumbre relativa en A1 es: {delt_A1/A1:.2%}") # Incertidumbre relativa
print(f"La incertidumbre relativa en B1 es: {delt_B1/B1:.2%}")
print(f"La incertidumbre relativa en C1 es: {delt_C1/C1:.2%}")

residuos1=P1_ajuste-isoterma_Van_der_Waals(V1_ajuste,A1,B1,C1) # Bondad del ajuste
ss_res1=np.sum(residuos1**2)
ss_tot1=np.sum((P1_ajuste-np.mean(P1_ajuste))**2)
R2_1=1-(ss_res1/ss_tot1) # Coeficiente de determinación
N1=len(P1_ajuste)
p1=3
gl_1=N1-p1
R2_1_ajustado=1-(1-R2_1)*(N1-1)/gl_1 # R² ajustado
sigma_1=delt_P
P1_ajuste_teorico=isoterma_Van_der_Waals(V1_ajuste,A1,B1,C1)
chi_1_ajuste=np.sum(((P1_ajuste-P1_ajuste_teorico)/sigma_1)**2) # Chi-cuadrado
chi_1_reducido=chi_1_ajuste/gl_1 # Chi-cuadrado reducido
sigma_modelo_ajuste_1=np.sqrt(np.sum(residuos1**2)/gl_1) # Error estándar del modelo ajustado

print(f"Número de datos del ajuste para T1: {N1}")
print(f"R² para T1: {R2_1:.4f}")
print(f"R² ajustado para T1: {R2_1_ajustado:.4f}")
print(f"Chi-cuadrado para T1: {chi_1_ajuste:.3f}")
print(f"Chi-cuadrado reducido para T1: {chi_1_reducido:.3f}")
print(f"Error estándar del modelo ajustado para T1: {sigma_modelo_ajuste_1:.3e} bar")

P1_coex=21.5 # Cálculo de las raíces de van der Waals para la isoterma 1 (T1):
coefs_1=np.array([P1_coex,-(P1_coex*B1+C1*T1),A1,-A1*B1])
raices_1=np.roots(coefs_1)
raices_reales_1=raices_1[np.isreal(raices_1)].real
print(f"Volúmenes reales a presión {P1_coex} bar: {raices_reales_1}")

V1_prueba=np.linspace(min(V1_ajuste),max(V1_ajuste),100)# Definición de intervalos del ajuste
P1_prueba=isoterma_Van_der_Waals(V1_prueba,A1,B1,C1)

plt.figure()
plt.tick_params(axis='both', labelsize=20)

plt.errorbar(V1,P1,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='blue',ecolor='blue',capsize=3,label='_nolegend_')
plt.scatter(V1,P1,color='blue',label=f"Datos $T1={T1:.1f}~K$")
plt.plot(V1_prueba,P1_prueba,color='blue',linestyle='-',linewidth=2,label=f"Ajuste VdW T1")

plt.xlabel(r'$V\;(\mathrm{mL})$', fontsize=25)
plt.ylabel(r'$P\;(\mathrm{bar})$', fontsize=25)
plt.legend(fontsize=12, ncol=2, loc='upper right')
plt.gca().set_facecolor('aliceblue')
plt.show()

V1_liq=np.linspace(min(V1_ajuste),raices_reales_1[2],100)# Definición de intervalos del ajuste
P1_liq=isoterma_Van_der_Waals(V1_liq,A1,B1,C1)

V1_coex=np.linspace(raices_reales_1[2],raices_reales_1[0],100)
P1_coex=np.full_like(V1_coex,P1_coex)

V1_vapor=np.linspace(raices_reales_1[0],max(V1_ajuste),100)
P1_vapor=isoterma_Van_der_Waals(V1_vapor,A1,B1,C1)

#%%

# Ajuste para la isoterma 2 (T2):

T=T2
params,cov=curve_fit(isoterma_Van_der_Waals,V2_ajuste,P2_ajuste,p0=p0,bounds=([0, 0, 0], [np.inf, 0.5*0.9, np.inf]))

A2,B2,C2=params
delt_A2=np.sqrt(np.diag(cov))[0]
delt_B2=np.sqrt(np.diag(cov))[1]
delt_C2=np.sqrt(np.diag(cov))[2]
print(f"Parámetros ajustados para T2: A={A2:.3e} ± {delt_A2:.3e} bar mL², B={B2:.3e} ± {delt_B2:.3e} mL, C={C2:.3e} ± {delt_C2:.3e} bar mL/K")

print(f"Covarianza para T2:\n{cov}")

print(f"La incertidumbre relativa en A2 es: {delt_A2/A2:.2%}")
print(f"La incertidumbre relativa en B2 es: {delt_B2/B2:.2%}")
print(f"La incertidumbre relativa en C2 es: {delt_C2/C2:.2%}")

residuos2=P2_ajuste-isoterma_Van_der_Waals(V2_ajuste,A2,B2,C2)
ss_res2=np.sum(residuos2**2)
ss_tot2=np.sum((P2_ajuste-np.mean(P2_ajuste))**2)
R2_2=1-(ss_res2/ss_tot2)
N2=len(P2_ajuste)
p2=3
gl_2=N2-p2
R2_2_ajustado=1-(1-R2_2)*(N2-1)/gl_2
sigma_2=delt_P
P2_ajuste_teorico=isoterma_Van_der_Waals(V2_ajuste,A2,B2,C2)
chi_2_ajuste=np.sum(((P2_ajuste-P2_ajuste_teorico)/sigma_2)**2)
chi_2_reducido=chi_2_ajuste/gl_2
sigma_modelo_ajuste_2=np.sqrt(np.sum(residuos2**2)/gl_2)

print(f"Número de datos del ajuste para T2: {N2}")
print(f"R² para T2: {R2_2:.4f}")
print(f"R² ajustado para T2: {R2_2_ajustado:.4f}")
print(f"Chi-cuadrado para T2: {chi_2_ajuste:.3f}")
print(f"Chi-cuadrado reducido para T2: {chi_2_reducido:.3f}")
print(f"Error estándar del modelo ajustado para T2: {sigma_modelo_ajuste_2:.3e} bar")

P2_coex=24
coefs_2=np.array([P2_coex,-(P2_coex*B2+C2*T2),A2,-A2*B2])
raices_2=np.roots(coefs_2)
raices_reales_2=raices_2[np.isreal(raices_2)].real
print(f"Volúmenes reales a presión {P2_coex} bar: {raices_reales_2}")

V2_liq=np.linspace(min(V2_ajuste),raices_reales_2[2],100)
P2_liq=isoterma_Van_der_Waals(V2_liq,A2,B2,C2)

V2_coex=np.linspace(raices_reales_2[2],raices_reales_2[0],100)
P2_coex=np.full_like(V2_coex,P2_coex)

V2_vapor=np.linspace(raices_reales_2[0],max(V2_ajuste),100)
P2_vapor=isoterma_Van_der_Waals(V2_vapor,A2,B2,C2)

#%%

# Ajuste para la isoterma 3 (T3):

T=T3
params,cov=curve_fit(isoterma_Van_der_Waals,V3_ajuste,P3_ajuste,p0=p0,bounds=([0, 0, 0], [np.inf, 0.5*0.9, np.inf]))

A3,B3,C3=params
delt_A3=np.sqrt(np.diag(cov))[0]
delt_B3=np.sqrt(np.diag(cov))[1]
delt_C3=np.sqrt(np.diag(cov))[2]
print(f"Parámetros ajustados para T3: A={A3:.3e} ± {delt_A3:.3e} bar mL², B={B3:.3e} ± {delt_B3:.3e} mL, C={C3:.3e} ± {delt_C3:.3e} bar mL/K")

print(f"Covarianza para T3:\n{cov}")

print(f"La incertidumbre relativa en A3 es: {delt_A3/A3:.2%}")
print(f"La incertidumbre relativa en B3 es: {delt_B3/B3:.2%}")
print(f"La incertidumbre relativa en C3 es: {delt_C3/C3:.2%}")

residuos3=P3_ajuste-isoterma_Van_der_Waals(V3_ajuste,A3,B3,C3)
ss_res3=np.sum(residuos3**2)
ss_tot3=np.sum((P3_ajuste-np.mean(P3_ajuste))**2)
R2_3=1-(ss_res3/ss_tot3)
N3=len(P3_ajuste)
p3=3
gl_3=N3-p3
R2_3_ajustado=1-(1-R2_3)*(N3-1)/gl_3
sigma_3=delt_P
P3_ajuste_teorico=isoterma_Van_der_Waals(V3_ajuste,A3,B3,C3)
chi_3_ajuste=np.sum(((P3_ajuste-P3_ajuste_teorico)/sigma_3)**2)
chi_3_reducido=chi_3_ajuste/gl_3
sigma_modelo_ajuste_3=np.sqrt(np.sum(residuos3**2)/gl_3)

print(f"Número de datos del ajuste para T3: {N3}")
print(f"R² para T3: {R2_3:.4f}")
print(f"R² ajustado para T3: {R2_3_ajustado:.4f}")
print(f"Chi-cuadrado para T3: {chi_3_ajuste:.3f}")
print(f"Chi-cuadrado reducido para T3: {chi_3_reducido:.3f}")
print(f"Error estándar del modelo ajustado para T3: {sigma_modelo_ajuste_3:.3e} bar")

P3_coex=26.5
coefs_3=np.array([P3_coex,-(P3_coex*B3+C3*T3),A3,-A3*B3])
raices_3=np.roots(coefs_3)
raices_reales_3=raices_3[np.isreal(raices_3)].real
print(f"Volúmenes reales a presión {P3_coex} bar: {raices_reales_3}")

V3_liq=np.linspace(min(V3_ajuste),raices_reales_3[2],100)
P3_liq=isoterma_Van_der_Waals(V3_liq,A3,B3,C3)

V3_coex=np.linspace(raices_reales_3[2],raices_reales_3[0],100)
P3_coex=np.full_like(V3_coex,P3_coex)

V3_vapor=np.linspace(raices_reales_3[0],max(V3_ajuste),100)
P3_vapor=isoterma_Van_der_Waals(V3_vapor,A3,B3,C3)

#%%

# Ajuste para la isoterma 4 (T4):

T=T4
params,cov=curve_fit(isoterma_Van_der_Waals,V4_ajuste,P4_ajuste,p0=p0,bounds=([0, 0, 0], [np.inf, 0.5*0.9, np.inf]))

A4,B4,C4=params
delt_A4=np.sqrt(np.diag(cov))[0]
delt_B4=np.sqrt(np.diag(cov))[1]
delt_C4=np.sqrt(np.diag(cov))[2]
print(f"Parámetros ajustados para T4: A={A4:.3e} ± {delt_A4:.3e} bar mL², B={B4:.3e} ± {delt_B4:.3e} mL, C={C4:.3e} ± {delt_C4:.3e} bar mL/K")

print(f"Covarianza para T4:\n{cov}")

print(f"La incertidumbre relativa en A4 es: {delt_A4/A4:.2%}")
print(f"La incertidumbre relativa en B4 es: {delt_B4/B4:.2%}")
print(f"La incertidumbre relativa en C4 es: {delt_C4/C4:.2%}")

residuos4=P4_ajuste-isoterma_Van_der_Waals(V4_ajuste,A4,B4,C4)
ss_res4=np.sum(residuos4**2)
ss_tot4=np.sum((P4_ajuste-np.mean(P4_ajuste))**2)
R2_4=1-(ss_res4/ss_tot4)
N4=len(P4_ajuste)
p4=3
gl_4=N4-p4
R2_4_ajustado=1-(1-R2_4)*(N4-1)/gl_4
sigma_4=delt_P
P4_ajuste_teorico=isoterma_Van_der_Waals(V4_ajuste,A4,B4,C4)
chi_4_ajuste=np.sum(((P4_ajuste-P4_ajuste_teorico)/sigma_4)**2)
chi_4_reducido=chi_4_ajuste/gl_4
sigma_modelo_ajuste_4=np.sqrt(np.sum(residuos4**2)/gl_4)

print(f"Número de datos del ajuste para T4: {N4}")
print(f"R² para T4: {R2_4:.4f}")
print(f"R² ajustado para T4: {R2_4_ajustado:.4f}")
print(f"Chi-cuadrado para T4: {chi_4_ajuste:.3f}")
print(f"Chi-cuadrado reducido para T4: {chi_4_reducido:.3f}")
print(f"Error estándar del modelo ajustado para T4: {sigma_modelo_ajuste_4:.3e} bar")

P4_coex=29.5
coefs_4=np.array([P4_coex,-(P4_coex*B4+C4*T4),A4,-A4*B4])
raices_4=np.roots(coefs_4)
raices_reales_4=raices_4[np.isreal(raices_4)].real
print(f"Volúmenes reales a presión {P4_coex} bar: {raices_reales_4}")

V4_liq=np.linspace(min(V4_ajuste),raices_reales_4[2],100)
P4_liq=isoterma_Van_der_Waals(V4_liq,A4,B4,C4)

V4_coex=np.linspace(raices_reales_4[2],raices_reales_4[0],100)
P4_coex=np.full_like(V4_coex,P4_coex)

V4_vapor=np.linspace(raices_reales_4[0],max(V4_ajuste),100)
P4_vapor=isoterma_Van_der_Waals(V4_vapor,A4,B4,C4)

#%%

# Ajuste para la isoterma 5 (T5):

T=T5
params,cov=curve_fit(isoterma_Van_der_Waals,V5_ajuste,P5_ajuste,p0=p0,bounds=([0, 0, 0], [np.inf, 0.5*0.9, np.inf]))

A5,B5,C5=params
delt_A5=np.sqrt(np.diag(cov))[0]
delt_B5=np.sqrt(np.diag(cov))[1]
delt_C5=np.sqrt(np.diag(cov))[2]
print(f"Parámetros ajustados para T5: A={A5:.3e} ± {delt_A5:.3e} bar mL², B={B5:.3e} ± {delt_B5:.3e} mL, C={C5:.3e} ± {delt_C5:.3e} bar mL/K")

print(f"Covarianza para T5:\n{cov}")

print(f"La incertidumbre relativa en A5 es: {delt_A5/A5:.2%}")
print(f"La incertidumbre relativa en B5 es: {delt_B5/B5:.2%}")
print(f"La incertidumbre relativa en C5 es: {delt_C5/C5:.2%}")

residuos5=P5_ajuste-isoterma_Van_der_Waals(V5_ajuste,A5,B5,C5)
ss_res5=np.sum(residuos5**2)
ss_tot5=np.sum((P5_ajuste-np.mean(P5_ajuste))**2)
R2_5=1-(ss_res5/ss_tot5)
N5=len(P5_ajuste)
p5=3
gl_5=N5-p5
R2_5_ajustado=1-(1-R2_5)*(N5-1)/gl_5
sigma_5=delt_P
P5_ajuste_teorico=isoterma_Van_der_Waals(V5_ajuste,A5,B5,C5)
chi_5_ajuste=np.sum(((P5_ajuste-P5_ajuste_teorico)/sigma_5)**2)
chi_5_reducido=chi_5_ajuste/gl_5
sigma_modelo_ajuste_5=np.sqrt(np.sum(residuos5**2)/gl_5)

print(f"Número de datos del ajuste para T5: {N5}")
print(f"R² para T5: {R2_5:.4f}")
print(f"R² ajustado para T5: {R2_5_ajustado:.4f}")
print(f"Chi-cuadrado para T5: {chi_5_ajuste:.3f}")
print(f"Chi-cuadrado reducido para T5: {chi_5_reducido:.3f}")
print(f"Error estándar del modelo ajustado para T5: {sigma_modelo_ajuste_5:.3e} bar")

P5_coex=32.5
coefs_5=np.array([P5_coex,-(P5_coex*B5+C5*T5),A5,-A5*B5])
raices_5=np.roots(coefs_5)
raices_reales_5=raices_5[np.isreal(raices_5)].real
print(f"Volúmenes reales a presión {P5_coex} bar: {raices_reales_5}")

V5_liq=np.linspace(min(V5_ajuste),raices_reales_5[2],100)
P5_liq=isoterma_Van_der_Waals(V5_liq,A5,B5,C5)

V5_coex=np.linspace(raices_reales_5[2],raices_reales_5[0],100)
P5_coex=np.full_like(V5_coex,P5_coex)

V5_vapor=np.linspace(raices_reales_5[0],max(V5_ajuste),100)
P5_vapor=isoterma_Van_der_Waals(V5_vapor,A5,B5,C5)

#%%

# Ajuste para la isoterma 6 (T6):

T=T6
params,cov=curve_fit(isoterma_Van_der_Waals,V6_ajuste,P6_ajuste,p0=p0,bounds=([0, 0, 0], [np.inf, 0.5*0.9, np.inf]))

A6,B6,C6=params
delt_A6=np.sqrt(np.diag(cov))[0]
delt_B6=np.sqrt(np.diag(cov))[1]
delt_C6=np.sqrt(np.diag(cov))[2]
print(f"Parámetros ajustados para T6: A={A6:.3e} ± {delt_A6:.3e} bar mL², B={B6:.3e} ± {delt_B6:.3e} mL, C={C6:.3e} ± {delt_C6:.3e} bar mL/K")

print(f"Covarianza para T6:\n{cov}")

print(f"La incertidumbre relativa en A6 es: {delt_A6/A6:.2%}")
print(f"La incertidumbre relativa en B6 es: {delt_B6/B6:.2%}")
print(f"La incertidumbre relativa en C6 es: {delt_C6/C6:.2%}")

residuos6=P6_ajuste-isoterma_Van_der_Waals(V6_ajuste,A6,B6,C6)
ss_res6=np.sum(residuos6**2)
ss_tot6=np.sum((P6_ajuste-np.mean(P6_ajuste))**2)
R2_6=1-(ss_res6/ss_tot6)
N6=len(P6_ajuste)
p6=3
gl_6=N6-p6
R2_6_ajustado=1-(1-R2_6)*(N6-1)/gl_6
sigma_6=delt_P
P6_ajuste_teorico=isoterma_Van_der_Waals(V6_ajuste,A6,B6,C6)
chi_6_ajuste=np.sum(((P6_ajuste-P6_ajuste_teorico)/sigma_6)**2)
chi_6_reducido=chi_6_ajuste/gl_6
sigma_modelo_ajuste_6=np.sqrt(np.sum(residuos6**2)/gl_6)

print(f"Número de datos del ajuste para T6: {N6}")
print(f"R² para T6: {R2_6:.4f}")
print(f"R² ajustado para T6: {R2_6_ajustado:.4f}")
print(f"Chi-cuadrado para T6: {chi_6_ajuste:.3f}")
print(f"Chi-cuadrado reducido para T6: {chi_6_reducido:.3f}")
print(f"Error estándar del modelo ajustado para T6: {sigma_modelo_ajuste_6:.3e} bar")

V6_total=np.linspace(min(V6_ajuste),max(V6_ajuste),100)
P6_total=isoterma_Van_der_Waals(V6_total,A6,B6,C6)

plt.figure()
plt.tick_params(axis='both', labelsize=20)

plt.errorbar(V1,P1,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='blue',ecolor='blue',capsize=3,label='_nolegend_')
plt.scatter(V1,P1,color='blue',label=f"Datos $T1={T1:.1f}~K$")
plt.plot(V1_liq,P1_liq,color='blue',linestyle='-',linewidth=2,label=f"Ajuste VdW líquido T1")
plt.plot(V1_coex,P1_coex,color='blue',linestyle='-',linewidth=2,label=f"Coexistencia T1")
plt.plot(V1_vapor,P1_vapor,color='blue',linestyle='-',linewidth=2,label=f"Ajuste VdW vapor T1")

plt.errorbar(V2,P2,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='green',ecolor='green',capsize=3,label='_nolegend_')
plt.scatter(V2,P2,color='green',label=f"Datos $T2={T2:.1f}~K$")
plt.plot(V2_liq,P2_liq,color='green',linestyle='-',linewidth=2,label=f"Ajuste VdW líquido T2")
plt.plot(V2_coex,P2_coex,color='green',linestyle='-',linewidth=2,label=f"Coexistencia T2")
plt.plot(V2_vapor,P2_vapor,color='green',linestyle='-',linewidth=2,label=f"Ajuste VdW vapor T2")

plt.errorbar(V3,P3,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='red',ecolor='red',capsize=3,label='_nolegend_')
plt.scatter(V3,P3,color='red',label=f"Datos $T3={T3:.1f}~K$")
plt.plot(V3_liq,P3_liq,color='red',linestyle='-',linewidth=2,label=f"Ajuste VdW líquido T3")
plt.plot(V3_coex,P3_coex,color='red',linestyle='-',linewidth=2,label=f"Coexistencia T3")
plt.plot(V3_vapor,P3_vapor,color='red',linestyle='-',linewidth=2,label=f"Ajuste VdW vapor T3")

plt.errorbar(V4,P4,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='cyan',ecolor='cyan',capsize=3,label='_nolegend_')
plt.scatter(V4,P4,color='cyan',label=f"Datos $T4={T4:.1f}~K$")
plt.plot(V4_liq,P4_liq,color='cyan',linestyle='-',linewidth=2,label=f"Ajuste VdW líquido T4")
plt.plot(V4_coex,P4_coex,color='cyan',linestyle='-',linewidth=2,label=f"Coexistencia T4")
plt.plot(V4_vapor,P4_vapor,color='cyan',linestyle='-',linewidth=2,label=f"Ajuste VdW vapor T4")

plt.errorbar(V5,P5,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='magenta',ecolor='magenta',capsize=3,label='_nolegend_')
plt.scatter(V5,P5,color='magenta',label=f"Datos $T5={T5:.1f}~K$")
plt.plot(V5_liq,P5_liq,color='magenta',linestyle='-',linewidth=2,label=f"Ajuste VdW líquido T5")
plt.plot(V5_coex,P5_coex,color='magenta',linestyle='-',linewidth=2,label=f"Coexistencia T5")
plt.plot(V5_vapor,P5_vapor,color='magenta',linestyle='-',linewidth=2,label=f"Ajuste VdW vapor T5")

plt.errorbar(V6,P6,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='orange',ecolor='orange',capsize=3,label='_nolegend_')
plt.scatter(V6,P6,color='orange',label=f"Datos $T6={T6:.1f}~K$")
plt.plot(V6_total,P6_total,color='orange',linestyle='-',linewidth=2,label=f"Ajuste VdW T6")

plt.title('Ajustes de Van der Waals para las 6 isotermas', fontsize=25)
plt.xlabel(r'$V\;(\mathrm{mL})$', fontsize=25)
plt.ylabel(r'$P\;(\mathrm{bar})$', fontsize=25)
plt.legend(fontsize=12, ncol=2, loc='upper right')
plt.gca().set_facecolor('aliceblue')
plt.show()

#%%

## AJUSTE REDLICH-KWONG:

# Ajuste para la isoterma 1 (T1):

T=T1
params_prima,cov_prima=curve_fit(Redlich_Kwong,V1_ajuste,P1_ajuste,p0=p0)

A1_prima,B1_prima,C1_prima=params_prima
delt_A1_prima=np.sqrt(np.diag(cov_prima))[0]
delt_B1_prima=np.sqrt(np.diag(cov_prima))[1]
delt_C1_prima=np.sqrt(np.diag(cov_prima))[2]
print(f"Parámetros ajustados RK para T1: A={A1_prima:.3e} ± {delt_A1_prima:.3e} bar mL² B={B1_prima:.3e} ± {delt_B1_prima:.3e} mL, C={C1_prima:.3e} ± {delt_C1_prima:.3e} bar mL/K")

print(f"Covarianza RK para T1:\n{cov_prima}")

print(f"La incertidumbre relativa en A1_prima es: {delt_A1_prima/A1_prima:.2%}")
print(f"La incertidumbre relativa en B1_prima es: {delt_B1_prima/B1_prima:.2%}")
print(f"La incertidumbre relativa en C1_prima es: {delt_C1_prima/C1_prima:.2%}")

residuos1_prima=P1_ajuste-Redlich_Kwong(V1_ajuste,A1_prima,B1_prima,C1_prima)
ss_res1_prima=np.sum(residuos1_prima**2)
ss_tot1_prima=np.sum((P1_ajuste-np.mean(P1_ajuste))**2)
R2_1_prima=1-(ss_res1_prima/ss_tot1_prima)
N1_prima=len(P1_ajuste)
p1_prima=3
gl_1_prima=N1_prima-p1_prima
R2_1_ajustado_prima=1-(1-R2_1_prima)*(N1_prima-1)/gl_1_prima
sigma_1_prima=delt_P
P1_ajuste_teorico_prima=Redlich_Kwong(V1_ajuste,A1_prima,B1_prima,C1_prima)
chi_1_ajuste_prima=np.sum(((P1_ajuste-P1_ajuste_teorico_prima)/sigma_1_prima)**2)
chi_1_reducido_prima=chi_1_ajuste_prima/gl_1_prima
sigma_modelo_ajuste_1_prima=np.sqrt(np.sum(residuos1_prima**2)/gl_1_prima)

print(f"Número de datos del ajuste RK para T1: {N1_prima}")
print(f"R² RK para T1: {R2_1_prima:.4f}")
print(f"R² ajustado RK para T1: {R2_1_ajustado_prima:.4f}")
print(f"Chi-cuadrado RK para T1: {chi_1_ajuste_prima:.3f}")
print(f"Chi-cuadrado reducido RK para T1: {chi_1_reducido_prima:.3f}")
print(f"Error estándar del modelo ajustado RK para T1: {sigma_modelo_ajuste_1_prima:.3e} bar")

P1_coex=21.5
coefs_1_prima=np.array([P1_coex,-(C1_prima*T1),(-P1_coex*(B1_prima**2)-C1_prima*T1*B1_prima + A1_prima/np.sqrt(T1)),(-A1_prima*B1_prima/np.sqrt(T1))])
raices_1_prima=np.roots(coefs_1_prima)
raices_reales_1_prima=raices_1_prima[np.isreal(raices_1_prima)].real
print(f"Volúmenes reales RK a presión {P1_coex} bar: {raices_reales_1_prima}")

V1_liq_prima=np.linspace(min(V1_ajuste),raices_reales_1_prima[2],100)
P1_liq_prima=Redlich_Kwong(V1_liq_prima,A1_prima,B1_prima,C1_prima)

V1_coex_prima=np.linspace(raices_reales_1_prima[2],raices_reales_1_prima[0],100)
P1_coex_prima=np.full_like(V1_coex_prima,P1_coex)

V1_vapor_prima=np.linspace(raices_reales_1_prima[0],max(V1_ajuste),100)
P1_vapor_prima=Redlich_Kwong(V1_vapor_prima,A1_prima,B1_prima,C1_prima)

# Ajuste para la isoterma 2 (T2):

T=T2
params_prima,cov_prima=curve_fit(Redlich_Kwong,V2_ajuste,P2_ajuste,p0=p0)

A2_prima,B2_prima,C2_prima=params_prima
delt_A2_prima=np.sqrt(np.diag(cov_prima))[0]
delt_B2_prima=np.sqrt(np.diag(cov_prima))[1]
delt_C2_prima=np.sqrt(np.diag(cov_prima))[2]
print(f"Parámetros ajustados RK para T2: A={A2_prima:.3e} ± {delt_A2_prima:.3e} bar mL², B={B2_prima:.3e} ± {delt_B2_prima:.3e} mL, C={C2_prima:.3e} ± {delt_C2_prima:.3e} bar mL/K")

print(f"Covarianza RK para T2:\n{cov_prima}")

print(f"La incertidumbre relativa en A2_prima es: {delt_A2_prima/A2_prima:.2%}")
print(f"La incertidumbre relativa en B2_prima es: {delt_B2_prima/B2_prima:.2%}")
print(f"La incertidumbre relativa en C2_prima es: {delt_C2_prima/C2_prima:.2%}")

residuos2_prima=P2_ajuste-Redlich_Kwong(V2_ajuste,A2_prima,B2_prima,C2_prima)
ss_res2_prima=np.sum(residuos2_prima**2)
ss_tot2_prima=np.sum((P2_ajuste-np.mean(P2_ajuste))**2)
R2_2_prima=1-(ss_res2_prima/ss_tot2_prima)
N2_prima=len(P2_ajuste)
p2_prima=3
gl_2_prima=N2_prima-p2_prima
R2_2_ajustado_prima=1-(1-R2_2_prima)*(N2_prima-1)/gl_2_prima
sigma_2_prima=delt_P
P2_ajuste_teorico_prima=Redlich_Kwong(V2_ajuste,A2_prima,B2_prima,C2_prima)
chi_2_ajuste_prima=np.sum(((P2_ajuste-P2_ajuste_teorico_prima)/sigma_2_prima)**2)
chi_2_reducido_prima=chi_2_ajuste_prima/gl_2_prima
sigma_modelo_ajuste_2_prima=np.sqrt(np.sum(residuos2_prima**2)/gl_2_prima)

print(f"Número de datos del ajuste RK para T2: {N2_prima}")
print(f"R² RK para T2: {R2_2_prima:.4f}")
print(f"R² ajustado RK para T2: {R2_2_ajustado_prima:.4f}")
print(f"Chi-cuadrado RK para T2: {chi_2_ajuste_prima:.3f}")
print(f"Chi-cuadrado reducido RK para T2: {chi_2_reducido_prima:.3f}")
print(f"Error estándar del modelo ajustado RK para T2: {sigma_modelo_ajuste_2_prima:.3e} bar")

P2_coex=24
coefs_2_prima=np.array([P2_coex,-(C2_prima*T2),(-P2_coex*(B2_prima**2)-C2_prima*T2*B2_prima + A2_prima/np.sqrt(T2)),(-A2_prima*B2_prima/np.sqrt(T2))])
raices_2_prima=np.roots(coefs_2_prima)
raices_reales_2_prima=raices_2_prima[np.isreal(raices_2_prima)].real
print(f"Volúmenes reales RK a presión {P2_coex} bar: {raices_reales_2_prima}")

V2_liq_prima=np.linspace(min(V2_ajuste),raices_reales_2_prima[2],100)
P2_liq_prima=Redlich_Kwong(V2_liq_prima,A2_prima,B2_prima,C2_prima)

V2_coex_prima=np.linspace(raices_reales_2_prima[2],raices_reales_2_prima[0],100)
P2_coex_prima=np.full_like(V2_coex_prima,P2_coex)

V2_vapor_prima=np.linspace(raices_reales_2_prima[0],max(V2_ajuste),100)
P2_vapor_prima=Redlich_Kwong(V2_vapor_prima,A2_prima,B2_prima,C2_prima)

# Ajuste para la isoterma 3 (T3):

T=T3
params_prima,cov_prima=curve_fit(Redlich_Kwong,V3_ajuste,P3_ajuste,p0=p0)

A3_prima,B3_prima,C3_prima=params_prima
delt_A3_prima=np.sqrt(np.diag(cov_prima))[0]
delt_B3_prima=np.sqrt(np.diag(cov_prima))[1]
delt_C3_prima=np.sqrt(np.diag(cov_prima))[2]
print(f"Parámetros ajustados RK para T3: A={A3_prima:.3e} ± {delt_A3_prima:.3e} bar mL², B={B3_prima:.3e} ± {delt_B3_prima:.3e} mL, C={C3_prima:.3e} ± {delt_C3_prima:.3e} bar mL/K")

print(f"Covarianza RK para T3:\n{cov_prima}")

print(f"La incertidumbre relativa en A3_prima es: {delt_A3_prima/A3_prima:.2%}")
print(f"La incertidumbre relativa en B3_prima es: {delt_B3_prima/B3_prima:.2%}")
print(f"La incertidumbre relativa en C3_prima es: {delt_C3_prima/C3_prima:.2%}")

residuos3_prima=P3_ajuste-Redlich_Kwong(V3_ajuste,A3_prima,B3_prima,C3_prima)
ss_res3_prima=np.sum(residuos3_prima**2)
ss_tot3_prima=np.sum((P3_ajuste-np.mean(P3_ajuste))**2)
R2_3_prima=1-(ss_res3_prima/ss_tot3_prima)
N3_prima=len(P3_ajuste)
p3_prima=3
gl_3_prima=N3_prima-p3_prima
R2_3_ajustado_prima=1-(1-R2_3_prima)*(N3_prima-1)/gl_3_prima
sigma_3_prima=delt_P
P3_ajuste_teorico_prima=Redlich_Kwong(V3_ajuste,A3_prima,B3_prima,C3_prima)
chi_3_ajuste_prima=np.sum(((P3_ajuste-P3_ajuste_teorico_prima)/sigma_3_prima)**2)
chi_3_reducido_prima=chi_3_ajuste_prima/gl_3_prima
sigma_modelo_ajuste_3_prima=np.sqrt(np.sum(residuos3_prima**2)/gl_3_prima)

print(f"Número de datos del ajuste RK para T3: {N3_prima}")
print(f"R² RK para T3: {R2_3_prima:.4f}")
print(f"R² ajustado RK para T3: {R2_3_ajustado_prima:.4f}")
print(f"Chi-cuadrado RK para T3: {chi_3_ajuste_prima:.3f}")
print(f"Chi-cuadrado reducido RK para T3: {chi_3_reducido_prima:.3f}")
print(f"Error estándar del modelo ajustado RK para T3: {sigma_modelo_ajuste_3_prima:.3e} bar")

P3_coex=26.5
coefs_3_prima=np.array([P3_coex,-(C3_prima*T3),(-P3_coex*(B3_prima**2)-C3_prima*T3*B3_prima + A3_prima/np.sqrt(T3)),(-A3_prima*B3_prima/np.sqrt(T3))])
raices_3_prima=np.roots(coefs_3_prima)
raices_reales_3_prima=raices_3_prima[np.isreal(raices_3_prima)].real
print(f"Volúmenes reales RK a presión {P3_coex} bar: {raices_reales_3_prima}")

V3_liq_prima=np.linspace(min(V3_ajuste),raices_reales_3_prima[2],100)
P3_liq_prima=Redlich_Kwong(V3_liq_prima,A3_prima,B3_prima,C3_prima)

V3_coex_prima=np.linspace(raices_reales_3_prima[2],raices_reales_3_prima[0],100)
P3_coex_prima=np.full_like(V3_coex_prima,P3_coex)

V3_vapor_prima=np.linspace(raices_reales_3_prima[0],max(V3_ajuste),100)
P3_vapor_prima=Redlich_Kwong(V3_vapor_prima,A3_prima,B3_prima,C3_prima)

# Ajuste para la isoterma 4 (T4):

T=T4
params_prima,cov_prima=curve_fit(Redlich_Kwong,V4_ajuste,P4_ajuste,p0=p0)

A4_prima,B4_prima,C4_prima=params_prima
delt_A4_prima=np.sqrt(np.diag(cov_prima))[0]
delt_B4_prima=np.sqrt(np.diag(cov_prima))[1]
delt_C4_prima=np.sqrt(np.diag(cov_prima))[2]
print(f"Parámetros ajustados RK para T4: A={A4_prima:.3e} ± {delt_A4_prima:.3e} bar mL², B={B4_prima:.3e} ± {delt_B4_prima:.3e} mL, C={C4_prima:.3e} ± {delt_C4_prima:.3e} bar mL/K")

print(f"Covarianza RK para T4:\n{cov_prima}")

print(f"La incertidumbre relativa en A4_prima es: {delt_A4_prima/A4_prima:.2%}")
print(f"La incertidumbre relativa en B4_prima es: {delt_B4_prima/B4_prima:.2%}")
print(f"La incertidumbre relativa en C4_prima es: {delt_C4_prima/C4_prima:.2%}")

residuos4_prima=P4_ajuste-Redlich_Kwong(V4_ajuste,A4_prima,B4_prima,C4_prima)
ss_res4_prima=np.sum(residuos4_prima**2)
ss_tot4_prima=np.sum((P4_ajuste-np.mean(P4_ajuste))**2)
R2_4_prima=1-(ss_res4_prima/ss_tot4_prima)
N4_prima=len(P4_ajuste)
p4_prima=3
gl_4_prima=N4_prima-p4_prima
R2_4_ajustado_prima=1-(1-R2_4_prima)*(N4_prima-1)/gl_4_prima
sigma_4_prima=delt_P
P4_ajuste_teorico_prima=Redlich_Kwong(V4_ajuste,A4_prima,B4_prima,C4_prima)
chi_4_ajuste_prima=np.sum(((P4_ajuste-P4_ajuste_teorico_prima)/sigma_4_prima)**2)
chi_4_reducido_prima=chi_4_ajuste_prima/gl_4_prima
sigma_modelo_ajuste_4_prima=np.sqrt(np.sum(residuos4_prima**2)/gl_4_prima)

print(f"Número de datos del ajuste RK para T4: {N4_prima}")
print(f"R² RK para T4: {R2_4_prima:.4f}")
print(f"R² ajustado RK para T4: {R2_4_ajustado_prima:.4f}")
print(f"Chi-cuadrado RK para T4: {chi_4_ajuste_prima:.3f}")
print(f"Chi-cuadrado reducido RK para T4: {chi_4_reducido_prima:.3f}")
print(f"Error estándar del modelo ajustado RK para T4: {sigma_modelo_ajuste_4_prima:.3e} bar")

P4_coex=29.5
coefs_4_prima=np.array([P4_coex,-(C4_prima*T4),(-P4_coex*(B4_prima**2)-C4_prima*T4*B4_prima + A4_prima/np.sqrt(T4)),(-A4_prima*B4_prima/np.sqrt(T4))])
raices_4_prima=np.roots(coefs_4_prima)
raices_reales_4_prima=raices_4_prima[np.isreal(raices_4_prima)].real
print(f"Volúmenes reales RK a presión {P4_coex} bar: {raices_reales_4_prima}")

V4_liq_prima=np.linspace(min(V4_ajuste),raices_reales_4_prima[2],100)
P4_liq_prima=Redlich_Kwong(V4_liq_prima,A4_prima,B4_prima,C4_prima)

V4_coex_prima=np.linspace(raices_reales_4_prima[2],raices_reales_4_prima[0],100)
P4_coex_prima=np.full_like(V4_coex_prima,P4_coex)

V4_vapor_prima=np.linspace(raices_reales_4_prima[0],max(V4_ajuste),100)
P4_vapor_prima=Redlich_Kwong(V4_vapor_prima,A4_prima,B4_prima,C4_prima)

# Ajuste para la isoterma 5 (T5):

T=T5
params_prima,cov_prima=curve_fit(Redlich_Kwong,V5_ajuste,P5_ajuste,p0=p0)

A5_prima,B5_prima,C5_prima=params_prima
delt_A5_prima=np.sqrt(np.diag(cov_prima))[0]
delt_B5_prima=np.sqrt(np.diag(cov_prima))[1]
delt_C5_prima=np.sqrt(np.diag(cov_prima))[2]
print(f"Parámetros ajustados RK para T5: A={A5_prima:.3e} ± {delt_A5_prima:.3e} bar mL², B={B5_prima:.3e} ± {delt_B5_prima:.3e} mL, C={C5_prima:.3e} ± {delt_C5_prima:.3e} bar mL/K")

print(f"Covarianza RK para T5:\n{cov_prima}")

print(f"La incertidumbre relativa en A5_prima es: {delt_A5_prima/A5_prima:.2%}")
print(f"La incertidumbre relativa en B5_prima es: {delt_B5_prima/B5_prima:.2%}")
print(f"La incertidumbre relativa en C5_prima es: {delt_C5_prima/C5_prima:.2%}")

residuos5_prima=P5_ajuste-Redlich_Kwong(V5_ajuste,A5_prima,B5_prima,C5_prima)
ss_res5_prima=np.sum(residuos5_prima**2)
ss_tot5_prima=np.sum((P5_ajuste-np.mean(P5_ajuste))**2)
R2_5_prima=1-(ss_res5_prima/ss_tot5_prima)
N5_prima=len(P5_ajuste)
p5_prima=3
gl_5_prima=N5_prima-p5_prima
R2_5_ajustado_prima=1-(1-R2_5_prima)*(N5_prima-1)/gl_5_prima
sigma_5_prima=delt_P
P5_ajuste_teorico_prima=Redlich_Kwong(V5_ajuste,A5_prima,B5_prima,C5_prima)
chi_5_ajuste_prima=np.sum(((P5_ajuste-P5_ajuste_teorico_prima)/sigma_5_prima)**2)
chi_5_reducido_prima=chi_5_ajuste_prima/gl_5_prima
sigma_modelo_ajuste_5_prima=np.sqrt(np.sum(residuos5_prima**2)/gl_5_prima)

print(f"Número de datos del ajuste RK para T5: {N5_prima}")
print(f"R² RK para T5: {R2_5_prima:.4f}")
print(f"R² ajustado RK para T5: {R2_5_ajustado_prima:.4f}")
print(f"Chi-cuadrado RK para T5: {chi_5_ajuste_prima:.3f}")
print(f"Chi-cuadrado reducido RK para T5: {chi_5_reducido_prima:.3f}")
print(f"Error estándar del modelo ajustado RK para T5: {sigma_modelo_ajuste_5_prima:.3e} bar")

P5_coex=32.5
coefs_5_prima=np.array([P5_coex,-(C5_prima*T5),(-P5_coex*(B5_prima**2)-C5_prima*T5*B5_prima + A5_prima/np.sqrt(T5)),(-A5_prima*B5_prima/np.sqrt(T5))])
raices_5_prima=np.roots(coefs_5_prima)
raices_reales_5_prima=raices_5_prima[np.isreal(raices_5_prima)].real
print(f"Volúmenes reales RK a presión {P5_coex} bar: {raices_reales_5_prima}")

V5_liq_prima=np.linspace(min(V5_ajuste),raices_reales_5_prima[2],100)
P5_liq_prima=Redlich_Kwong(V5_liq_prima,A5_prima,B5_prima,C5_prima)

V5_coex_prima=np.linspace(raices_reales_5_prima[2],raices_reales_5_prima[0],100)
P5_coex_prima=np.full_like(V5_coex_prima,P5_coex)

V5_vapor_prima=np.linspace(raices_reales_5_prima[0],max(V5_ajuste),100)
P5_vapor_prima=Redlich_Kwong(V5_vapor_prima,A5_prima,B5_prima,C5_prima)

# Ajuste para la isoterma 6 (T6):

T=T6
params_prima,cov_prima=curve_fit(Redlich_Kwong,V6_ajuste,P6_ajuste,p0=p0)

A6_prima,B6_prima,C6_prima=params_prima
delt_A6_prima=np.sqrt(np.diag(cov_prima))[0]
delt_B6_prima=np.sqrt(np.diag(cov_prima))[1]
delt_C6_prima=np.sqrt(np.diag(cov_prima))[2]
print(f"Parámetros ajustados RK para T6: A={A6_prima:.3e} ± {delt_A6_prima:.3e} bar mL², B={B6_prima:.3e} ± {delt_B6_prima:.3e} mL, C={C6_prima:.3e} ± {delt_C6_prima:.3e} bar mL/K")

print(f"Covarianza RK para T6:\n{cov_prima}")

print(f"La incertidumbre relativa en A6_prima es: {delt_A6_prima/A6_prima:.2%}")
print(f"La incertidumbre relativa en B6_prima es: {delt_B6_prima/B6_prima:.2%}")
print(f"La incertidumbre relativa en C6_prima es: {delt_C6_prima/C6_prima:.2%}")

residuos6_prima=P6_ajuste-Redlich_Kwong(V6_ajuste,A6_prima,B6_prima,C6_prima)
ss_res6_prima=np.sum(residuos6_prima**2)
ss_tot6_prima=np.sum((P6_ajuste-np.mean(P6_ajuste))**2)
R2_6_prima=1-(ss_res6_prima/ss_tot6_prima)
N6_prima=len(P6_ajuste)
p6_prima=3
gl_6_prima=N6_prima-p6_prima
R2_6_ajustado_prima=1-(1-R2_6_prima)*(N6_prima-1)/gl_6_prima
sigma_6_prima=delt_P
P6_ajuste_teorico_prima=Redlich_Kwong(V6_ajuste,A6_prima,B6_prima,C6_prima)
chi_6_ajuste_prima=np.sum(((P6_ajuste-P6_ajuste_teorico_prima)/sigma_6_prima)**2)
chi_6_reducido_prima=chi_6_ajuste_prima/gl_6_prima
sigma_modelo_ajuste_6_prima=np.sqrt(np.sum(residuos6_prima**2)/gl_6_prima)

print(f"Número de datos del ajuste RK para T6: {N6_prima}")
print(f"R² RK para T6: {R2_6_prima:.4f}")
print(f"R² ajustado RK para T6: {R2_6_ajustado_prima:.4f}")
print(f"Chi-cuadrado RK para T6: {chi_6_ajuste_prima:.3f}")
print(f"Chi-cuadrado reducido RK para T6: {chi_6_reducido_prima:.3f}")
print(f"Error estándar del modelo ajustado RK para T6: {sigma_modelo_ajuste_6_prima:.3e} bar")

V6_total_prima=np.linspace(min(V6_ajuste),max(V6_ajuste),100)
P6_total_prima=Redlich_Kwong(V6_total_prima,A6_prima,B6_prima,C6_prima)

plt.figure()
plt.tick_params(axis='both', labelsize=20)

plt.errorbar(V1,P1,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='blue',ecolor='blue',capsize=3,label='_nolegend_')
plt.scatter(V1,P1,color='blue',label=f"Datos $T1={T1:.1f}~K$")
plt.plot(V1_liq_prima,P1_liq_prima,color='blue',linestyle='-',linewidth=2,label=f"Ajuste RK lí­quido T1")
plt.plot(V1_coex_prima,P1_coex_prima,color='blue',linestyle='-',linewidth=2,label=f"Coexistencia T1")
plt.plot(V1_vapor_prima,P1_vapor_prima,color='blue',linestyle='-',linewidth=2,label=f"Ajuste RK vapor T1")

plt.errorbar(V2,P2,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='green',ecolor='green',capsize=3,label='_nolegend_')
plt.scatter(V2,P2,color='green',label=f"Datos $T2={T2:.1f}~K$")
plt.plot(V2_liq_prima,P2_liq_prima,color='green',linestyle='-',linewidth=2,label=f"Ajuste RK lí­quido T2")
plt.plot(V2_coex_prima,P2_coex_prima,color='green',linestyle='-',linewidth=2,label=f"Coexistencia T2")
plt.plot(V2_vapor_prima,P2_vapor_prima,color='green',linestyle='-',linewidth=2,label=f"Ajuste RK vapor T2")

plt.errorbar(V3,P3,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='red',ecolor='red',capsize=3,label='_nolegend_')
plt.scatter(V3,P3,color='red',label=f"Datos $T3={T3:.1f}~K$")
plt.plot(V3_liq_prima,P3_liq_prima,color='red',linestyle='-',linewidth=2,label=f"Ajuste RK lí­quido T3")
plt.plot(V3_coex_prima,P3_coex_prima,color='red',linestyle='-',linewidth=2,label=f"Coexistencia T3")
plt.plot(V3_vapor_prima,P3_vapor_prima,color='red',linestyle='-',linewidth=2,label=f"Ajuste RK vapor T3")

plt.errorbar(V4,P4,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='cyan',ecolor='cyan',capsize=3,label='_nolegend_')
plt.scatter(V4,P4,color='cyan',label=f"Datos $T4={T4:.1f}~K$")
plt.plot(V4_liq_prima,P4_liq_prima,color='cyan',linestyle='-',linewidth=2,label=f"Ajuste RK lí­quido T4")
plt.plot(V4_coex_prima,P4_coex_prima,color='cyan',linestyle='-',linewidth=2,label=f"Coexistencia T4")
plt.plot(V4_vapor_prima,P4_vapor_prima,color='cyan',linestyle='-',linewidth=2,label=f"Ajuste RK vapor T4")

plt.errorbar(V5,P5,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='magenta',ecolor='magenta',capsize=3,label='_nolegend_')
plt.scatter(V5,P5,color='magenta',label=f"Datos $T5={T5:.1f}~K$")
plt.plot(V5_liq_prima,P5_liq_prima,color='magenta',linestyle='-',linewidth=2,label=f"Ajuste RK lí­quido T5")
plt.plot(V5_coex_prima,P5_coex_prima,color='magenta',linestyle='-',linewidth=2,label=f"Coexistencia T5")
plt.plot(V5_vapor_prima,P5_vapor_prima,color='magenta',linestyle='-',linewidth=2,label=f"Ajuste RK vapor T5")

plt.errorbar(V6,P6,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='orange',ecolor='orange',capsize=3,label='_nolegend_')
plt.scatter(V6,P6,color='orange',label=f"Datos $T6={T6:.1f}~K$")
plt.plot(V6_total_prima,P6_total_prima,color='orange',linestyle='-',linewidth=2,label=f"Ajuste RK T6")

plt.title('Ajustes de Redlich-Kwong para las 6 isotermas', fontsize=25)
plt.xlabel(r'$V\;(\mathrm{mL})$', fontsize=25)
plt.ylabel(r'$P\;(\mathrm{bar})$', fontsize=25)
plt.legend(fontsize=12, ncol=2, loc='upper right')
plt.gca().set_facecolor('aliceblue')
plt.show()

## Cálculo de a, b, n y sus incertidumbres para cada isoterma ajustada:

R=83.14 # bar mL/ K mol

def abn(A,B,C,delt_A,delt_B,delt_C,R):
    n=C/R
    a=A/n**2
    b=B/n
    delt_n=delt_C/R
    delt_a=np.sqrt((delt_A/n**2)**2 + (2*delt_n*A/n**3)**2)
    delt_b=np.sqrt((delt_B/n)**2 + (delt_n*B/n**2)**2)
    return a,b,n,delt_a,delt_b,delt_n

# Con el ajuste de RK:

a1_prima,b1_prima,n1_prima,delt_a1_prima,delt_b1_prima,delt_n1_prima=abn(A1_prima,B1_prima,C1_prima,delt_A1_prima,delt_B1_prima,delt_C1_prima,R)
a2_prima,b2_prima,n2_prima,delt_a2_prima,delt_b2_prima,delt_n2_prima=abn(A2_prima,B2_prima,C2_prima,delt_A2_prima,delt_B2_prima,delt_C2_prima,R)
a3_prima,b3_prima,n3_prima,delt_a3_prima,delt_b3_prima,delt_n3_prima=abn(A3_prima,B3_prima,C3_prima,delt_A3_prima,delt_B3_prima,delt_C3_prima,R)
a4_prima,b4_prima,n4_prima,delt_a4_prima,delt_b4_prima,delt_n4_prima=abn(A4_prima,B4_prima,C4_prima,delt_A4_prima,delt_B4_prima,delt_C4_prima,R)
a5_prima,b5_prima,n5_prima,delt_a5_prima,delt_b5_prima,delt_n5_prima=abn(A5_prima,B5_prima,C5_prima,delt_A5_prima,delt_B5_prima,delt_C5_prima,R)
a6_prima,b6_prima,n6_prima,delt_a6_prima,delt_b6_prima,delt_n6_prima=abn(A6_prima,B6_prima,C6_prima,delt_A6_prima,delt_B6_prima,delt_C6_prima,R)

# Con el ajuste de VdW:

a1,b1,n1,delt_a1,delt_b1,delt_n1=abn(A1,B1,C1,delt_A1,delt_B1,delt_C1,R)
a2,b2,n2,delt_a2,delt_b2,delt_n2=abn(A2,B2,C2,delt_A2,delt_B2,delt_C2,R)
a3,b3,n3,delt_a3,delt_b3,delt_n3=abn(A3,B3,C3,delt_A3,delt_B3,delt_C3,R)
a4,b4,n4,delt_a4,delt_b4,delt_n4=abn(A4,B4,C4,delt_A4,delt_B4,delt_C4,R)
a5,b5,n5,delt_a5,delt_b5,delt_n5=abn(A5,B5,C5,delt_A5,delt_B5,delt_C5,R)
a6,b6,n6,delt_a6,delt_b6,delt_n6=abn(A6,B6,C6,delt_A6,delt_B6,delt_C6,R)

def media(x,n,delt_x):
    return np.sum(x)/n,np.sqrt(np.sum(delt_x**2)/n**2)

a_vdw,delt_a_vdw=media(np.array([a1,a2,a3,a4,a5,a6]),6,np.array([delt_a1,delt_a2,delt_a3,delt_a4,delt_a5,delt_a6]))
b_vdw,delt_b_vdw=media(np.array([b1,b2,b3,b4,b5,b6]),6,np.array([delt_b1,delt_b2,delt_b3,delt_b4,delt_b5,delt_b6]))
n_vdw,delt_n_vdw=media(np.array([n1,n2,n3,n4,n5,n6]),6,np.array([delt_n1,delt_n2,delt_n3,delt_n4,delt_n5,delt_n6]))

a_rk,delt_a_rk=media(np.array([a1_prima,a2_prima,a3_prima,a4_prima,a5_prima,a6_prima]),6,np.array([delt_a1_prima,delt_a2_prima,delt_a3_prima,delt_a4_prima,delt_a5_prima,delt_a6_prima]))
b_rk,delt_b_rk=media(np.array([b1_prima,b2_prima,b3_prima,b4_prima,b5_prima,b6_prima]),6,np.array([delt_b1_prima,delt_b2_prima,delt_b3_prima,delt_b4_prima,delt_b5_prima,delt_b6_prima]))
n_rk,delt_n_rk=media(np.array([n1_prima,n2_prima,n3_prima,n4_prima,n5_prima,n6_prima]),6,np.array([delt_n1_prima,delt_n2_prima,delt_n3_prima,delt_n4_prima,delt_n5_prima,delt_n6_prima]))

print(f"Valores medios para el ajuste de VdW:")
print(f"a = {a_vdw:.4f} ± {delt_a_vdw:.4f}")
print(f"b = {b_vdw:.4f} ± {delt_b_vdw:.4f}")
print(f"n = {n_vdw:.4f} ± {delt_n_vdw:.4f}")

print(f"Valores medios para el ajuste de RK:")
print(f"a = {a_rk:.4f} ± {delt_a_rk:.4f}")
print(f"b = {b_rk:.4f} ± {delt_b_rk:.4f}")
print(f"n = {n_rk:.4f} ± {delt_n_rk:.4f}")

## CÁLCULO DE LAS COORDENADAS CRÍTICAS:

def criticas_vdw(a,b,n,delt_a,delt_b,delt_n,R):
    P_c=a/(27*b**2)
    T_c=8*a/(27*R*b)
    V_c=3*n*b
    delt_P_c=np.sqrt((delt_a/(27*b**2))**2 + ((2*a*delt_b)/(27*b**3))**2)
    delt_T_c=np.sqrt(((8*delt_a)/(27*R*b))**2 + ((8*a*delt_b)/(27*R*b**2))**2)
    delt_V_c=np.sqrt((3*b*delt_n)**2 + (3*n*delt_b)**2)
    return P_c,T_c,V_c,delt_P_c,delt_T_c,delt_V_c

def criticas_rk(a,b,n,delt_a,delt_b,delt_n,R):
    k1=0.42748
    k2=0.08664
    Zc=1/3
    k=k1/k2
    T_c=(a/(k*R*b))**(2/3)
    P_c=k2*R*T_c/b
    Vc_molar=Zc*R*T_c/P_c
    V_c=n*Vc_molar
    dTc_da=(2/3)*T_c/a
    dTc_db=-(2/3)*T_c/b
    dPc_da=(k2*R/b)*dTc_da
    dPc_db=(k2*R/b)*dTc_db-(k2*R*T_c)/(b**2)
    dVc_da=Zc*R*((dTc_da*P_c - T_c*dPc_da)/(P_c**2))*n
    dVc_db=Zc*R*((dTc_db*P_c - T_c*dPc_db)/(P_c**2))*n
    dVc_dn=Vc_molar
    delt_T_c=np.sqrt((dTc_da*delt_a)**2 + (dTc_db*delt_b)**2)
    delt_P_c=np.sqrt((dPc_da*delt_a)**2 + (dPc_db*delt_b)**2)
    delt_V_c=np.sqrt((dVc_da*delt_a)**2 + (dVc_db*delt_b)**2 + (dVc_dn*delt_n)**2)
    return P_c,T_c,V_c,delt_P_c,delt_T_c,delt_V_c

P_c_vdw,T_c_vdw,V_c_vdw,delt_P_c_vdw,delt_T_c_vdw,delt_V_c_vdw=criticas_vdw(a_vdw,b_vdw,n_vdw,delt_a_vdw,delt_b_vdw,delt_n_vdw,R)
P_c_rk,T_c_rk,V_c_rk,delt_P_c_rk,delt_T_c_rk,delt_V_c_rk=criticas_rk(a_rk,b_rk,n_rk,delt_a_rk,delt_b_rk,delt_n_rk,R)

print(f"Coordenadas criticas VdW:")
print(f"P_c = {P_c_vdw:.4f} +/- {delt_P_c_vdw:.4f}")
print(f"T_c = {T_c_vdw:.4f} +/- {delt_T_c_vdw:.4f}")
print(f"V_c = {V_c_vdw:.4f} +/- {delt_V_c_vdw:.4f}")

print(f"Coordenadas criticas RK:")
print(f"P_c = {P_c_rk:.4f} +/- {delt_P_c_rk:.4f}")
print(f"T_c = {T_c_rk:.4f} +/- {delt_T_c_rk:.4f}")
print(f"V_c = {V_c_rk:.4f} +/- {delt_V_c_rk:.4f}")



