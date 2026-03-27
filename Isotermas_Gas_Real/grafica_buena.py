## RESULTADOS PARA LAS ISOTERMAS DE UN GAS REAL:

#%%

import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
plt.close('all')
plt.rcParams['lines.markersize'] = 5

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

isoterma_1_ajuste=datos_original['Isoterma 1']
isoterma_2_ajuste=datos_original['Isoterma 2']
isoterma_3_ajuste=datos_original['Isoterma 3']
isoterma_4_ajuste=datos_original['Isoterma 4']
isoterma_5_ajuste=datos_original['Isoterma 5']
isoterma_6_ajuste=datos_original['Isoterma 6']

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

def p0_desde_datos(V,P,T):
    V0 = np.median(V)
    P0 = np.median(P)
    C0 = (P0*V0)/T
    A0 = P0*V0**2
    B0 = 0.05*V0
    return [A0, B0, C0]

# Presiones de coexistencia (usar para pesos/intervalos experimentales)
P1_coex=21.5
P2_coex=24
P3_coex=26.5
P4_coex=29.5
P5_coex=32.5
P6_coex=0
delta_tramos=0.5

def sigma_experimental(P, Pcoex, delta, sigma_liq, sigma_vap, sigma_coex):
    mask_coex = (P >= (Pcoex - delta)) & (P <= (Pcoex + delta))
    mask_liq = P > (Pcoex + delta)
    mask_vap = P < (Pcoex - delta)
    sigma = np.where(mask_coex, sigma_coex, np.where(mask_vap, sigma_vap, sigma_liq))
    return sigma

def tramos_experimentales(V, P, Pcoex, delta, npts=100, tail_points=0, tail_frac=0.03, enforce_monotone_liq=False):
    mask_liq = P > (Pcoex + delta)
    mask_vap = P < (Pcoex - delta)
    mask_coex = ~(mask_liq | mask_vap)
    if enforce_monotone_liq and np.any(mask_liq):
        idx_liq = np.where(mask_liq)[0]
        order = np.argsort(V[idx_liq])
        p_liq = P[idx_liq][order]
        i_min = int(np.argmin(p_liq))
        keep_idx = idx_liq[order][: i_min + 1]
        new_mask = np.zeros_like(mask_liq, dtype=bool)
        new_mask[keep_idx] = True
        mask_liq = new_mask
    V_liq = None
    V_vap = None
    V_coex = None
    if np.any(mask_coex):
        v_coex_min = np.min(V[mask_coex])
        v_coex_max = np.max(V[mask_coex])
        V_coex = np.linspace(v_coex_min, v_coex_max, npts)
    if np.any(mask_liq):
        v_liq_min = np.min(V[mask_liq])
        v_liq_max = np.max(V[mask_liq])
        if V_coex is not None:
            v_liq_max = v_coex_min
        if v_liq_max >= v_liq_min:
            V_liq = np.linspace(v_liq_min, v_liq_max, npts)
    if np.any(mask_vap):
        v_vap_min = np.min(V[mask_vap])
        v_vap_max = np.max(V[mask_vap])
        if V_coex is not None:
            v_vap_min = v_coex_max
        if v_vap_max >= v_vap_min:
            V_vap = np.linspace(v_vap_min, v_vap_max, npts)
    if V_vap is not None and tail_points > 0:
        V_all_max = np.max(V)
        if V_all_max > V_vap[-1]:
            span = max(V_all_max - V_vap[-1], tail_frac * V_all_max)
            V_tail = np.linspace(V_vap[-1], min(V_vap[-1] + span, V_all_max), tail_points, endpoint=True)
            V_vap = np.unique(np.concatenate([V_vap, V_tail]))
    return V_liq, V_coex, V_vap

def cola_visual(V_coex, V_pts, P_pts, Pcoex, delt_P, npts=60):
    if V_coex is None or V_pts is None or P_pts is None or len(V_pts) == 0:
        return None, None
    v0 = V_coex[-1]
    v1 = np.max(V_pts)
    if v1 <= v0:
        return None, None
    p1 = P_pts[np.argmax(V_pts)]
    shift = 1.2 * delt_P * (-1 if p1 < Pcoex else 1)
    p1s = p1 + shift
    V_tail = np.linspace(v0, v1, npts)
    P_tail = np.linspace(Pcoex, p1s, npts)
    return V_tail, P_tail

def ajusta_union_coex(Pcoex, V_liq, P_liq, V_coex, V_vap, P_vap, V_exp, P_exp, npts=100, nblend=12, p_margin=1.0):
    v_coex_min = np.min(V_coex) if V_coex is not None else None
    v_coex_max = np.max(V_coex) if V_coex is not None else None
    if V_liq is not None and P_liq is not None:
        i_min = int(np.argmin(P_liq))
        min_val = P_liq[i_min]
        if min_val < Pcoex:
            # Caso: minimo por debajo de Pcoex -> corta al cruzar Pcoex y une en seco
            V_liq = V_liq[: i_min + 1]
            P_liq = P_liq[: i_min + 1]
            idx = np.where(P_liq <= Pcoex)[0]
            if len(idx) > 0:
                j = idx[0]
                if j == 0:
                    V_cross = V_liq[0]
                    P_cross = P_liq[0]
                else:
                    V0, V1 = V_liq[j - 1], V_liq[j]
                    P0, P1 = P_liq[j - 1], P_liq[j]
                    V_cross = V0 + (P0 - Pcoex) * (V1 - V0) / (P0 - P1)
                    P_cross = Pcoex
                V_liq = np.concatenate([V_liq[:j], [V_cross]])
                P_liq = np.concatenate([P_liq[:j], [P_cross]])
            v_coex_min = V_liq[-1] if v_coex_min is None else min(v_coex_min, V_liq[-1])
        else:
            # Caso: minimo por encima de Pcoex -> corta antes del minimo e interpola suave hasta coex
            i_anchor = max(0, i_min - nblend)
            V_liq = V_liq[: i_anchor + 1]
            P_liq = P_liq[: i_anchor + 1]
            v_l_end = V_liq[-1]
            if v_coex_min is None:
                v_coex_min = v_l_end
            # Alarga el tramo horizontal hacia la izquierda hasta que corte la interpolacion en P = Pcoex
            v_sorted = np.sort(V_exp)
            dv = np.median(np.diff(v_sorted)) if len(v_sorted) > 1 else 0.02
            if v_coex_min is None:
                v_coex_min = v_l_end + max(dv, 0.02)
            v_coex_min = max(v_coex_min, v_l_end + max(dv, 0.02))
            if v_coex_min > v_l_end:
                nblend_use = nblend
                V_ext = np.linspace(v_l_end, v_coex_min, nblend_use, endpoint=True)
                # Interpolacion suave con pendiente coherente con el tramo liquido y pendiente 0 al llegar a coexistencia
                if len(V_liq) >= 2:
                    P0 = P_liq[-1]
                    V0 = V_liq[-1]
                    P1 = Pcoex
                    V1 = v_coex_min
                    dV = V1 - V0
                    m0 = (P_liq[-1] - P_liq[-2]) / (V_liq[-1] - V_liq[-2]) if (V_liq[-1] - V_liq[-2]) != 0 else 0.0
                    m1 = 0.0
                    t = np.linspace(0.0, 1.0, nblend_use)
                    h00 = (1 + 2*t) * (1 - t) ** 2
                    h10 = t * (1 - t) ** 2
                    h01 = t ** 2 * (3 - 2*t)
                    h11 = t ** 2 * (t - 1)
                    P_ext = h00 * P0 + h10 * m0 * dV + h01 * P1 + h11 * m1 * dV
                else:
                    t = np.linspace(0.0, 1.0, nblend_use)
                    P_start = P_liq[-1]
                    P_ext = Pcoex + (P_start - Pcoex) * (1 - t) ** 2
                idx_cut = np.where(P_ext <= Pcoex)[0]
                k = idx_cut[0] if len(idx_cut) > 0 else (len(P_ext) - 1)
                V_liq = np.concatenate([V_liq, V_ext[1:k+1]])
                P_liq = np.concatenate([P_liq, P_ext[1:k+1]])
                v_coex_min = V_liq[-1]
            v_coex_min = V_liq[-1] if v_coex_min is None else min(v_coex_min, V_liq[-1])
    if V_vap is not None and P_vap is not None:
        i_join = int(np.argmin(np.abs(P_vap - Pcoex)))
        V_vap = V_vap[i_join:]
        P_vap = P_vap[i_join:]
        v_coex_max = V_vap[0] if v_coex_max is None else max(v_coex_max, V_vap[0])
    P_coex_line = None
    if v_coex_min is not None and v_coex_max is not None:
        if v_coex_max < v_coex_min:
            v_coex_min, v_coex_max = v_coex_max, v_coex_min
        V_coex = np.linspace(v_coex_min, v_coex_max, npts)
        P_coex_line = np.full_like(V_coex, Pcoex)
    return V_liq, P_liq, V_coex, P_coex_line, V_vap, P_vap


# Ajuste para la isoterma 1 (T1):

T=T1 # Ajuste
p0 = p0_desde_datos(V1_ajuste,P1_ajuste,T1)
params,cov=curve_fit(isoterma_Van_der_Waals,V1_ajuste,P1_ajuste,p0=p0,bounds = ([0, 0, 0], [np.inf, 0.9*np.min(V1_ajuste), np.inf]),sigma=sigma_experimental(P1_ajuste,P1_coex,delta_tramos,0.5*delt_P,0.3*delt_P,3*delt_P),absolute_sigma=True,maxfev=100000)

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
raices_reales_1.sort()
print(f"Volúmenes reales a presión {P1_coex} bar: {raices_reales_1}")

V1_liq, V1_coex, V1_vapor = tramos_experimentales(V1, P1, P1_coex, delta_tramos, tail_points=8, tail_frac=0.04, enforce_monotone_liq=True)
P1_liq = isoterma_Van_der_Waals(V1_liq,A1,B1,C1) if V1_liq is not None else None
P1_coex_line = np.full_like(V1_coex,P1_coex) if V1_coex is not None else None
P1_vapor = isoterma_Van_der_Waals(V1_vapor,A1,B1,C1) if V1_vapor is not None else None
mask_vap_1 = P1 < (P1_coex - delta_tramos)
V1_vap_pts = V1[mask_vap_1]
P1_vap_pts = P1[mask_vap_1]
usar_cola_visual_t1 = len(V1_vap_pts) < 3
if usar_cola_visual_t1:
    V1_vapor, P1_vapor = cola_visual(V1_coex, V1_vap_pts, P1_vap_pts, P1_coex, delt_P, npts=60)
V1_liq, P1_liq, V1_coex, P1_coex_line, V1_vapor, P1_vapor = ajusta_union_coex(P1_coex, V1_liq, P1_liq, V1_coex, V1_vapor, P1_vapor, V1, P1)

# Ajuste para la isoterma 2 (T2):

T=T2
p0 = p0_desde_datos(V2_ajuste,P2_ajuste,T2)
params,cov=curve_fit(isoterma_Van_der_Waals,V2_ajuste,P2_ajuste,p0=p0,bounds=([0, 0, 0], [np.inf, 0.9*np.min(V2_ajuste), np.inf]),sigma=sigma_experimental(P2_ajuste,P2_coex,delta_tramos,0.5*delt_P,0.3*delt_P,3*delt_P),absolute_sigma=True,maxfev=100000)

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
raices_reales_2.sort()
print(f"Volúmenes reales a presión {P2_coex} bar: {raices_reales_2}")

V2_liq, V2_coex, V2_vapor = tramos_experimentales(V2, P2, P2_coex, delta_tramos, enforce_monotone_liq=True)
P2_liq = isoterma_Van_der_Waals(V2_liq,A2,B2,C2) if V2_liq is not None else None
P2_coex_line = np.full_like(V2_coex,P2_coex) if V2_coex is not None else None
P2_vapor = isoterma_Van_der_Waals(V2_vapor,A2,B2,C2) if V2_vapor is not None else None
V2_liq, P2_liq, V2_coex, P2_coex_line, V2_vapor, P2_vapor = ajusta_union_coex(P2_coex, V2_liq, P2_liq, V2_coex, V2_vapor, P2_vapor, V2, P2)

# Ajuste para la isoterma 3 (T3):

T=T3
p0 = p0_desde_datos(V3_ajuste,P3_ajuste,T3)
params,cov=curve_fit(isoterma_Van_der_Waals,V3_ajuste,P3_ajuste,p0=p0,bounds=([0, 0, 0], [np.inf, 0.9*np.min(V3_ajuste), np.inf]),sigma=sigma_experimental(P3_ajuste,P3_coex,delta_tramos,0.5*delt_P,0.3*delt_P,3*delt_P),absolute_sigma=True,maxfev=100000)

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
raices_reales_3.sort()
print(f"Volúmenes reales a presión {P3_coex} bar: {raices_reales_3}")

V3_liq, V3_coex, V3_vapor = tramos_experimentales(V3, P3, P3_coex, delta_tramos, enforce_monotone_liq=True)
P3_liq = isoterma_Van_der_Waals(V3_liq,A3,B3,C3) if V3_liq is not None else None
P3_coex_line = np.full_like(V3_coex,P3_coex) if V3_coex is not None else None
P3_vapor = isoterma_Van_der_Waals(V3_vapor,A3,B3,C3) if V3_vapor is not None else None
V3_liq, P3_liq, V3_coex, P3_coex_line, V3_vapor, P3_vapor = ajusta_union_coex(P3_coex, V3_liq, P3_liq, V3_coex, V3_vapor, P3_vapor, V3, P3)



# Ajuste para la isoterma 4 (T4):

T=T4
p0 = p0_desde_datos(V4_ajuste,P4_ajuste,T4)
params,cov=curve_fit(isoterma_Van_der_Waals,V4_ajuste,P4_ajuste,p0=p0,bounds=([0, 0, 0], [np.inf, 0.9*np.min(V4_ajuste), np.inf]),sigma=sigma_experimental(P4_ajuste,P4_coex,delta_tramos,0.5*delt_P,0.3*delt_P,3*delt_P),absolute_sigma=True,maxfev=100000)

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
raices_reales_4.sort()
print(f"Volúmenes reales a presión {P4_coex} bar: {raices_reales_4}")

V4_liq, V4_coex, V4_vapor = tramos_experimentales(V4, P4, P4_coex, delta_tramos, enforce_monotone_liq=True)
P4_liq = isoterma_Van_der_Waals(V4_liq,A4,B4,C4) if V4_liq is not None else None
P4_coex_line = np.full_like(V4_coex,P4_coex) if V4_coex is not None else None
P4_vapor = isoterma_Van_der_Waals(V4_vapor,A4,B4,C4) if V4_vapor is not None else None
V4_liq, P4_liq, V4_coex, P4_coex_line, V4_vapor, P4_vapor = ajusta_union_coex(P4_coex, V4_liq, P4_liq, V4_coex, V4_vapor, P4_vapor, V4, P4)



# Ajuste para la isoterma 5 (T5):

T=T5
p0 = p0_desde_datos(V5_ajuste,P5_ajuste,T5)
params,cov=curve_fit(isoterma_Van_der_Waals,V5_ajuste,P5_ajuste,p0=p0,bounds=([0, 0, 0], [np.inf, 0.9*np.min(V5_ajuste), np.inf]),sigma=sigma_experimental(P5_ajuste,P5_coex,delta_tramos,0.5*delt_P,0.3*delt_P,3*delt_P),absolute_sigma=True,maxfev=100000)

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
P6_coex=0
coefs_5=np.array([P5_coex,-(P5_coex*B5+C5*T5),A5,-A5*B5])
raices_5=np.roots(coefs_5)
raices_reales_5=raices_5[np.isreal(raices_5)].real
raices_reales_5.sort()
print(f"Volúmenes reales a presión {P5_coex} bar: {raices_reales_5}")

V5_liq, V5_coex, V5_vapor = tramos_experimentales(V5, P5, P5_coex, delta_tramos, enforce_monotone_liq=True)
P5_liq = isoterma_Van_der_Waals(V5_liq,A5,B5,C5) if V5_liq is not None else None
P5_coex_line = np.full_like(V5_coex,P5_coex) if V5_coex is not None else None
P5_vapor = isoterma_Van_der_Waals(V5_vapor,A5,B5,C5) if V5_vapor is not None else None
V5_liq, P5_liq, V5_coex, P5_coex_line, V5_vapor, P5_vapor = ajusta_union_coex(P5_coex, V5_liq, P5_liq, V5_coex, V5_vapor, P5_vapor, V5, P5)



# Ajuste para la isoterma 6 (T6):

T=T6
p0 = p0_desde_datos(V6_ajuste,P6_ajuste,T6)
params,cov=curve_fit(isoterma_Van_der_Waals,V6_ajuste,P6_ajuste,p0=p0,bounds=([0, 0, 0], [np.inf, 0.9*np.min(V6_ajuste), np.inf]),sigma=sigma_experimental(P6_ajuste,P6_coex,delta_tramos,0.5*delt_P,0.3*delt_P,3*delt_P),absolute_sigma=True,maxfev=100000)

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
if V1_liq is not None:
    plt.plot(V1_liq,P1_liq,color='blue',linestyle='-',linewidth=2,label=f"Ajuste VdW líquido T1")
if V1_coex is not None:
    plt.plot(V1_coex,P1_coex_line,color='blue',linestyle='--',linewidth=2,label=f"Coexistencia T1")
if V1_vapor is not None:
    plt.plot(V1_vapor,P1_vapor,color='blue',linestyle='-',linewidth=2,label=f"Ajuste VdW vapor T1")

plt.errorbar(V2,P2,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='green',ecolor='green',capsize=3,label='_nolegend_')
plt.scatter(V2,P2,color='green',label=f"Datos $T2={T2:.1f}~K$")
if V2_liq is not None:
    plt.plot(V2_liq,P2_liq,color='green',linestyle='-',linewidth=2,label=f"Ajuste VdW líquido T2")
if V2_coex is not None:
    plt.plot(V2_coex,P2_coex_line,color='green',linestyle='--',linewidth=2,label=f"Coexistencia T2")
if V2_vapor is not None:
    plt.plot(V2_vapor,P2_vapor,color='green',linestyle='-',linewidth=2,label=f"Ajuste VdW vapor T2")

plt.errorbar(V3,P3,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='red',ecolor='red',capsize=3,label='_nolegend_')
plt.scatter(V3,P3,color='red',label=f"Datos $T3={T3:.1f}~K$")
if V3_liq is not None:
    plt.plot(V3_liq,P3_liq,color='red',linestyle='-',linewidth=2,label=f"Ajuste VdW líquido T3")
if V3_coex is not None:
    plt.plot(V3_coex,P3_coex_line,color='red',linestyle='--',linewidth=2,label=f"Coexistencia T3")
if V3_vapor is not None:
    plt.plot(V3_vapor,P3_vapor,color='red',linestyle='-',linewidth=2,label=f"Ajuste VdW vapor T3")

plt.errorbar(V4,P4,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='cyan',ecolor='cyan',capsize=3,label='_nolegend_')
plt.scatter(V4,P4,color='cyan',label=f"Datos $T4={T4:.1f}~K$")
if V4_liq is not None:
    plt.plot(V4_liq,P4_liq,color='cyan',linestyle='-',linewidth=2,label=f"Ajuste VdW líquido T4")
if V4_coex is not None:
    plt.plot(V4_coex,P4_coex_line,color='cyan',linestyle='--',linewidth=2,label=f"Coexistencia T4")
if V4_vapor is not None:
    plt.plot(V4_vapor,P4_vapor,color='cyan',linestyle='-',linewidth=2,label=f"Ajuste VdW vapor T4")

plt.errorbar(V5,P5,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='magenta',ecolor='magenta',capsize=3,label='_nolegend_')
plt.scatter(V5,P5,color='magenta',label=f"Datos $T5={T5:.1f}~K$")
if V5_liq is not None:
    plt.plot(V5_liq,P5_liq,color='magenta',linestyle='-',linewidth=2,label=f"Ajuste VdW líquido T5")
if V5_coex is not None:
    plt.plot(V5_coex,P5_coex_line,color='magenta',linestyle='--',linewidth=2,label=f"Coexistencia T5")
if V5_vapor is not None:
    plt.plot(V5_vapor,P5_vapor,color='magenta',linestyle='-',linewidth=2,label=f"Ajuste VdW vapor T5")

plt.errorbar(V6,P6,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='orange',ecolor='orange',capsize=3,label='_nolegend_')
plt.scatter(V6,P6,color='orange',label=f"Datos $T6={T6:.1f}~K$")
plt.plot(V6_total,P6_total,color='orange',linestyle='-',linewidth=2,label=f"Ajuste VdW T6")

plt.xlabel(r'$V\;(\mathrm{mL})$', fontsize=25)
plt.ylabel(r'$P\;(\mathrm{bar})$', fontsize=25)
plt.legend(fontsize=12, ncol=2, loc='upper right')
ax_vdw = plt.gca()
plt.gca().set_facecolor('aliceblue')
plt.show(block=False)



## AJUSTE REDLICH-KWONG:

# Ajuste para la isoterma 1 (T1):

T=T1
params_prima,cov_prima=curve_fit(Redlich_Kwong,V1_ajuste,P1_ajuste,p0=p0,sigma=sigma_experimental(P1_ajuste,P1_coex,delta_tramos,0.5*delt_P,0.3*delt_P,3*delt_P),absolute_sigma=True,maxfev=100000)

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
raices_reales_1_prima.sort()
print(f"Volúmenes reales RK a presión {P1_coex} bar: {raices_reales_1_prima}")

V1_liq_prima, V1_coex_prima, V1_vapor_prima = tramos_experimentales(V1, P1, P1_coex, delta_tramos, tail_points=8, tail_frac=0.04, enforce_monotone_liq=True)
P1_liq_prima = Redlich_Kwong(V1_liq_prima,A1_prima,B1_prima,C1_prima) if V1_liq_prima is not None else None
P1_coex_prima_line = np.full_like(V1_coex_prima,P1_coex) if V1_coex_prima is not None else None
P1_vapor_prima = Redlich_Kwong(V1_vapor_prima,A1_prima,B1_prima,C1_prima) if V1_vapor_prima is not None else None
if usar_cola_visual_t1:
    V1_vapor_prima, P1_vapor_prima = cola_visual(V1_coex_prima, V1_vap_pts, P1_vap_pts, P1_coex, delt_P, npts=60)
V1_liq_prima, P1_liq_prima, V1_coex_prima, P1_coex_prima_line, V1_vapor_prima, P1_vapor_prima = ajusta_union_coex(P1_coex, V1_liq_prima, P1_liq_prima, V1_coex_prima, V1_vapor_prima, P1_vapor_prima, V1, P1)

# Ajuste para la isoterma 2 (T2):

T=T2
p0 = p0_desde_datos(V2_ajuste,P2_ajuste,T2)
params_prima,cov_prima=curve_fit(Redlich_Kwong,V2_ajuste,P2_ajuste,p0=p0,sigma=sigma_experimental(P2_ajuste,P2_coex,delta_tramos,0.5*delt_P,0.3*delt_P,3*delt_P),absolute_sigma=True,maxfev=100000)

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
raices_reales_2_prima.sort()
print(f"Volúmenes reales RK a presión {P2_coex} bar: {raices_reales_2_prima}")

V2_liq_prima, V2_coex_prima, V2_vapor_prima = tramos_experimentales(V2, P2, P2_coex, delta_tramos, enforce_monotone_liq=True)
P2_liq_prima = Redlich_Kwong(V2_liq_prima,A2_prima,B2_prima,C2_prima) if V2_liq_prima is not None else None
P2_coex_prima_line = np.full_like(V2_coex_prima,P2_coex) if V2_coex_prima is not None else None
P2_vapor_prima = Redlich_Kwong(V2_vapor_prima,A2_prima,B2_prima,C2_prima) if V2_vapor_prima is not None else None
V2_liq_prima, P2_liq_prima, V2_coex_prima, P2_coex_prima_line, V2_vapor_prima, P2_vapor_prima = ajusta_union_coex(P2_coex, V2_liq_prima, P2_liq_prima, V2_coex_prima, V2_vapor_prima, P2_vapor_prima, V2, P2)

# Ajuste para la isoterma 3 (T3):

T=T3
p0 = p0_desde_datos(V3_ajuste,P3_ajuste,T3)
params_prima,cov_prima=curve_fit(Redlich_Kwong,V3_ajuste,P3_ajuste,p0=p0,sigma=sigma_experimental(P3_ajuste,P3_coex,delta_tramos,0.5*delt_P,0.3*delt_P,3*delt_P),absolute_sigma=True,maxfev=100000)

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
raices_reales_3_prima.sort()
print(f"Volúmenes reales RK a presión {P3_coex} bar: {raices_reales_3_prima}")

V3_liq_prima, V3_coex_prima, V3_vapor_prima = tramos_experimentales(V3, P3, P3_coex, delta_tramos, enforce_monotone_liq=True)
P3_liq_prima = Redlich_Kwong(V3_liq_prima,A3_prima,B3_prima,C3_prima) if V3_liq_prima is not None else None
P3_coex_prima_line = np.full_like(V3_coex_prima,P3_coex) if V3_coex_prima is not None else None
P3_vapor_prima = Redlich_Kwong(V3_vapor_prima,A3_prima,B3_prima,C3_prima) if V3_vapor_prima is not None else None
V3_liq_prima, P3_liq_prima, V3_coex_prima, P3_coex_prima_line, V3_vapor_prima, P3_vapor_prima = ajusta_union_coex(P3_coex, V3_liq_prima, P3_liq_prima, V3_coex_prima, V3_vapor_prima, P3_vapor_prima, V3, P3)
if V3_coex_prima is not None and P3_liq_prima is not None and len(P3_liq_prima) > 0:
    P3_coex_prima_line = np.full_like(V3_coex_prima, P3_liq_prima[-1])

# Ajuste para la isoterma 4 (T4):

T=T4
p0 = p0_desde_datos(V4_ajuste,P4_ajuste,T4)
params_prima,cov_prima=curve_fit(Redlich_Kwong,V4_ajuste,P4_ajuste,p0=p0,sigma=sigma_experimental(P4_ajuste,P4_coex,delta_tramos,0.5*delt_P,0.3*delt_P,3*delt_P),absolute_sigma=True,maxfev=100000)

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
raices_reales_4_prima.sort()
print(f"Volúmenes reales RK a presión {P4_coex} bar: {raices_reales_4_prima}")

V4_liq_prima, V4_coex_prima, V4_vapor_prima = tramos_experimentales(V4, P4, P4_coex, delta_tramos, enforce_monotone_liq=True)
P4_liq_prima = Redlich_Kwong(V4_liq_prima,A4_prima,B4_prima,C4_prima) if V4_liq_prima is not None else None
P4_coex_prima_line = np.full_like(V4_coex_prima,P4_coex) if V4_coex_prima is not None else None
P4_vapor_prima = Redlich_Kwong(V4_vapor_prima,A4_prima,B4_prima,C4_prima) if V4_vapor_prima is not None else None
V4_liq_prima, P4_liq_prima, V4_coex_prima, P4_coex_prima_line, V4_vapor_prima, P4_vapor_prima = ajusta_union_coex(P4_coex, V4_liq_prima, P4_liq_prima, V4_coex_prima, V4_vapor_prima, P4_vapor_prima, V4, P4)

# Ajuste para la isoterma 5 (T5):

T=T5
p0 = p0_desde_datos(V5_ajuste,P5_ajuste,T5)
params_prima,cov_prima=curve_fit(Redlich_Kwong,V5_ajuste,P5_ajuste,p0=p0,sigma=sigma_experimental(P5_ajuste,P5_coex,delta_tramos,0.5*delt_P,0.3*delt_P,3*delt_P),absolute_sigma=True,maxfev=100000)

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
P6_coex=0
coefs_5_prima=np.array([P5_coex,-(C5_prima*T5),(-P5_coex*(B5_prima**2)-C5_prima*T5*B5_prima + A5_prima/np.sqrt(T5)),(-A5_prima*B5_prima/np.sqrt(T5))])
raices_5_prima=np.roots(coefs_5_prima)
raices_reales_5_prima=raices_5_prima[np.isreal(raices_5_prima)].real
raices_reales_5_prima.sort()
print(f"Volúmenes reales RK a presión {P5_coex} bar: {raices_reales_5_prima}")

V5_liq_prima, V5_coex_prima, V5_vapor_prima = tramos_experimentales(V5, P5, P5_coex, delta_tramos, enforce_monotone_liq=True)
P5_liq_prima = Redlich_Kwong(V5_liq_prima,A5_prima,B5_prima,C5_prima) if V5_liq_prima is not None else None
P5_coex_prima_line = np.full_like(V5_coex_prima,P5_coex) if V5_coex_prima is not None else None
P5_vapor_prima = Redlich_Kwong(V5_vapor_prima,A5_prima,B5_prima,C5_prima) if V5_vapor_prima is not None else None
V5_liq_prima, P5_liq_prima, V5_coex_prima, P5_coex_prima_line, V5_vapor_prima, P5_vapor_prima = ajusta_union_coex(P5_coex, V5_liq_prima, P5_liq_prima, V5_coex_prima, V5_vapor_prima, P5_vapor_prima, V5, P5)

# Ajuste para la isoterma 6 (T6):

T=T6
p0 = p0_desde_datos(V6_ajuste,P6_ajuste,T6)
params_prima,cov_prima=curve_fit(Redlich_Kwong,V6_ajuste,P6_ajuste,p0=p0,sigma=sigma_experimental(P6_ajuste,P6_coex,delta_tramos,0.5*delt_P,0.3*delt_P,3*delt_P),absolute_sigma=True,maxfev=100000)

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
if V1_liq_prima is not None:
    plt.plot(V1_liq_prima,P1_liq_prima,color='blue',linestyle='-',linewidth=2,label=f"Ajuste RK líquido T1")
if V1_coex_prima is not None:
    plt.plot(V1_coex_prima,P1_coex_prima_line,color='blue',linestyle='--',linewidth=2,label=f"Coexistencia T1")
if V1_vapor_prima is not None:
    plt.plot(V1_vapor_prima,P1_vapor_prima,color='blue',linestyle='-',linewidth=2,label=f"Ajuste RK vapor T1")

plt.errorbar(V2,P2,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='green',ecolor='green',capsize=3,label='_nolegend_')
plt.scatter(V2,P2,color='green',label=f"Datos $T2={T2:.1f}~K$")
if V2_liq_prima is not None:
    plt.plot(V2_liq_prima,P2_liq_prima,color='green',linestyle='-',linewidth=2,label=f"Ajuste RK líquido T2")
if V2_coex_prima is not None:
    plt.plot(V2_coex_prima,P2_coex_prima_line,color='green',linestyle='--',linewidth=2,label=f"Coexistencia T2")
if V2_vapor_prima is not None:
    plt.plot(V2_vapor_prima,P2_vapor_prima,color='green',linestyle='-',linewidth=2,label=f"Ajuste RK vapor T2")

plt.errorbar(V3,P3,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='red',ecolor='red',capsize=3,label='_nolegend_')
plt.scatter(V3,P3,color='red',label=f"Datos $T3={T3:.1f}~K$")
if V3_liq_prima is not None:
    plt.plot(V3_liq_prima,P3_liq_prima,color='red',linestyle='-',linewidth=2,label=f"Ajuste RK líquido T3")
if V3_coex_prima is not None:
    plt.plot(V3_coex_prima,P3_coex_prima_line,color='red',linestyle='--',linewidth=2,label=f"Coexistencia T3")
if V3_vapor_prima is not None:
    plt.plot(V3_vapor_prima,P3_vapor_prima,color='red',linestyle='-',linewidth=2,label=f"Ajuste RK vapor T3")

plt.errorbar(V4,P4,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='cyan',ecolor='cyan',capsize=3,label='_nolegend_')
plt.scatter(V4,P4,color='cyan',label=f"Datos $T4={T4:.1f}~K$")
if V4_liq_prima is not None:
    plt.plot(V4_liq_prima,P4_liq_prima,color='cyan',linestyle='-',linewidth=2,label=f"Ajuste RK líquido T4")
if V4_coex_prima is not None:
    plt.plot(V4_coex_prima,P4_coex_prima_line,color='cyan',linestyle='--',linewidth=2,label=f"Coexistencia T4")
if V4_vapor_prima is not None:
    plt.plot(V4_vapor_prima,P4_vapor_prima,color='cyan',linestyle='-',linewidth=2,label=f"Ajuste RK vapor T4")

plt.errorbar(V5,P5,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='magenta',ecolor='magenta',capsize=3,label='_nolegend_')
plt.scatter(V5,P5,color='magenta',label=f"Datos $T5={T5:.1f}~K$")
if V5_liq_prima is not None:
    plt.plot(V5_liq_prima,P5_liq_prima,color='magenta',linestyle='-',linewidth=2,label=f"Ajuste RK líquido T5")
if V5_coex_prima is not None:
    plt.plot(V5_coex_prima,P5_coex_prima_line,color='magenta',linestyle='--',linewidth=2,label=f"Coexistencia T5")
if V5_vapor_prima is not None:
    plt.plot(V5_vapor_prima,P5_vapor_prima,color='magenta',linestyle='-',linewidth=2,label=f"Ajuste RK vapor T5")

plt.errorbar(V6,P6,yerr=delt_P,xerr=delt_V,marker='o',linestyle='',color='orange',ecolor='orange',capsize=3,label='_nolegend_')
plt.scatter(V6,P6,color='orange',label=f"Datos $T6={T6:.1f}~K$")
plt.plot(V6_total_prima,P6_total_prima,color='orange',linestyle='-',linewidth=2,label=f"Ajuste RK T6")

plt.xlabel(r'$V\;(\mathrm{mL})$', fontsize=25)
plt.ylabel(r'$P\;(\mathrm{bar})$', fontsize=25)
plt.legend(fontsize=12, ncol=2, loc='upper right')
ax_rk = plt.gca()
plt.gca().set_facecolor('aliceblue')
plt.show(block=False)

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
    # Expresiones indicadas
    P_c=((k2**5 * a**2 * R)/(k1**2 * b**5))**(1/3)
    T_c=((k2*a)/(k1*R*b))**(2/3)
    V_c=(n*b)/(3*k2)
    # Propagacion de errores
    dPc_da=(2/3)*P_c/a
    dPc_db=-(5/3)*P_c/b
    dTc_da=(2/3)*T_c/a
    dTc_db=-(2/3)*T_c/b
    dVc_db=(1/3)*n/k2
    dVc_dn=(1/3)*b/k2
    delt_P_c=np.sqrt((dPc_da*delt_a)**2 + (dPc_db*delt_b)**2)
    delt_T_c=np.sqrt((dTc_da*delt_a)**2 + (dTc_db*delt_b)**2)
    delt_V_c=np.sqrt((dVc_db*delt_b)**2 + (dVc_dn*delt_n)**2)
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

def compresibilidad(P_c, V_c, T_c, n, delt_P_c, delt_V_c, delt_T_c, delt_n, R):
    Z_c = (P_c * V_c) / (n * R * T_c)
    dZ_dP = Z_c / P_c
    dZ_dV = Z_c / V_c
    dZ_dT = -Z_c / T_c
    dZ_dn = -Z_c / n
    delt_Z_c = np.sqrt((dZ_dP * delt_P_c)**2 + (dZ_dV * delt_V_c)**2 +
                       (dZ_dT * delt_T_c)**2 + (dZ_dn * delt_n)**2)
    return Z_c, delt_Z_c

Z_c_vdw, delt_Z_c_vdw = compresibilidad(P_c_vdw, V_c_vdw, T_c_vdw, n_vdw,
                                        delt_P_c_vdw, delt_V_c_vdw, delt_T_c_vdw, delt_n_vdw, R)
Z_c_rk, delt_Z_c_rk = compresibilidad(P_c_rk, V_c_rk, T_c_rk, n_rk,
                                      delt_P_c_rk, delt_V_c_rk, delt_T_c_rk, delt_n_rk, R)

print(f"VdW: Z_c = {Z_c_vdw:.4f} +/- {delt_Z_c_vdw:.4f} (adimensional)")
print(f"RK: Z_c = {Z_c_rk:.4f} +/- {delt_Z_c_rk:.4f} (adimensional)")

# Puntos criticos en las graficas de ajuste
if 'ax_vdw' in globals():
    P_c_vdw_plot = P_c_vdw
    V_c_vdw_plot = V_c_vdw
    ax_vdw.errorbar(V_c_vdw_plot, P_c_vdw_plot, xerr=delt_V_c_vdw, yerr=delt_P_c_vdw,
                    color='black', ecolor='black', capsize=2,
                    linestyle='none', zorder=6)
    ax_vdw.scatter(V_c_vdw_plot, P_c_vdw_plot, color='black', marker='o', s=40, label='Punto crítico', zorder=7)
    ax_vdw.legend(fontsize=12, ncol=2, loc='upper right')
    ax_vdw.figure.canvas.draw()
if 'ax_rk' in globals():
    P_c_rk_plot = P_c_rk
    V_c_rk_plot = V_c_rk
    ax_rk.errorbar(V_c_rk_plot, P_c_rk_plot, xerr=delt_V_c_rk, yerr=delt_P_c_rk,
                   color='black', ecolor='black', capsize=2,
                   linestyle='none', zorder=6)
    ax_rk.scatter(V_c_rk_plot, P_c_rk_plot, color='black', marker='o', s=40, label='Punto crítico', zorder=7)
    ax_rk.legend(fontsize=12, ncol=2, loc='upper right')
    ax_rk.figure.canvas.draw()
plt.show()

## CREO ARRAY DE DATOS PARA ACCEDER MEJOR A ELLOS:

Temp=np.array([T1,T2,T3,T4,T5,T6])
Aes=np.array([A1,A2,A3,A4,A5,A6])
delt_Aes=np.array([delt_A1,delt_A2,delt_A3,delt_A4,delt_A5,delt_A6])
Bs=np.array([B1,B2,B3,B4,B5,B6])
delt_Bs=np.array([delt_B1,delt_B2,delt_B3,delt_B4,delt_B5,delt_B6])
Cs=np.array([C1,C2,C3,C4,C5,C6])
delt_Cs=np.array([delt_C1,delt_C2,delt_C3,delt_C4,delt_C5,delt_C6])
Erres=np.array([R2_1,R2_2,R2_3,R2_4,R2_5,R2_6])
Erres_ajustadas=np.array([R2_1_ajustado,R2_2_ajustado,R2_3_ajustado,R2_4_ajustado,R2_5_ajustado,R2_6_ajustado])
Ji_red=np.array([chi_1_reducido,chi_2_reducido,chi_3_reducido,chi_4_reducido,chi_5_reducido,chi_6_reducido])
Sigmas=np.array([sigma_modelo_ajuste_1,sigma_modelo_ajuste_2,sigma_modelo_ajuste_3,sigma_modelo_ajuste_4,sigma_modelo_ajuste_5,sigma_modelo_ajuste_6])

print('Temperaturas:',Temp,'A',Aes,'B:',Bs,'C:',Cs,'R2:',Erres,'R2_aj:',Erres_ajustadas,'CHI:',Ji_red,'Sigma:',Sigmas,'Delta_A:',delt_Aes,'Delta_b:',delt_Bs,'Delta_C:',delt_Cs)

Aes_prima=np.array([A1_prima,A2_prima,A3_prima,A4_prima,A5_prima,A6_prima])
delt_Aes_prima=np.array([delt_A1_prima,delt_A2_prima,delt_A3_prima,delt_A4_prima,delt_A5_prima,delt_A6_prima])
Bs_prima=np.array([B1_prima,B2_prima,B3_prima,B4_prima,B5_prima,B6_prima])
delt_Bs_prima=np.array([delt_B1_prima,delt_B2_prima,delt_B3_prima,delt_B4_prima,delt_B5_prima,delt_B6_prima])
Cs_prima=np.array([C1_prima,C2_prima,C3_prima,C4_prima,C5_prima,C6_prima])
delt_Cs_prima=np.array([delt_C1_prima,delt_C2_prima,delt_C3_prima,delt_C4_prima,delt_C5_prima,delt_C6_prima])
Erres_prima=np.array([R2_1_prima,R2_2_prima,R2_3_prima,R2_4_prima,R2_5_prima,R2_6_prima])
Erres_ajustadas_prima=np.array([R2_1_ajustado_prima,R2_2_ajustado_prima,R2_3_ajustado_prima,R2_4_ajustado_prima,R2_5_ajustado_prima,R2_6_ajustado_prima])
Ji_red_prima=np.array([chi_1_reducido_prima,chi_2_reducido_prima,chi_3_reducido_prima,chi_4_reducido_prima,chi_5_reducido_prima,chi_6_reducido_prima])
Sigmas_prima=np.array([sigma_modelo_ajuste_1_prima,sigma_modelo_ajuste_2_prima,sigma_modelo_ajuste_3_prima,sigma_modelo_ajuste_4_prima,sigma_modelo_ajuste_5_prima,sigma_modelo_ajuste_6_prima])

print('Temperaturas:',Temp,'A',Aes_prima,'B:',Bs_prima,'C:',Cs_prima,'R2:',Erres_prima,'R2_aj:',Erres_ajustadas_prima,'CHI:',Ji_red_prima,'Sigma:',Sigmas_prima,'Delta_A:',delt_Aes_prima,'Delta_b:',delt_Bs_prima,'Delta_C:',delt_Cs_prima)

aes=np.array([a1,a2,a3,a4,a5,a6])
delt_aes=np.array([delt_a1,delt_a2,delt_a3,delt_a4,delt_a5,delt_a6])
bs=np.array([b1,b2,b3,b4,b5,b6])
delt_bs=np.array([delt_b1,delt_b2,delt_b3,delt_b4,delt_b5,delt_b6])
ns=np.array([n1,n2,n3,n4,n5,n6])
delt_ns=np.array([delt_n1,delt_n2,delt_n3,delt_n4,delt_n5,delt_n6])
medias=np.array([a_vdw,b_vdw,n_vdw])
delt_medias=np.array([delt_a_vdw,delt_b_vdw,delt_n_vdw])

print(['aes:',aes,'delt_aes:',delt_aes,'bs:',bs,'delt_bs:',delt_bs,'ns:',ns,'delt_ns:',delt_ns,'medias:',medias,'delt_medias:',delt_medias])

aes_prima=np.array([a1_prima,a2_prima,a3_prima,a4_prima,a5_prima,a6_prima])
delt_aes_prima=np.array([delt_a1_prima,delt_a2_prima,delt_a3_prima,delt_a4_prima,delt_a5_prima,delt_a6_prima])
bs_prima=np.array([b1_prima,b2_prima,b3_prima,b4_prima,b5_prima,b6_prima])
delt_bs_prima=np.array([delt_b1_prima,delt_b2_prima,delt_b3_prima,delt_b4_prima,delt_b5_prima,delt_b6_prima])
ns_prima=np.array([n1_prima,n2_prima,n3_prima,n4_prima,n5_prima,n6_prima])
delt_ns_prima=np.array([delt_n1_prima,delt_n2_prima,delt_n3_prima,delt_n4_prima,delt_n5_prima,delt_n6_prima])
medias_prima=np.array([a_rk,b_rk,n_rk])
delt_medias_prima=np.array([delt_a_rk,delt_b_rk,delt_n_rk])

print(['aes:',aes_prima,'delt_aes:',delt_aes_prima,'bs:',bs_prima,'delt_bs:',delt_bs_prima,'ns:',ns_prima,'delt_ns:',delt_ns_prima,'medias:',medias_prima,'delt_medias:',delt_medias_prima])