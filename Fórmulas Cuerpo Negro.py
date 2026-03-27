#%%

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import t as tst

def mín_cuadrados_LFI(x,y):
# =============================================================================
#   n es el número de medidas experimentales
# =============================================================================
    n=x.shape[0]
    t=tst.ppf(1-0.05/2,n-2)
# =============================================================================
#   Calculo la suma de los datos
# =============================================================================
    x_suma=0
    y_suma=0
    for i in range(n):
        x_suma+=x[i]
        y_suma+=y[i]
# =============================================================================
#   Calculo el promedio de los datos
# =============================================================================
    x_promedio=x_suma/n
    y_promedio=y_suma/n
# =============================================================================
#   Calculo el elemento E para sacar la pendiente
# =============================================================================
    E=0    
    for i in range(n):
        E+=x[i]*y[i]
    E=E-n*y_promedio*x_promedio
# =============================================================================
#   Calculo el elemento D para sacar la pendiente
# =============================================================================
    D=0
    for i in range(n):
        D+=x[i]**2
    D=D-n*x_promedio**2
# =============================================================================
#   Calculo la pendiente de la recta de regresión
# =============================================================================
    m=E/D
# =============================================================================
#   Calculo la ordenada en el origen de la recta de regresión
# =============================================================================
    c=y_promedio-m*x_promedio
# =============================================================================
#   Calculo las incertidumbres
# =============================================================================
    s_res_2=0
    for i in range(n):
        s_res_2+=(y[i]-m*x[i]-c)**2
    s_res_2=s_res_2*(1/(n-2))
    s_m_2=s_res_2/D
    s_c_2=s_res_2*((1/n)+((x_promedio**2)/D))
    delt_m=t*np.sqrt(s_m_2)
    delt_c=t*np.sqrt(s_c_2)
    return m,c,delt_m,delt_c

I=np.array([0.65,0.90,1.00,1.08,1.16,1.24,1.30,1.37,1.46,1.51,1.53,1.59,1.66,1.72,1.77,1.84,1.89,
            1.95,2.00,2.03,2.09,2.14,2.19,2.23,2.28,2.31,2.36,2.41,2.46,2.49,2.54,2.57,
            2.62,2.66,2.70,2.75,2.79,2.85,2.90,2.93,2.99])
V=np.array([0.273,0.599,0.888,1.16,1.48,1.79,2.02,2.29,2.66,2.88,3.00,3.27,3.58,3.85,4.10,4.43,
            4.69,5.00,5.28,5.43,5.76,6.00,6.29,6.53,6.78,7.01,7.29,7.56,7.86,8.05,8.37,8.58,
            8.91,9.19,9.45,9.76,10.05,10.46,10.79,11.02,11.44])
I_2=np.array([0.34,0.74,0.94,1.05,1.14,1.24,1.30,1.38,1.47,1.55,1.60,1.66,1.71,1.75,1.83,1.89,1.94,
            1.99,2.01,2.07,2.11,2.17,2.22,2.28,2.33,2.37,2.41,2.45,2.48,2.52,2.57,2.63,2.68,2.73,2.77,
            2.78,2.82,2.87,2.92,2.96,3])
V_2=np.array([0.11,0.364,0.691,1.04,1.39,1.75,2.01,2.34,2.71,3.04,3.28,3.54,3.79,4,4.39,4.68,4.91,5.17,
            5.31,5.60,5.85,6.16,6.44,6.80,7.12,7.35,7.59,7.85,8.02,8.30,8.57,8.94,9.33,9.60,9.89,10.00,
            10.27,10.60,10.93,11.22,11.49])
delt_V=np.array([0.001,0.001,0.001,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,
            0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,
            0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01])
delt_I=np.array([0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,
            0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,
            0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01])

print(np.shape(I),np.shape(V),np.shape(delt_V),np.shape(delt_I))
print(np.shape(I_2),np.shape(V_2),np.shape(delt_V),np.shape(delt_I))

def En(I,V):
    E=I*V
    delt_E=np.sqrt((delt_I*V)**2+(delt_V*I)**2)
    return E,delt_E

E,delt_E=En(I,V)
E_2,delt_E_2=En(I_2,V_2)

def Res(V,I,delt_V,delt_I):
    R=V/I
    delt_R=np.sqrt((delt_V/I)**2+(delt_I*V/I**2)**2)
    return R,delt_R

R,delt_R=Res(V,I,delt_V,delt_I)
R_2,delt_R_2=Res(V_2,I_2,delt_V,delt_I)


def Temp_1(R,delt_R,R_300=0.3):
    T=50.25+260.7*(R/R_300)-11*(R/R_300)**2
    delt_T=delt_R*np.abs(260.7/R_300-2*11*R/R_300**2)
    return T,delt_T

T_1,delt_T1=Temp_1(R[0],delt_R[0],R_300=0.3)
T_12,delt_T12=Temp_1(R_2[0],delt_R_2[0],R_300=0.3)

def Temp_2(R,delt_R,R_300=0.3):
    T=47.657+258.72*(R/R_300)**0.9-2.7381*(R/R_300)**1.6
    delt_T=delt_R*np.abs(232.848/((R**0.1)*(R_300**0.9))-4.38096*(R**0.6/R_300**1.6))
    return T,delt_T

T_2,delt_T2=Temp_2(R[0:40],delt_R[0:40],R_300=0.3)
T_22,delt_T22=Temp_2(R_2[0:40],delt_R_2[0:40],R_300=0.3)

def lamb_max(T,delt_T,b=2.8978*10**-3):
    lamb_max=b/T
    delt_lamb_max=np.abs(b*delt_T/T**2)
    return lamb_max,delt_lamb_max

lamb_max_1,delt_lamb_max_1=lamb_max(T_1,delt_T1)
lamb_max_2,delt_lamb_max_2=lamb_max(T_2,delt_T2)

lamb_max_12,delt_lamb_max_12=lamb_max(T_12,delt_T12)
lamb_max_22,delt_lamb_max_22=lamb_max(T_22,delt_T22)

# VOLTAJE FRENTE A CORRIENTE

plt.figure()
plt.tick_params(axis='both', labelsize=20)
plt.errorbar(I,V,yerr=delt_V,xerr=delt_I,fmt='o',color='red',
             ecolor='red',capsize=5,label='_nolegend_')
plt.scatter(I,V,color='red',label=r'Datos experimentales')
plt.xlabel(r'$I\;(\mathrm{A})$', fontsize=25)
plt.ylabel(r'$V\;(\mathrm{V})$', fontsize=25)
plt.grid('True')
plt.legend(fontsize=20)
plt.gca().set_facecolor('aliceblue')
plt.show()

plt.figure()
plt.tick_params(axis='both', labelsize=20)
plt.errorbar(I_2,V_2,yerr=delt_V,xerr=delt_I,fmt='o',color='red',
             ecolor='red',capsize=5,label='_nolegend_')
plt.scatter(I_2,V_2,color='red',label=r'Datos experimentales')
plt.xlabel(r'$I\;(\mathrm{A})$', fontsize=25)
plt.ylabel(r'$V\;(\mathrm{V})$', fontsize=25)
plt.grid('True')
plt.legend(fontsize=20)
plt.gca().set_facecolor('aliceblue')
plt.show()

# RESISTENCIA FRENTE A LA TEMPERATURA

plt.figure()
plt.tick_params(axis='both', labelsize=20)
plt.errorbar(T_1,R[0],yerr=delt_R[0],xerr=delt_T1,fmt='o',color='red',
             ecolor='red',capsize=5,label='_nolegend_')
plt.errorbar(T_2,R[0:40],yerr=delt_R[0:40],xerr=delt_T2,fmt='o',color='red',
             ecolor='red',capsize=5,label='_nolegend_')
plt.scatter(T_2,R[0:40],color='red',label=r'Datos experimentales')
plt.xlabel(r'$T\;(\mathrm{K})$', fontsize=25)
plt.ylabel(r'$R\;(\mathrm{\Omega})$', fontsize=25)
plt.grid('True')
plt.legend(fontsize=20)
plt.gca().set_facecolor('aliceblue')
plt.show()

plt.figure()
plt.tick_params(axis='both', labelsize=20)
plt.errorbar(T_12,R[0],yerr=delt_R[0],xerr=delt_T12,fmt='o',color='red',
             ecolor='red',capsize=5,label='_nolegend_')
plt.errorbar(T_22,R[0:40],yerr=delt_R[0:40],xerr=delt_T22,fmt='o',color='red',
             ecolor='red',capsize=5,label='_nolegend_')
plt.scatter(T_22,R[0:40],color='red',label=r'Datos experimentales')
plt.xlabel(r'$T\;(\mathrm{K})$', fontsize=25)
plt.ylabel(r'$R\;(\mathrm{\Omega})$', fontsize=25)
plt.grid('True')
plt.legend(fontsize=20)
plt.gca().set_facecolor('aliceblue')
plt.show()

# ENERGÍA FRENTE A LA TEMPERATURA

m,c,delt_m,delt_c=mín_cuadrados_LFI(np.log(T_2),np.log(E))
log_T2=np.linspace(np.min(np.log(T_2)),np.max(np.log(T_2)),1000)
log_E_fit=m*log_T2+c

plt.figure()
plt.tick_params(axis='both', labelsize=20)
plt.errorbar(np.log(T_1),np.log(E)[0],yerr=delt_E[0]/E[0],xerr=delt_T1/T_1,fmt='o',color='red',
             ecolor='red',capsize=5,label='_nolegend_')
plt.errorbar(np.log(T_2),np.log(E)[0:40],yerr=delt_E[0:40]/E[0:40],xerr=delt_T2/T_2,fmt='o',color='red',
             ecolor='red',capsize=5,label='_nolegend_')
plt.scatter(np.log(T_2),np.log(E)[0:40],color='red',label=r'Datos experimentales')
plt.plot(log_T2,log_E_fit,color='blue',label=r'Ajuste lineal',zorder=3)
plt.xlabel(r'$log(T)$', fontsize=25)
plt.ylabel(r'$log(E)$', fontsize=25)
plt.grid('True')
plt.legend(fontsize=20)
plt.gca().set_facecolor('aliceblue')
plt.show()

print(f'La pendiente de la primera recta de regresión mala es: {m} ± {delt_m} que es la potencia de la temperatura.')

m_2,c_2,delt_m_2,delt_c_2=mín_cuadrados_LFI(np.log(T_2)[29:39],np.log(E)[30:40])
log_T2_2=np.linspace(np.log(T_2)[29],np.max(np.log(T_2)),1000)
log_E_fit_2=m_2*log_T2_2+c_2

plt.figure()
plt.tick_params(axis='both', labelsize=20)
plt.errorbar(np.log(T_1),np.log(E)[0],yerr=delt_E[0]/E[0],xerr=delt_T1/T_1,fmt='o',color='red',
             ecolor='red',capsize=5,label='_nolegend_')
plt.errorbar(np.log(T_2),np.log(E)[0:40],yerr=delt_E[0:40]/E[0:40],xerr=delt_T2/T_2,fmt='o',color='red',
             ecolor='red',capsize=5,label='_nolegend_')
plt.scatter(np.log(T_2),np.log(E)[0:40],color='red',label=r'Datos experimentales')
plt.plot(log_T2_2,log_E_fit_2,color='blue',label=r'Ajuste lineal',zorder=3)
plt.xlabel(r'$log(T)$', fontsize=25)
plt.ylabel(r'$log(E)$', fontsize=25)
plt.grid('True')
plt.legend(fontsize=20)
plt.gca().set_facecolor('aliceblue')
plt.show()

print(f'La pendiente de la primera recta de regresión buena es: {m_2} ± {delt_m_2} que es la potencia de la temperatura.')

m_21,c_21,delt_m_21,delt_c_21=mín_cuadrados_LFI(np.log(T_22),np.log(E_2)[0:40])
log_T22_1=np.linspace(np.log(T_22)[0],np.max(np.log(T_22)),1000)
log_E_2_fit_1=m_21*log_T22_1+c_21

plt.figure()
plt.tick_params(axis='both', labelsize=20)
plt.errorbar(np.log(T_12),np.log(E_2)[0],yerr=delt_E_2[0]/E_2[0],xerr=delt_T12/T_12,fmt='o',color='red',
             ecolor='red',capsize=5,label='_nolegend_')
plt.errorbar(np.log(T_22),np.log(E_2)[0:40],yerr=delt_E_2[0:40]/E_2[0:40],xerr=delt_T22/T_22,fmt='o',color='red',
             ecolor='red',capsize=5,label='_nolegend_')
plt.scatter(np.log(T_22),np.log(E_2)[0:40],color='red',label=r'Datos experimentales')
plt.plot(log_T22_1,log_E_2_fit_1,'b-',label=r'Ajuste lineal',zorder=3)
plt.xlabel(r'$log(T)$', fontsize=25)
plt.ylabel(r'$log(E)$', fontsize=25)
plt.grid('True')
plt.legend(fontsize=20)
plt.gca().set_facecolor('aliceblue')
plt.show()

print(f'La pendiente de la segunda recta de regresión mala es: {m_21} ± {delt_m_21} que es la potencia de la temperatura.')

m_22,c_22,delt_m_22,delt_c_22=mín_cuadrados_LFI(np.log(T_22)[29:39],np.log(E_2)[30:40])
log_T22_2=np.linspace(np.log(T_22)[29],np.max(np.log(T_22)),1000)
log_E_2_fit_2=m_22*log_T22_2+c_22

plt.figure()
plt.tick_params(axis='both', labelsize=20)
plt.errorbar(np.log(T_12),np.log(E_2)[0],yerr=delt_E_2[0]/E_2[0],xerr=delt_T12/T_12,fmt='o',color='red',
             ecolor='red',capsize=5,label='_nolegend_')
plt.errorbar(np.log(T_22),np.log(E_2)[0:40],yerr=delt_E_2[0:40]/E_2[0:40],xerr=delt_T22/T_22,fmt='o',color='red',
             ecolor='red',capsize=5,label='_nolegend_')
plt.scatter(np.log(T_22),np.log(E_2)[0:40],color='red',label=r'Datos experimentales')
plt.plot(log_T22_2,log_E_2_fit_2,'b-',label=r'Ajuste lineal',zorder=3)
plt.xlabel(r'$log(T)$', fontsize=25)
plt.ylabel(r'$log(E)$', fontsize=25)
plt.grid('True')
plt.legend(fontsize=20)
plt.gca().set_facecolor('aliceblue')
plt.show()

print(f'La pendiente de la segunda recta de regresión buena es: {m_22} ± {delt_m_22} que es la potencia de la temperatura.')

# LONGITUD DE ONDA DE MÁXIMA EMISIÓN FRENTE A LA TEMPERATURA

plt.figure()
plt.tick_params(axis='both', labelsize=20)
plt.errorbar(T_1,lamb_max_1*1e6,yerr=delt_lamb_max_1*1e6,xerr=delt_T1,fmt='o',color='red',
             ecolor='red',capsize=5,label='_nolegend_')
plt.errorbar(T_2,lamb_max_2*1e6,yerr=delt_lamb_max_2*1e6,xerr=delt_T2,fmt='o',color='red',
             ecolor='red',capsize=5,label='_nolegend_')
plt.scatter(T_2,lamb_max_2*1e6,color='red',label=r'Datos experimentales')
plt.xlabel(r'$T\;(\mathrm{K})$', fontsize=25)
plt.ylabel(r'$\lambda_{\max}\;(\mathrm{\mu m})$', fontsize=25)
plt.grid('True')
plt.legend(fontsize=20)
plt.gca().set_facecolor('aliceblue')
plt.show()

plt.figure()
plt.tick_params(axis='both', labelsize=20)
plt.errorbar(T_12,lamb_max_12*1e6,yerr=delt_lamb_max_12*1e6,xerr=delt_T12,fmt='o',color='red',
             ecolor='red',capsize=5,label='_nolegend_')
plt.errorbar(T_22,lamb_max_22*1e6,yerr=delt_lamb_max_22*1e6,xerr=delt_T22,fmt='o',color='red',
             ecolor='red',capsize=5,label='_nolegend_')
plt.scatter(T_22,lamb_max_22*1e6,color='red',label=r'Datos experimentales')
plt.xlabel(r'$T\;(\mathrm{K})$', fontsize=25)
plt.ylabel(r'$\lambda_{\max}\;(\mathrm{\mu m})$', fontsize=25)
plt.grid('True')
plt.legend(fontsize=20)
plt.gca().set_facecolor('aliceblue')
plt.show()

# RESISTENCIA DE INCANDESCENCIA

V_inc=1.09
I_inc=1.05
delt_V_inc=0.001
delt_I_Inc=0.01

R_inc,delt_R_inc=Res(V_inc,I_inc,delt_V_inc,delt_I_Inc)
T_inc_2,delt_T_inc_2=Temp_2(R_inc,delt_R_inc)

print(f'La resistencia de incandescencia es: {R_inc} ± {delt_R_inc}')
print(f'La temperatura de incandescencia es: {T_inc_2} ± {delt_T_inc_2}')

# TABLAS LATEX (PANDAS)
try:
    import pandas as pd

    def _format_sig(x, sig=2):
        if x == 0:
            return "0"
        s = f"{x:.{sig}g}"
        return s.replace(".", ",")

    def _format_sci_latex(x, sig=2):
        if x == 0:
            return "0"
        s = f"{x:.{sig}e}"
        mantissa, exp = s.split("e")
        mantissa = mantissa.replace(".", ",")
        exp = exp.lstrip("+")
        exp = exp.lstrip("0") or "0"
        return f"{mantissa}~10^{{{exp}}}"

    def _format_pair(val, err, use_sci=False, err_sig=2):
        if use_sci:
            return _format_sci_latex(val, err_sig), _format_sci_latex(err, err_sig)

        err_str = _format_sig(err, err_sig)
        if "e" in err_str or "E" in err_str:
            return _format_sci_latex(val, err_sig), _format_sci_latex(err, err_sig)

        if "," in err_str:
            decimals = len(err_str.split(",")[1])
        else:
            decimals = 0
        val_str = f"{val:.{decimals}f}".replace(".", ",")
        return val_str, err_str

    def _add_hlines(latex_str):
        lines = latex_str.splitlines()
        out = []
        for line in lines:
            out.append(line)
            if line.strip().endswith("\\\\") and not line.strip().startswith("\\"):
                out.append("\\hline")
        return "\n".join(out)

    # 1) Tabla V - dV - I - dI
    use_sci_v = np.any((np.abs(V) > 0) & (np.abs(V) < 1e-3)) or np.any((np.abs(delt_V) > 0) & (np.abs(delt_V) < 1e-3))
    use_sci_i = np.any((np.abs(I) > 0) & (np.abs(I) < 1e-3)) or np.any((np.abs(delt_I) > 0) & (np.abs(delt_I) < 1e-3))

    v_vals, v_errs = [], []
    i_vals, i_errs = [], []
    for vv, dv, ii, di in zip(V, delt_V, I, delt_I):
        v_str, dv_str = _format_pair(vv, dv, use_sci=use_sci_v)
        i_str, di_str = _format_pair(ii, di, use_sci=use_sci_i)
        v_vals.append(v_str)
        v_errs.append(dv_str)
        i_vals.append(i_str)
        i_errs.append(di_str)

    df_v_i = pd.DataFrame({
        r"$V\;(\mathrm{V})$": v_vals,
        r"$\Delta V\;(\mathrm{V})$": v_errs,
        r"$I\;(\mathrm{A})$": i_vals,
        r"$\Delta I\;(\mathrm{A})$": i_errs,
    })

    # 2) Tabla E - dE - R - dR - T - dT - lambda_max - dlambda_max
    # T y dT empiezan con T_1/delt_T1 y luego continúan con T_2/delt_T2
    E_all = E
    delt_E_all = delt_E
    R_all = R
    delt_R_all = delt_R
    T_all = np.concatenate((np.array([T_1]), T_2))
    delt_T_all = np.concatenate((np.array([delt_T1]), delt_T2))
    lamb_all = 1e6 * np.concatenate((np.array([lamb_max_1]), lamb_max_2))
    delt_lamb_all = 1e6 * np.concatenate((np.array([delt_lamb_max_1]), delt_lamb_max_2))

    use_sci_E = np.any((np.abs(E_all) > 0) & (np.abs(E_all) < 1e-3)) or np.any((np.abs(delt_E_all) > 0) & (np.abs(delt_E_all) < 1e-3))
    use_sci_R = np.any((np.abs(R_all) > 0) & (np.abs(R_all) < 1e-3)) or np.any((np.abs(delt_R_all) > 0) & (np.abs(delt_R_all) < 1e-3))
    use_sci_T = np.any((np.abs(T_all) > 0) & (np.abs(T_all) < 1e-3)) or np.any((np.abs(delt_T_all) > 0) & (np.abs(delt_T_all) < 1e-3))
    use_sci_L = np.any((np.abs(lamb_all) > 0) & (np.abs(lamb_all) < 1e-3)) or np.any((np.abs(delt_lamb_all) > 0) & (np.abs(delt_lamb_all) < 1e-3))

    E_vals, E_errs = [], []
    R_vals, R_errs = [], []
    T_vals, T_errs = [], []
    L_vals, L_errs = [], []

    for e, de, r, dr, t, dt, lmb, dlmb in zip(E_all, delt_E_all, R_all, delt_R_all, T_all, delt_T_all, lamb_all, delt_lamb_all):
        e_str, de_str = _format_pair(e, de, use_sci=use_sci_E)
        r_str, dr_str = _format_pair(r, dr, use_sci=use_sci_R)
        t_str, dt_str = _format_pair(t, dt, use_sci=use_sci_T)
        l_str, dl_str = _format_pair(lmb, dlmb, use_sci=use_sci_L)
        E_vals.append(e_str)
        E_errs.append(de_str)
        R_vals.append(r_str)
        R_errs.append(dr_str)
        T_vals.append(t_str)
        T_errs.append(dt_str)
        L_vals.append(l_str)
        L_errs.append(dl_str)

    df_ertl = pd.DataFrame({
        r"$E\;(\mathrm{W})$": E_vals,
        r"$\Delta E\;(\mathrm{W})$": E_errs,
        r"$R\;(\Omega)$": R_vals,
        r"$\Delta R\;(\Omega)$": R_errs,
        r"$T\;(\mathrm{K})$": T_vals,
        r"$\Delta T\;(\mathrm{K})$": T_errs,
        r"$\lambda_{\max}\;(\mathrm{\\mu m})$": L_vals,
        r"$\Delta \lambda_{\max}\;(\mathrm{\\mu m})$": L_errs,
    })

    latex_v_i = _add_hlines(df_v_i.to_latex(index=False))
    latex_ertl = _add_hlines(df_ertl.to_latex(index=False))

    print("\nTabla LaTeX - V, dV, I, dI:\n")
    print(latex_v_i)
    print("\nTabla LaTeX - E, dE, R, dR, T, dT, lambda_max, dlambda_max:\n")
    print(latex_ertl)

    # 3) Segunda curva: V_2 - dV - I_2 - dI
    use_sci_v2 = np.any((np.abs(V_2) > 0) & (np.abs(V_2) < 1e-3)) or np.any((np.abs(delt_V) > 0) & (np.abs(delt_V) < 1e-3))
    use_sci_i2 = np.any((np.abs(I_2) > 0) & (np.abs(I_2) < 1e-3)) or np.any((np.abs(delt_I) > 0) & (np.abs(delt_I) < 1e-3))

    v2_vals, v2_errs = [], []
    i2_vals, i2_errs = [], []
    for vv, dv, ii, di in zip(V_2, delt_V, I_2, delt_I):
        v_str, dv_str = _format_pair(vv, dv, use_sci=use_sci_v2)
        i_str, di_str = _format_pair(ii, di, use_sci=use_sci_i2)
        v2_vals.append(v_str)
        v2_errs.append(dv_str)
        i2_vals.append(i_str)
        i2_errs.append(di_str)

    df_v_i_2 = pd.DataFrame({
        r"$V\;(\mathrm{V})$": v2_vals,
        r"$\Delta V\;(\mathrm{V})$": v2_errs,
        r"$I\;(\mathrm{A})$": i2_vals,
        r"$\Delta I\;(\mathrm{A})$": i2_errs,
    })

    # 4) Segunda curva: E_2 - dE_2 - R_2 - dR_2 - T_12/T_22 - dT_12/dT_22 - lambda_12/22
    E2_all = E_2
    delt_E2_all = delt_E_2
    R2_all = R_2
    delt_R2_all =delt_R_2
    T2_all = np.concatenate((np.array([T_12]), T_22))
    delt_T2_all = np.concatenate((np.array([delt_T12]), delt_T22))
    lamb2_all = 1e6 * np.concatenate((np.array([lamb_max_12]), lamb_max_22))
    delt_lamb2_all = 1e6 * np.concatenate((np.array([delt_lamb_max_12]), delt_lamb_max_22))

    use_sci_E2 = np.any((np.abs(E2_all) > 0) & (np.abs(E2_all) < 1e-3)) or np.any((np.abs(delt_E2_all) > 0) & (np.abs(delt_E2_all) < 1e-3))
    use_sci_R2 = np.any((np.abs(R2_all) > 0) & (np.abs(R2_all) < 1e-3)) or np.any((np.abs(delt_R2_all) > 0) & (np.abs(delt_R2_all) < 1e-3))
    use_sci_T2 = np.any((np.abs(T2_all) > 0) & (np.abs(T2_all) < 1e-3)) or np.any((np.abs(delt_T2_all) > 0) & (np.abs(delt_T2_all) < 1e-3))
    use_sci_L2 = np.any((np.abs(lamb2_all) > 0) & (np.abs(lamb2_all) < 1e-3)) or np.any((np.abs(delt_lamb2_all) > 0) & (np.abs(delt_lamb2_all) < 1e-3))

    E2_vals, E2_errs = [], []
    R2_vals, R2_errs = [], []
    T2_vals, T2_errs = [], []
    L2_vals, L2_errs = [], []

    for e, de, r, dr, t, dt, lmb, dlmb in zip(E2_all, delt_E2_all, R2_all, delt_R2_all, T2_all, delt_T2_all, lamb2_all, delt_lamb2_all):
        e_str, de_str = _format_pair(e, de, use_sci=use_sci_E2)
        r_str, dr_str = _format_pair(r, dr, use_sci=use_sci_R2)
        t_str, dt_str = _format_pair(t, dt, use_sci=use_sci_T2)
        l_str, dl_str = _format_pair(lmb, dlmb, use_sci=use_sci_L2)
        E2_vals.append(e_str)
        E2_errs.append(de_str)
        R2_vals.append(r_str)
        R2_errs.append(dr_str)
        T2_vals.append(t_str)
        T2_errs.append(dt_str)
        L2_vals.append(l_str)
        L2_errs.append(dl_str)

    df_ertl_2 = pd.DataFrame({
        r"$E\;(\mathrm{W})$": E2_vals,
        r"$\Delta E\;(\mathrm{W})$": E2_errs,
        r"$R\;(\Omega)$": R2_vals,
        r"$\Delta R\;(\Omega)$": R2_errs,
        r"$T\;(\mathrm{K})$": T2_vals,
        r"$\Delta T\;(\mathrm{K})$": T2_errs,
        r"$\lambda_{\max}\;(\mathrm{\\mu m})$": L2_vals,
        r"$\Delta \lambda_{\max}\;(\mathrm{\\mu m})$": L2_errs,
    })

    latex_v_i_2 = _add_hlines(df_v_i_2.to_latex(index=False))
    latex_ertl_2 = _add_hlines(df_ertl_2.to_latex(index=False))

    print("\nTabla LaTeX - V_2, dV, I_2, dI:\n")
    print(latex_v_i_2)
    print("\nTabla LaTeX - E_2, dE_2, R_2, dR_2, T_2, dT_2, lambda_max, dlambda_max:\n")
    print(latex_ertl_2)
except ImportError:
    print("\nPandas no está instalado. Instálalo con: pip install pandas\n")
