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


I=np.array([181.2,165.6,147.2,133.6,118.8,103.1,89,76.3,64.9,54.8,43.2,35.2,26.7,19.8,17.3,14.3,11.5,8.9,7.2,5.2,3.7,
            2.0,1.2,0.07,0.03,0.00,-0.03,-0.07,-1.2,-2.0,-3.7,-5.2,-7.0,-7.3])*1e-9
V=np.array([0.052,0.1,0.156,0.200,0.249,0.303,0.355,0.405,0.452,0.499,0.555,0.600,0.653,
            0.702,0.724,0.750,0.777,0.805,0.825,0.851,0.876,0.903,0.920,0.930,0.940,0.946,0.955,0.963,0.977,1,
            1.061,1.144,1.385,1.520])*-1
delt_V=0.001
delt_I=0.1e-9

print(np.shape(I),np.shape(V))

plt.figure()
plt.tick_params(axis='both', labelsize=20)
plt.errorbar(V,I*1e9,yerr=delt_I*1e9,xerr=delt_V,fmt='o',color='red',
             ecolor='red',capsize=5,label='_nolegend_')
plt.scatter(V,I*1e9,color='red',label=r'Datos experimentales')
plt.ylabel(r'$I\;(\mathrm{nA})$', fontsize=25)
plt.xlabel(r'$V\;(\mathrm{V})$', fontsize=25)
plt.grid('True')
plt.legend(fontsize=20)
plt.gca().set_facecolor('aliceblue')
plt.show()

V_f1=np.array([0.255,0.324,0.461,0.664,0.946])*-1
V_f2=np.array([0.244,0.323,0.462,0.665,0.946])*-1
lamb_f1=np.array([656,634.77,592.06,525.96,468.86])*1e-9
lamb_f2=lamb_f1
delt_lamb_f1=5e-9
delt_lamb_f2=delt_lamb_f1
lamb=np.linspace(450,670,1000)*1e-9
nu_f1=3e8/lamb_f1
nu_f2=3e8/lamb_f2
delt_nu_f1=np.abs(3e8/lamb_f1**2*delt_lamb_f1)
delt_nu_f2=np.abs(3e8/lamb_f2**2*delt_lamb_f2)
nu=3e8/lamb

m_1,c_1,delt_m_1,delt_c_1=mín_cuadrados_LFI(nu_f1,V_f1)
m_2,c_2,delt_m_2,delt_c_2=mín_cuadrados_LFI(nu_f2,V_f2)

V_f1_teo=m_1*nu+c_1
V_f2_teo=m_2*nu+c_2

plt.figure()
plt.tick_params(axis='both', labelsize=20)
plt.errorbar(nu_f1,V_f1,yerr=delt_V,xerr=0,fmt='o',color='red',
             ecolor='red',capsize=5,label='_nolegend_')
plt.scatter(nu_f1,V_f1,color='red',label=r'Datos experimentales')
plt.plot(nu,V_f1_teo,color='blue',label=r'Ajuste lineal')
plt.xlabel(r'$\nu\;(\mathrm{Hz})$', fontsize=25)
plt.ylabel(r'$V\;(\mathrm{V})$', fontsize=25)
plt.grid('True')
plt.legend(fontsize=20)
plt.gca().set_facecolor('aliceblue')
plt.show()

plt.figure()
plt.tick_params(axis='both', labelsize=20)
plt.errorbar(nu_f2,V_f2,yerr=delt_V,xerr=0,fmt='o',color='red',
             ecolor='red',capsize=5,label='_nolegend_')
plt.scatter(nu_f2,V_f2,color='red',label=r'Datos experimentales')
plt.plot(nu,V_f2_teo,color='blue',label=r'Ajuste lineal')
plt.xlabel(r'$\nu\;(\mathrm{Hz})$', fontsize=25)
plt.ylabel(r'$V\;(\mathrm{V})$', fontsize=25)
plt.grid('True')
plt.legend(fontsize=20)
plt.gca().set_facecolor('aliceblue')
plt.show()

def h(m,delt_m):
    return m*-1.6e-19,delt_m*1.6e-19

def nu_0(c,delt_c,h,delt_h):
    return -c*-1.6e-19/h,np.sqrt((delt_c*1.6e-19/h)**2+(c*1.6e-19*delt_h/h**2)**2)

h_f1,delt_f1=h(m_1,delt_m_1)
h_f2,delt_f2=h(m_2,delt_m_2)
nu_0_f1,delt_nu_0_f1=nu_0(c_1,delt_c_1,h_f1,delt_f1)
nu_0_f2,delt_nu_0_f2=nu_0(c_2,delt_c_2,h_f2,delt_f2)

print(f'Constante de Planck para altas Potencias: {h_f1} ± {delt_f1} J*s')
print(f'Constante de Planck para bajas Potencias: {h_f2} ± {delt_f2} J*s')

print(f'Frecuencia umbral para altas Potencias: {nu_0_f1} ± {delt_nu_0_f1} Hz')
print(f'Frecuencia umbral para bajas Potencias: {nu_0_f2} ± {delt_nu_0_f2} Hz')

