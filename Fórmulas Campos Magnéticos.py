# =============================================================================
# Sesión 3: Campos magnéticos
# =============================================================================

#%%

# =============================================================================
# Librerías y funciones comunes:
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import t
plt.close('all')

def min_cuadrados(x,y,t):
# =============================================================================
#   n es el número de medidas experimentales
# =============================================================================
    n=x.shape[0]
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

# =============================================================================
# 1) Campo magnético de la bobina en función de la corriente:
# =============================================================================

L=162e-3
d=26e-3
N=300
I_teo=np.linspace(0.2,1.0,10000)
I=np.array([0.2,0.4,0.6,0.8,1.0])
B_exp=np.array([])
mu_0=np.pi*4e-7
B_teo=mu_0*I_teo*N/L
delt_I=0
delt_B_exp=0
ts1=t.ppf(1-0.05/2,3)

m1e,c1e,delt_me,delt_ce=min_cuadrados(I,B_exp,ts1)

plt.figure()
plt.tick_params(axis='both',labelsize=20)  
plt.plot(I_teo,B_teo,label='$B_{teórico}$',color='blue')
plt.plot(I_teo,m1e*I_teo+c1e,'-',color='red')
plt.errorbar(I,B_exp,yerr=delt_B_exp,xerr=delt_I,fmt='o',color='red',ecolor='red',capsize=5)
plt.scatter(I,B_exp,c='red',marker='o',label='$B_{experimental}$')
plt.xlabel(r'$I~(A)$',fontsize=25)
plt.ylabel(r'$B~(T)$',fontsize=25)
plt.grid('True')
plt.legend(fontsize=20)
plt.gca().set_facecolor('aliceblue')
plt.show()

print('Ejercicio 1')
print('La pendiente de la regresión es: ',m1e,'±',delt_me)
print('La pendiente de la regresión teórica es: ',mu_0*N/L)

# =============================================================================
# 2) Campo magnético de la bobina en función del número de espiras:
# =============================================================================

N_teo=np.linspace(75,300,10000)
N_2=np.array([75,150,300])
B_teo2=mu_0*1*N_teo/(2*np.sqrt((L/2)**2+(d/2)**2))
B_exp2=np.array([])
ts2=t.ppf(1-0.05/2,1)

m2e,c2e,delt_m2e,delt_c2e=min_cuadrados(N_2,B_exp2,ts2)

plt.figure()
plt.tick_params(axis='both',labelsize=20)  
plt.plot(N_teo,B_teo2,label='$B_{teórico}$',color='blue')
plt.plot(N_teo,m2e*N_teo+c2e,'-',color='red')
plt.errorbar(N_2,B_exp2,yerr=delt_B_exp,xerr=0,fmt='o',color='red',ecolor='red',capsize=5)
plt.scatter(N_2,B_exp2,c='red',marker='o',label='$B_{experimental}$')
plt.xlabel(r'$N$',fontsize=25)
plt.ylabel(r'$B~(T)$',fontsize=25)
plt.grid('True')
plt.legend(fontsize=20)
plt.gca().set_facecolor('aliceblue')
plt.show()

print('Ejercicio 2')
print('La pendiente de la regresión es: ',m2e,'±',delt_m2e)
print('La pendiente de la regresión teórica es: ',mu_0*1/(2*np.sqrt((L/2)**2+(d/2)**2)))

# =============================================================================
# 3) Campo magnético de la bobina a lo largo de su eje:
# =============================================================================

B_exp3=np.array([])
x=np.array([])
delt_x=0
x_teo=np.linspace(x[0],x[-1],10000)
a=x_teo+L/2
b=x_teo-L/2
B_teo3=mu_0*1*N/(2*L)*(a/(np.sqrt((d/2)**2+a**2))-b/(np.sqrt((d/2)**2+b**2)))

plt.figure()
plt.tick_params(axis='both',labelsize=20)  
plt.plot(x_teo,B_teo3,label='$B_{teórico}$',color='blue')
plt.errorbar(x,B_exp3,yerr=delt_B_exp,xerr=delt_x,fmt='o',color='red',ecolor='red',capsize=5)
plt.scatter(x,B_exp3,c='red',marker='o',label='$B_{experimental}$')
plt.xlabel(r'$x~(m)$',fontsize=25)
plt.ylabel(r'$B~(T)$',fontsize=25)
plt.grid('True')
plt.legend(fontsize=20)
plt.gca().set_facecolor('aliceblue')
plt.show()

# =============================================================================
# 4) Determinación de la fem inducida en función del campo magnético:
# =============================================================================

def Ohm(I,delt_V,R=10):
    V=I*R
    delt_I=np.sqrt((delt_V/R)**2)
    return V,delt_I    

I_R=np.linspace(10e-3,30e-3,5)
delt_V_R=0
V_R,delt_IR=Ohm(I_R,delt_V_R)
fem_exp=np.array([])
delt_fem_exp=0
f_4=10e3
delt_f4=0
w_4=2*np.pi*f_4
delt_w4=2*np.pi*delt_f4
N_sec=300
A=np.pi*(0.5*41e-3)**2
N_prima=0
l_prima=0
delt_l=0

def fem_teo(N,A,N_prima,w,I,l,delt_I,delt_l,delt_w):
    fem_teo=mu_0*N*A*N_prima*w*I/l
    delt_fem_teo=np.sqrt((mu_0*N*A*N_prima*w*delt_I/l)**2+(mu_0*N*A*N_prima*w*I*delt_l/l**2)**2+(mu_0*N*A*N_prima*delt_w*I/l)**2)
    return fem_teo,delt_fem_teo

I_teo4=np.linspace(10e-3,30e-3,10000)
fem_teo4,delt_fem_teo4=fem_teo(N_sec,A,N_prima,w_4,I_teo4,l_prima,delt_IR,delt_l)

m3e,c3e,delt_m3e,delt_c3e=min_cuadrados(I_R,fem_exp,ts1)

plt.figure()
plt.tick_params(axis='both',labelsize=20)  
plt.plot(I_teo4,fem_teo4,label='$\epsilon_{teórico}$',color='blue')
plt.plot(I_teo4,m3e*I_teo4+c3e,'-',color='red')
plt.errorbar(I_R,fem_exp,yerr=delt_fem_exp,xerr=delt_IR,fmt='o',color='red',ecolor='red',capsize=5)
plt.scatter(I_R,fem_exp,c='red',marker='o',label='$\epsilon_{experimental}$')
plt.xlabel(r'$I~(A)$',fontsize=25)
plt.ylabel(r'$\epsilon~(V)$',fontsize=25)
plt.grid('True')
plt.legend(fontsize=20)
plt.gca().set_facecolor('aliceblue')
plt.show()

print('Ejercicio 4')
print('La pendiente de la regresión es: ',m3e,'±',delt_m3e) 
print('La pendiente de la regresión teórica es: ',mu_0*N_sec*A*N_prima*w_4/l_prima)

# =============================================================================
# 5) Determinación de la fem inducida en función de la frecuencia:
# =============================================================================

I_5=30e-3
V_5,delt_I_5=Ohm(I_5,delt_V_R)
f_5=np.linspace(1000,12000,5)
delt_f_5=0
w_5=2*np.pi*f_5
delt_w_5=2*np.pi*delt_f_5
fem_exp5=np.array([])
delt_fem_exp5=0
f_teo=np.linspace(1000,12000,10000)
w_teo=2*np.pi*f_teo
fem_teo5,delt_fem_teo5=fem_teo(N_sec,A,N_prima,w_teo,I_5,l_prima,delt_I_5,delt_l,delt_w_5)
m4e,c4e,delt_m4e,delt_c4e=min_cuadrados(f_5,fem_exp5,ts1)

plt.figure()
plt.tick_params(axis='both',labelsize=20)  
plt.plot(f_teo,fem_teo5,label='$\epsilon_{teórico}$',color='blue')
plt.plot(f_teo,m4e*f_teo+c4e,'-',color='red')
plt.errorbar(f_5,fem_exp5,yerr=delt_fem_exp5,xerr=delt_f_5,fmt='o',color='red',ecolor='red',capsize=5)
plt.scatter(f_5,fem_exp5,c='red',marker='o',label='$\epsilon_{experimental}$')
plt.xlabel(r'$f~(Hz)$',fontsize=25)
plt.ylabel(r'$\epsilon~(V)$',fontsize=25)
plt.grid('True')
plt.legend(fontsize=20)
plt.gca().set_facecolor('aliceblue')
plt.show()

print('Ejercicio 5')
print('La pendiente de la regresión es: ',m4e,'±',delt_m4e)
print('La pendiente de la regresión teórica es: ',mu_0*N_sec*A*N_prima*I_5*l_prima**(-1)*2*np.pi)

# =============================================================================
# 5) Determinación de la fem inducida en función de la frecuencia:
# =============================================================================

f_6=1e4
delt_f6=0
w_6=2*np.pi*f_6
delt_w6=2*np.pi*delt_f6
N_6=np.array([100,200,300])
fem_exp6=np.array([])
delt_fem_exp6=0
N_teo=np.linspace(100,300,10000)
fem_teo6,delt_fem_teo6=fem_teo(N_teo,A,N_prima,w_6,I_5,l_prima,delt_I_5,delt_l,delt_w6)
m5e,c5e,delt_m5e,delt_c5e=min_cuadrados(N_6,fem_exp6,ts2)

plt.figure()
plt.tick_params(axis='both',labelsize=20)
plt.plot(N_teo,fem_teo6,label='$\epsilon_{teórico}$',color='blue')
plt.plot(N_teo,m5e*N_teo+c5e,'-',color='red')
plt.errorbar(N_6,fem_exp6,yerr=delt_fem_exp6,xerr=0,fmt='o',color='red',ecolor='red',capsize=5)
plt.scatter(N_6,fem_exp6,c='red',marker='o',label='$\epsilon_{experimental}$')
plt.xlabel(r'$N$',fontsize=25)
plt.ylabel(r'$\epsilon~(V)$',fontsize=25)
plt.grid('True')
plt.legend(fontsize=20)
plt.gca().set_facecolor('aliceblue')
plt.show()

print('Ejercicio 6')
print('La pendiente de la regresión es: ',m5e,'±',delt_m5e)
print('La pendiente de la regresión teórica es: ',mu_0*A*N_prima*w_6*I_5/l_prima)