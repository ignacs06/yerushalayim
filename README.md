# yerushalayim
Para agregar códigos de prácticas en común y de proyectos personales que tengamos

L=162e-3
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import t
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

ts1=t.ppf(1-0.05/2,3)
n1=300/(162e-3)
mu_0=np.pi*4e-7
i=np.array([0.2,0.4,0.6,0.8,1.0])
bt=mu_0*i*n1
error_a=0.01
error_b=1e-5
error_t=2.3e-5
bexp=np.array([0.00049,0.00095,0.00143,0.00191,0.00232])
m1t,c1t,dm1,dc1=min_cuadrados(i,bt,ts1)
m1e,c1e,dme,dce=min_cuadrados(i,bexp,ts1)
plt.figure()
plt.grid('True')
plt.plot(i,bt,'o',ms=3,label='B teórico',color='blue')
plt.plot(i,bexp,'o',ms=3,label='B experimental',color='red')
plt.plot(i,c1t+i*m1t,linestyle='--',label='Regresión B teórico')
plt.plot(i,c1e+i*m1e,linestyle='--',label='Regresión B experimental',color='pink')
plt.errorbar(i, bexp, yerr=error_b, fmt='o', xerr=error_a, ms=1, color='red',ecolor='black')
plt.errorbar(i, bt, yerr=error_t, fmt='o', xerr=error_a, ms=1, color='red',ecolor='black')
plt.legend()
plt.show()



#2 VARIACION NUMERO ESPIRAS

ts2=t.ppf(1-0.05/2,1)
nesp=np.array([75,150,300])
bst=mu_0*1*nesp/(162e-3)
bse=np.array([0.00077,0.00122,0.00232])
errorbst=np.array([5.7e-6,1.1e-5,2.3e-5])
mbst,cbst,dmbst,dcbst=min_cuadrados(nesp,bst,ts2)
mbse,cbse,dmbse,dcbse=min_cuadrados(nesp, bse, ts2)
plt.figure()
plt.grid('True')
plt.plot(nesp,bst,'o',ms=3,color='blue',label='B teórico')
plt.plot(nesp,bse,'o',ms=3,color='red',label='B experimental')
plt.plot(nesp,cbst+nesp*mbst,linestyle='--',color='skyblue',label='Regresión teórica')
plt.plot(nesp,cbse+nesp*mbse,linestyle='--',color='red',label='Regresión experimental')
plt.errorbar(nesp, bst, yerr=errorbst, fmt='o', ms=1, color='red',ecolor='black')
plt.legend()
plt.show()


#3 15 POSIÇAOES CARMEN ES EL ARRAY DE DISTANCIAS RESPECTO DEL CENTRO
bquince=np.array([19e-5,14.7e-4,21.4e-4,22.3e-4,22.8e-4,23e-4,23.1e-4,23.2e-4,23.3e-4,23.1e-4,23e-4,22.7e-4,22e-4,19e-4,2.3e-4])
carmen=np.array([-0.097,-0.075,-0.06,-0.05,-0.04,-0.03,-0.02,0,0.02,0.03,0.04,0.05,0.06,0.075,0.097])
carmenlinspace=np.linspace(-0.097,0.097,10000)
R=13e-3
def B3(x):
    a=x+L/2
    b=x-L/2
    return (a/np.sqrt(R**2+a**2)-b/np.sqrt(R**2+b**2))*mu_0*300/(2*L)
#Samu no te fíes mucho (NADA) de la fórmula de aquí abajo (dB3)
def dB3(x):
    return np.sqrt((B3(carmen)**2)*(0.01**2+(0.0001/0.162)**2+0.0022))
plt.figure()
plt.grid(True)
plt.plot(carmen,bquince,'o',ms=3)
plt.errorbar(carmen, bquince, yerr=dB3(carmen), fmt='o', ms=1, color='red',ecolor='black')
plt.plot(carmenlinspace,B3(carmenlinspace),linestyle='--',color='skyblue')

def ohm(V):
    R=0.3
    I=V/R
    return I

#Ley de Faraday

'''Defino la función que da el voltaje inducido en la bobina
secundaria a partir de la ecuación (7) del guión AYUDA'''
def xi(I,f,N):
    df=10
    dI=1e-3
    dl=1e-3 
    dd=1e-3
    d=0.041
    l=0.75
    w=2*np.pi*f
    Np=485
    A=np.pi*(d/2)**2
    xi=mu_0	*N*A*Np*w*(I/l)
    ixi=mu_0*N*Np*np.sqrt((2*A*np.pi*I*df/l)**2+(2*np.pi*A*f*dI/l)**2+
                          (A*I*f*dl/(l**2))**2+(np.pi**2*f*I*d*dd/l)**2)
    return xi,ixi

'Samu de esta función de mierda (xi) sí te puedes fiar porque llevo veinte minutos'
'para calcular la incertidumbre de 5 putos datos' '''SÍ SAMUEL ESTÁ TODO EN SISTEMA'''
'''INTERNACIONAL'''

'''WE ARE CHARLIE KIRK'''    



I0=np.array([0.0146,0.0180,0.02168,0.02595,0.02988])
I0l=np.linspace(0.0146,0.02988,10000)
Aie=np.array([0.2830,0.3618,0.4226,0.5076,0.5572])
wearecharliekirk=min_cuadrados(I0, Aie, ts1)
plt.figure()
plt.grid(True)
plt.xlabel(r'$I_0$ (A)', fontsize=12)
plt.ylabel(r'$Amplitud_{ind}$ (V)', fontsize=12)
plt.plot(I0l,xi(I0l,1e4,300)[0],color='blue')
plt.plot(I0,xi(I0,1e4,300)[0],ms=3,color='blue',label='Valores teóricos')
plt.plot(I0,Aie,'o',ms=3,color='red',label='Valores experimentales')
plt.plot(I0l,wearecharliekirk[1]+wearecharliekirk[0]*I0l,linestyle='--')
plt.errorbar(I0,xi(I0,1e4,300)[0],xerr=1e-3,yerr=xi(I0,1e4,300)[1],fmt='o',ms=1,color='blue',ecolor='black')
plt.errorbar(I0, Aie, yerr=1e-4, fmt='o', ms=1, color='red',ecolor='black',label='Barras de error')
plt.legend()

#Apartado 5
I5=0.03
f5=np.array([3000,5010,7010,10000,12020])
f5l=np.linspace(3000,12020,1000)
A5e=np.array([0.2075,0.2621,0.4588,0.5076,0.5826])
jeffreyepstein=min_cuadrados(f5,A5e,ts1)
plt.figure()
plt.grid(True)
plt.xlabel(r'$Frecuencia$ (Hz)', fontsize=12)
plt.ylabel(r'$Amplitud_{ind}$ (V)', fontsize=12)
plt.plot(f5l,xi(0.03,f5l,300)[0],color='blue')
plt.plot(f5,xi(0.03,f5,300)[0],ms=3,color='blue',label='Valores teóricos')
plt.plot(f5,A5e,'o',ms=3,color='red',label='Valores experimentales')
plt.plot(f5l,jeffreyepstein[1]+jeffreyepstein[0]*f5l,linestyle='--')
plt.errorbar(f5,xi(0.03,f5,300)[0],xerr=10,yerr=xi(0.03,f5,300)[1],fmt='o',ms=1,color='blue',ecolor='black')
plt.errorbar(f5, A5e, yerr=1e-4, fmt='o', ms=1, color='red',ecolor='black')
plt.errorbar(f5,A5e,xerr=10,yerr=1e-4,fmt='o',ms=1,color='blue',ecolor='purple',label='Barras de error')
plt.legend()

#Apartado 6
f6=1e4
N6l=np.linspace(100,300,200)
N6=np.array([100,200,300])
A6e=np.array([0.2090,0.4079,0.7178])
netanyahu=min_cuadrados(N6,A6e,ts2)
plt.figure()
plt.grid(True)
plt.xlabel(r'Nº de espiras', fontsize=12)
plt.ylabel(r'$Amplitud_{ind}$ (V)', fontsize=12)
plt.plot(N6l,xi(I5,f6,N6l)[0],color='blue')
plt.plot(N6,xi(I5,f6,N6)[0],ms=3,color='blue',label='Valores teóricos')
plt.plot(N6,A6e,'o',ms=3,color='red',label='Valores experimentales')
plt.plot(N6l,netanyahu[1]+netanyahu[0]*N6l,linestyle='--')
plt.errorbar(N6,xi(0.03,f6,N6)[0],xerr=10,yerr=xi(0.03,f6,N6l)[1],fmt='o',ms=1,color='blue',ecolor='black')
plt.errorbar(N6, A6e, yerr=1e-4, fmt='o', ms=1, color='red',ecolor='black')
plt.errorbar(N6,A6e,xerr=1,yerr=1e-4,fmt='o',ms=1,color='blue',ecolor='purple',label='Barras de error')
plt.legend()

