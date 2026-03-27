# -*- coding: utf-8 -*-
"""
Created on Fri Feb 13 02:18:27 2026

@author: Nacho
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

delt_x=0.2
delt_a=0.2
x_1=np.array([2*5])
a_1=np.array([2.2*5])
w_1=np.array([1000])*2*np.pi
R=22e3
C=4.7e-9
V_R=np.array([0.2*1.7,2.2*0.2,3.6*0.2,2.5*0.5,1.9*1,2.9*1,1.9*2,2.3*2,2.5*2,2.6*2])
delt_V=0.2
V_T=5
x=np.array([1.1*5,1.1*5,1.1*5,1.1*5,1*5,0.9*5,0.7*5,0.5*5,0.375*5,0.25*5])
a=np.array([1.1*5,1.1*5,1.1*5,1.2*5,1.1*5,1.1*5,1.1*5,1.05*5,1.05*5,1.05*5])
f=np.array([132,167,278,473,771,1280,2140,3560,5940,10020])
w=2*np.pi*f

def delta_exp(x,a):
    delt=np.arcsin(x/a)*(180/np.pi)
    delt_delt_exp=np.zeros(np.shape(a))
    for i in range(np.shape(a)[0]):
        if a[i]!=x[i]:
            delt_delt_exp[i]=np.sqrt(((180/np.pi)*delt_a*(x[i]/(np.sqrt(1-(x[i]**2/a[i]**2))*a[i]**2)))**2+((180/np.pi)*delt_x/(a[i]*np.sqrt(1-(x[i]**2/a[i]**2))))**2)
    return delt,delt_delt_exp

delt_exp1,delt_delt_exp1=delta_exp(x_1,a_1)
delt_exp,delt_delt_exp=delta_exp(x,a)

def delta_teo(w,R,C):
    delt=np.arctan(1/(w*R*C))*(180/np.pi)
    return delt

delt_teo1=delta_teo(w_1,R,C)

def Ohm(V_R,R):
    I=V_R/R
    delt_I=delt_V/R
    return I,delt_I

I,delt_I=Ohm(V_R,R)

def imped(V_T,I,R):
    Z=(V_T/I)
    delt_Z=np.sqrt((delt_V/I)**2+(V_T*delt_I/I**2)**2)
    return Z,delt_Z

Z,delt_Z=imped(V_T,I,R)
delt_Z[0]=0
delt_Z[1]=0

def Real(Z,delt_exp):
    Re=Z*np.cos(delt_exp/(180/np.pi))
    delt_Re=np.sqrt((delt_Z*np.cos(delt_exp/(180/np.pi)))**2+(Z*np.sin(delt_exp/(180/np.pi))*delt_delt_exp/(180/np.pi))**2)
    return Re,delt_Re

Re,delt_Re=Real(Z,delt_exp)
delt_Re[0]=0
delt_Re[1]=0

def Imag(Z,Re):
    Im=np.sqrt(Z**2-Re**2)
    delt_Im=np.sqrt((Z*delt_Z/np.sqrt(Z**2-Re**2))**2+(Re*delt_Re/np.sqrt(Z**2-Re**2))**2)
    return Im,delt_Im

Im,delt_Im=Imag(Z,Re)
delt_Im[0]=0
delt_Im[1]=0

def capacidad(Im,w):
    C_exp=1/(w*Im)
    delt_C=(delt_Im/(w*Im**2))
    return C_exp,delt_C

C_exp,delt_C=capacidad(Im,w)

f_exp=w/(2*np.pi)
f_teo=np.linspace(60,10000,10000)
w_teo=f_teo*2*np.pi
Z_teo=np.sqrt(R**2+(1/(w_teo*C))**2)
delt_teo=delta_teo(w_teo,R,C)


plt.figure()
plt.tick_params(axis='both',labelsize=20)  
plt.plot(f_teo,Z_teo,label='$Z_{teórica}$',color='blue')
plt.errorbar(f_exp[2:],Z[2:],yerr=delt_Z[2:],xerr=0,fmt='o',color='red',
             ecolor='red',capsize=5)
plt.scatter(f_exp[0:2],Z[0:2],c='red',marker='o',label='$Z_{experimental}$')
plt.xlabel(r'$f~(Hz)$',fontsize=25)
plt.ylabel(r'$|Z|~(\Omega)$',fontsize=25)
plt.xscale('log')
plt.yscale('log')
plt.grid('True')
plt.legend(fontsize=20)
plt.gca().set_facecolor('aliceblue')
plt.show()

plt.figure()
plt.tick_params(axis='both',labelsize=20)  
plt.plot(f_teo,delt_teo,label='$delta_{teórico}$',color='blue')
plt.errorbar(f_exp[2:],delt_exp[2:],yerr=delt_delt_exp[2:],xerr=0,fmt='o',color='red',
             ecolor='red',capsize=5)
plt.scatter(f_exp[0:2],delt_exp[0:2],c='red',marker='o',label='$delta_{experimental}$')
plt.xlabel(r'$f~(Hz)$',fontsize=25)
plt.ylabel(r'$\delta~(^\circ)$',fontsize=25)
plt.xscale('log')
plt.grid('True')
plt.legend(fontsize=20)
plt.gca().set_facecolor('aliceblue')
plt.show()

def tau(R,C):
    tau=R*C
    return tau

C_p2=4.7e-9
R_p2=np.array([1e6,220e3,22e3])
tau_p2=tau(R_p2,C_p2)

C_13=4.7e-9
R_13=np.array([22e3,330])
tau_13=tau(R_13,C_13)

C_16=np.array([4.7e-9,22e-9])
R_16=1e6
tau_16=tau(R_16,C_16)

data = {
    r'$f (Hz)$': f,
    r'$V_{R} (V)$': V_R,
    r'$x$': x,
    r'$a$': a,
    r'$\delta_{exp}(^{\circ})$': delt_exp
}


df = pd.DataFrame(data)

df[ r'$V_{R} (V)$'] = [f"${v:.3g}$" for v in df[ r'$V_{R} (V)$']]
df[r'$x$'] = [f"${v:.3g}$" for v in df[r'$x$']]
df[r'$a$'] = [f"${v:.3g}$" for v in df[r'$a$']]
df[r'$\delta_{exp}(^{\circ})$'] = [f"${v:.1f} \pm {e:.1f}$" for v,
                                   e in zip(df[r'$\delta_{exp}(^{\circ})$'], delt_delt_exp)]

print(df.to_latex(
    index=False,
    caption="a",
    label="tab:cuestion4",
    column_format='|c|c|c|c|c|', 
    position='H',
    bold_rows=True,
    escape=False
))

calc= {
       r'$f (Hz)$': f,
       r'$I (A)$': I,
       r'$|Z| (\Omega)$': Z,
       r'$\Re(Z) (\Omega)$': Re,
       r'$\Im(Z) (\Omega)$': Im
       }

cl= pd.DataFrame(calc)


cl[r'$|Z| (\Omega)$']=[f"${v:.6g}$" for v in Z]
cl[r'$\Im(Z) (\Omega)$']=[f"${v:.6g}$" for v in Im]
cl[r'$\Re(Z) (\Omega)$']=[f"${v:.6g}$" for v in Re]


print(cl.to_latex(
    index=False,
    caption="b",
    label="tab:cuestion4.2",
    column_format='|c|c|c|c|c|', 
    position='H',
    bold_rows=True,
    escape=False
))

capacidad={r'$f (Hz)$':f,
     r'$\Im(Z) (\Omega)$':Im,
     r'$C_{exp} (nF)$': C_exp*10**9     }

cap=pd.DataFrame(capacidad)
cap[r'$C_{exp} (nF)$'] = [f"${v:.1f} \pm {e:.1f}$" for v,
                   e in zip(cap[r'$C_{exp} (nF)$'], delt_C*10**9)]
cap[r'$\Im(Z) (\Omega)$']=[f"${v:.6g}$" for v in Im]
print(cap.to_latex(
    index=False,
    caption="b",
    label="tab:cuestion4.2",
    column_format='|c|c|c|', 
    position='H',
    bold_rows=True,
    escape=False
))

