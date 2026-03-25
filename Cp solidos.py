#%%

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit
import pandas as pd 

cp=75.3/18.02
#J/molK
mc=234.24
Tc=352.06
Tm=316.86
mf=250.8
#g
Tf=293.16



def C(cp,mc,Tc,Tm,mf,Tf):
    a=mc*(Tc-Tm)
    c=Tm-Tf
    b=mf*c
    C=cp*(a-b)/c
    return C

C1=C(cp,mc,Tc,Tm,mf,Tf) #J/K
 

Teb=370.96 

mal=88.62
mfal=246.90
Tfal=294.06
Tmal=298.56

mlat=90.76
mflat=233.78
Tflat=292.46
Tmlat=295.16

mac=89.12
mfac=214.87
Tfac=294.26
Tmac=295.66


def cp_solido(mf2,cp,Tm2,Tf2,C,ms,Teb):
    a=Tm2-Tf2
    b=Teb-Tm2
    cp=((mf2*cp*a)+(C*a))/(ms*b)
    return cp
    
cp_al=cp_solido(mfal,cp,Tmal,Tfal,C1,mal,Teb)
cp_lat=cp_solido(mflat,cp,Tmlat,Tflat,C1,mlat,Teb)
cp_ac=cp_solido(mfac,cp,Tmac,Tfac,C1,mac,Teb)

# %%
