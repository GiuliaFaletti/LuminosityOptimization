import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import LoadData as ld
import LuminosityOptimization as lo
import scipy.integrate as integrate
from lmfit import Model
import sympy as sp
from scipy.optimize import minimize, LinearConstraint
from scipy.integrate import quad

from DataModel import *



#2018
#initial guess
n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps= lo.Parameters2018() 
t_opt, L_tot_opt, L_int_opt, t_a18, t_opt18, L_int18, L_tot18=lo.L_optimal_18(data_ta18_sec) 
t_opt18 = lo.t_opt_eval(N_i, n_c, Xi, S_int, t_a18)

#x018=t_opt18*np.ones(len(a_18), dtype=np.float128)
#x018=np.linspace(np.amin(data_tf18_sec), np.amax(data_tf18_sec), len(data_tf18_sec))
x018=data_tf18_sec

#define the function to be optimize 
a_18=np.array(a_18, dtype=np.float128)
b_18=np.array(b_18, dtype=np.float128)
c_18=np.array(c_18, dtype=np.float128)
d_18=np.array(d_18, dtype=np.float128)

T018=np.array(T018)

print(T018)

tot=sum(data_tf18_sec)
print(tot)

def f_18(t3):
    result=np.empty(len(a_18))
    for i in range(len(a_18)):
        lam=lambda x1: a_18[i]*np.exp(-b_18[i]*x1)+c_18[i]*np.exp(-d_18[i]*x1)
        result[i]=quad(lam, T018[i], T018[i]+t3[i])[0]
    
    result = np.sum(result)
    return result

def cons(t3):
    res = np.sum(t3) - (tot+(3600*3))
    return res
    

res18 = minimize(f_18, x018, options={'disp': True}, constraints={'type':'ineq', 'fun': cons}) 
#res18 = minimize(f_18, x018, options={'disp': True}, method='BFGS')
print(res18)
print(x018, res18.x)        
  
opt_data_sum=np.sum(res18.x)/(3600*24) 
real_data_sum= tot/(3600*24)
print(opt_data_sum, real_data_sum)

w18=np.empty(len(x018))
for i in range(len(x018)):
    w18[i]=res18.x[i]-x018[i]

print(w18, len(w18))
