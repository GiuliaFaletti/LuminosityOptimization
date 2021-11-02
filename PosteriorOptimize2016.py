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
#print(f_16)   
#constraint = LinearConstraint(np.ones(z), lb=tot, ub=tot)
#res16 = minimize(func, x0,  constraints=constraint, method='BFGS') #method='BFGS', options={'disp': True},
#print(x0, res16.x)


#2016
#initial guess
n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps= lo.Parameters2016() 
t_opt, L_tot_opt, L_int_opt, t_a16, t_opt16, L_int16, L_tot16=lo.L_optimal_16(data_ta16_sec)  
t_opt16 = lo.t_opt_eval(N_i, n_c, Xi, S_int, t_a16)
#x016=t_opt16*np.ones(len(a_16), dtype=np.float128)
#x016=np.linspace(np.amin(data_tf16_sec), np.amax(data_tf16_sec), len(data_tf16_sec))
x016=data_tf16_sec

#define the function to be optimize 
a_16=np.array(a_16, dtype=np.float128)
b_16=np.array(b_16, dtype=np.float128)
c_16=np.array(c_16, dtype=np.float128)
d_16=np.array(d_16, dtype=np.float128)

T016=np.array(T016)

tot=sum(data_tf16_sec)
print(tot)

def f_16(t1):
    result=np.empty(len(a_16))
    for i in range(len(a_16)):
        lam=lambda x1: a_16[i]*np.exp(-b_16[i]*x1)+c_16[i]*np.exp(-d_16[i]*x1)
        result[i]=quad(lam, T016[i], T016[i]+t1[i])[0]
    
    result = np.sum(result)
    return result

def cons(t1):
    res = np.sum(t1) - (tot+3600)
    return res
    

res16 = minimize(f_16, x016, options={'disp': True}, constraints={'type':'ineq', 'fun': cons}) 
#res16 = minimize(f_16, x0, options={'disp': True}, method='BFGS')
print(res16)
print(x016, res16.x)        
  
opt_data_sum=np.sum(res16.x)/(3600*24) 
real_data_sum= tot/(3600*24)
print(opt_data_sum, real_data_sum)

w16=np.empty(len(x016))
for i in range(len(x016)):
    w16[i]=res16.x[i]-x016[i]

print(w16, len(w16))
    


