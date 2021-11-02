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
#print(f_17)   
#constraint = LinearConstraint(np.ones(z), lb=tot, ub=tot)
#res17 = minimize(func, x0,  constraints=constraint, method='BFGS') #method='BFGS', options={'disp': True},
#print(x0, res17.x)

#2017
#initial guess
n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps= lo.Parameters2017() 
t_opt, L_tot_opt, L_int_opt, t_a17, t_opt17, L_int17, L_tot17=lo.L_optimal_17(data_ta17_sec)  
t_opt17 = lo.t_opt_eval(N_i, n_c, Xi, S_int, t_a17)

#x017=t_opt17*np.ones(len(a_17), dtype=np.float128)
#x017=np.linspace(np.amin(data_tf17_sec), np.amax(data_tf17_sec), len(data_tf17_sec))
x017=data_tf17_sec

#define the function to be optimize 
a_17=np.array(a_17, dtype=np.float128)
b_17=np.array(b_17, dtype=np.float128)
c_17=np.array(c_17, dtype=np.float128)
d_17=np.array(d_17, dtype=np.float128)

T017=np.array(T017)


tot=sum(data_tf17_sec)
print(tot)

def f_17(t2):
    result=np.empty(len(a_17))
    for i in range(len(a_17)):
        lam=lambda x1: a_17[i]*np.exp(-b_17[i]*x1)+c_17[i]*np.exp(-d_17[i]*x1)
        result[i]=quad(lam, T017[i], T017[i]+t2[i])[0]
    
    result = np.sum(result)
    return result

def cons(t2):
    res = np.sum(t2) - (tot+3600)
    return res
    

res17 = minimize(f_17, x017, options={'disp': True}, constraints={'type':'ineq', 'fun': cons}) 
#res17 = minimize(f_17, x0, options={'disp': True}, method='BFGS')
print(res17)
print(x017, res17.x)        
  
opt_data_sum=np.sum(res17.x)/(3600*24) 
real_data_sum= tot/(3600*24)
print(opt_data_sum, real_data_sum)

w17=np.empty(len(x017))
for i in range(len(x017)):
    w17[i]=res17.x[i]-x017[i]

print(w17, len(w17))