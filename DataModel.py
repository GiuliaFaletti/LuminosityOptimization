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

#Evaluating the double exponential data model considering the correct data sample

#loading data
data_16, data_17, data_18, array16, array17, array18 = ld.Data()
data_tot, dataTot, arrayTot = ld.TotalDataSet(data_16, data_17, data_18)
data_ta16, data_tf16, data_ta17, data_tf17, data_ta18, data_tf18 = ld.loadFill()     
FillNumber16, FillNumber17, FillNumber18 = ld.FillNumber()
data_16_sec, data_ta16_sec, data_tf16_sec, data_17_sec, data_ta17_sec,\
   data_tf17_sec, data_18_sec, data_ta18_sec, data_tf18_sec=ld.Data_sec(array16,\
      data_ta16, data_tf16, array17, data_ta17, data_tf17, array18, data_ta18, data_tf18)
   
   
   
#defining the fit function
def fit(x, a, b, c, d):
    return (a*np.exp((-b)*x))+(c*np.exp((-d)*x))

model=Model(fit)


#2016
L_int_opt2016=[]
L_intfit16=[]
a_16=[]
b_16=[]
c_16=[]
d_16=[]
T016=[]
for i in range(len(FillNumber16)):
    #plotting results 
    plt.close("all")
    fig1,  ax1 = plt.subplots()
    print("######################################################################FILL",int(FillNumber16[i]),"#########################################################")
    text = str(int(FillNumber16[i])) #number of current fill
    #obtain the Times and luminosity evolution values for that fill
    f=open('ATLAS/ATLAS_fill_2016/{}_lumi_ATLAS.txt'.format(text),"r")
    lines=f.readlines()
    L_evolx=[]
    times=[]
    for x in lines:
        times.append(int(x.split(' ')[0]))  
        L_evolx.append(float(x.split(' ')[2]))
          
    f.close()
    Times = np.array(times)
    L_evol = np.array(L_evolx)
    
    T016.append(Times[0])
    
    #deleting the null values of the luminosity
    zero=np.where(L_evol<100)
    L_zero=np.delete(L_evol, zero)
    T_zero=np.delete(Times, zero)
        
    #check for enough points
    if len(L_zero)<10:
        zero=np.where(L_evol<5)
        L_zero=np.delete(L_evol, zero)
        T_zero=np.delete(Times, zero)

    #defining the derivative 
    dy = np.zeros(L_zero.shape)
    dy[0:-1] = np.diff(L_zero)/np.diff(T_zero)

     #start to slim down the fit interval       
    L_tofit=[]
    T_tofit=[]
    for idx in range(len(L_zero)):
        #cancelling too strong derivative points
        if dy[idx]<0 and dy[idx]>-1.5:
            L_tofit.append(L_zero[idx])
            T_tofit.append(T_zero[idx])
        if dy[idx]>0 or dy[idx]<-1.5:
            continue     
        
    #evaluating the differences between two subsequent points
    diff=np.diff(L_tofit)
        
    #deleting the discrepancies
    thr=np.max(abs(diff))*0.05
    idx_diff= np.where(abs(diff)>thr)[0]+1
        
    #new slim down of data
    L_tofit2=np.delete(L_tofit, idx_diff)
    T_tofit2=np.delete(T_tofit, idx_diff)
        
    #check for enough points
    if len(L_tofit2) < 30:
        L_tofit2=L_tofit
        T_tofit2=T_tofit
        
    L_fit=L_tofit2
    T_fit=T_tofit2            
    #normalization of the fit interval    
    norm_T_fit=[]
    norm_T_fit=np.array(norm_T_fit)
    for element in T_fit:
        z=(element-np.amin(T_fit))/(np.amax(T_fit)-np.amin(T_fit))
        norm_T_fit=np.append(norm_T_fit, z)
         
    #performing fit of last segments of data
    model.set_param_hint('b', value=0.2, min=0, max=100)
    model.set_param_hint('d', value=0.2, min=0, max=100)
    fit_result=model.fit(L_fit, x=norm_T_fit, a=1, b=0.2, c=1, d=0.2)
    print(fit_result.params['a'].value, fit_result.params['b'].value, fit_result.params['c'].value, fit_result.params['d'].value)
    ax1.plot(T_fit, L_fit, "b.", label='Smoothed data', markersize=4)
    ax1.plot(T_fit, fit_result.best_fit, 'r-', label='Double exponential fit')
    ax1.plot([], [], 'kx ', label='Reduced Chi-Square ={:.5f}'.format(fit_result.redchi))
    ax1.set_xlabel('Times [s]')
    ax1.set_ylabel('Luminosity evolution [Hz/\u03BCb]')
    plt.legend(loc='best')
    L_i=integrate.simps(fit_result.best_fit, T_fit)
    L_intfit16.append(L_i)
    a_16.append(fit_result.params['a'].value)
    b_16.append(fit_result.params['b'].value)
    c_16.append(fit_result.params['c'].value)
    d_16.append(fit_result.params['d'].value)
    ax1.set_title('{}'.format(text)) 
    plt.savefig('FitModel/{}_fitModel.pdf'.format(text)) 
       

plt.close('all')
fig3, ax3=plt.subplots()
ax3.hist(a_16, bins=11, density=True)
ax3.set_title('amp_1 double exp 2016')
ax3.set_ylabel('Normalized Frequencies')
ax3.set_xlabel('Parameter Values')
plt.savefig('FitModel/a_2016.pdf') 
#plt.show()
plt.close()
fig3, ax3=plt.subplots()
ax3.hist(b_16, bins=13, density=True)
ax3.set_title('lambda_1 double exp 2016')
ax3.set_ylabel('Normalized Frequencies')
ax3.set_xlabel('Parameter Values')
plt.savefig('FitModel/b_2016.pdf') 
#plt.show()
plt.close()
fig3, ax3=plt.subplots()
ax3.hist(c_16, bins=13, density=True)
ax3.set_title('amp_2 double exp 2016')
ax3.set_ylabel('Normalized Frequencies')
ax3.set_xlabel('Parameter Values')
plt.savefig('FitModel/c_2016.pdf')
#plt.show()
plt.close()
fig3, ax3=plt.subplots()
ax3.hist(d_16, bins=15, density=True)
ax3.set_title('lambda_2 double exp 2016')
ax3.set_ylabel('Normalized Frequencies')
ax3.set_xlabel('Parameter Values')
plt.savefig('FitModel/d_2016.pdf') 
#plt.show()
plt.close()


#2017
a_17=[]
b_17=[]
c_17=[]
d_17=[]
L_intfit17=[]
L_int_opt2017=[]
T017=[]
for i in range(len(FillNumber17)):
    #plotting results 
    plt.close("all")
    fig1,  ax1 = plt.subplots()
    print("######################################################################FILL",int(FillNumber17[i]),"#########################################################")
    text = str(int(FillNumber17[i])) #number of current fill
    #obtain the Times and luminosity evolution values for that fill
    f=open('ATLAS/ATLAS_fill_2017/{}_lumi_ATLAS.txt'.format(text),"r")
    lines=f.readlines()
    L_evolx=[]
    times=[]
    for x in lines:
        times.append(int(x.split(' ')[0]))  
        L_evolx.append(float(x.split(' ')[2]))
          
    f.close()
    Times = np.array(times)
    L_evol = np.array(L_evolx)
    
    T017.append(Times[0])
    
    #deleting the null values of the luminosity
    zero=np.where(L_evol<100)
    L_zero=np.delete(L_evol, zero)
    T_zero=np.delete(Times, zero)
        
    #check for enough points
    if len(L_zero)<10:
        zero=np.where(L_evol<5)
        L_zero=np.delete(L_evol, zero)
        T_zero=np.delete(Times, zero)

    #defining the derivative 
    dy = np.zeros(L_zero.shape)
    dy[0:-1] = np.diff(L_zero)/np.diff(T_zero)

    
     #start to slim down the fit interval       
    L_tofit=[]
    T_tofit=[]
    for idx in range(len(L_zero)):
        #cancelling too strong derivative points
        if dy[idx]<0 and dy[idx]>-1.5:
            L_tofit.append(L_zero[idx])
            T_tofit.append(T_zero[idx])
        if dy[idx]>0 or dy[idx]<-1.5:
            continue     
        
    #evaluating the differences between two subsequent points
    diff=np.diff(L_tofit)
        
    #deleting the discrepancies
    thr=np.max(abs(diff))*0.05
    idx_diff= np.where(abs(diff)>thr)[0]+1
        
    #new slim down of data
    L_tofit2=np.delete(L_tofit, idx_diff)
    T_tofit2=np.delete(T_tofit, idx_diff)
        
    #check for enough points
    if len(L_tofit2) < 30:
        L_tofit2=L_tofit
        T_tofit2=T_tofit
        
    L_fit=L_tofit2
    T_fit=T_tofit2     
    
 
      
    #normalization of the fit interval    
    norm_T_fit=[]
    norm_T_fit=np.array(norm_T_fit)
    for element in T_fit:
        z=(element-np.amin(T_fit))/(np.amax(T_fit)-np.amin(T_fit))
        norm_T_fit=np.append(norm_T_fit, z)
         
    #performing fit of last segments of data
    model.set_param_hint('b', value=0.2, min=0, max=100)
    model.set_param_hint('d', value=0.2, min=0, max=100)
    fit_result=model.fit(L_fit, x=norm_T_fit, a=1, b=0.2, c=1, d=0.2)
    
    a_17.append(fit_result.params['a'].value)
    b_17.append(fit_result.params['b'].value)
    c_17.append(fit_result.params['c'].value)
    d_17.append(fit_result.params['d'].value)
    L_i=integrate.simps(fit_result.best_fit, T_fit)
    L_intfit17.append(L_i)
    ax1.plot(T_fit, L_fit, "b.", label='Smoothed data', markersize=4)
    ax1.plot(T_fit, fit_result.best_fit, 'r-', label='Best fit')
    ax1.plot([], [], 'kx ', label='Reduced Chi-Square ={:.5f}'.format(fit_result.redchi))
    ax1.set_xlabel('Times [s]')
    ax1.set_ylabel('Luminosity evolution [Hz/\u03BCb]')
    plt.legend(loc='best')
        
    ax1.set_title('{}'.format(text)) 
    plt.savefig('FitModel/{}_fitModel.pdf'.format(text)) 
    ##plt.show()     

plt.close('all')
fig3, ax3=plt.subplots()
n3, bins3, patches3 = ax3.hist(a_17, bins=20, density=True)
ax3.set_title('amp_1 double exp 2017')
ax3.set_ylabel('Normalized Frequencies')
ax3.set_xlabel('Parameter Values')
plt.savefig('FitModel/a_2017.pdf')
#plt.show()
plt.close()
fig3, ax3=plt.subplots()
ax3.hist(b_17, bins=18, density=True)
ax3.set_title('lambda_1 double exp 2017')
ax3.set_ylabel('Normalized Frequencies')
ax3.set_xlabel('Parameter Values')
plt.savefig('FitModel/b_2017.pdf')
#plt.show()
plt.close()
fig3, ax3=plt.subplots()
ax3.hist(c_17, bins=20, density=True)
ax3.set_title('amp_2 double exp 2017')
plt.savefig('FitModel/c_2017.pdf') 
ax3.set_ylabel('Normalized Frequencies')
ax3.set_xlabel('Parameter Values')
#plt.show()
plt.close()
fig3, ax3=plt.subplots()
ax3.hist(d_17, bins=15, density=True)
ax3.set_title('lambda_2 double exp 2017')
ax3.set_ylabel('Normalized Frequencies')
ax3.set_xlabel('Parameter Values')
plt.savefig('FitModel/d_2017.pdf') 
#plt.show()
plt.close()

#2018
a_18=[]
b_18=[]
c_18=[]
d_18=[]
L_intfit18=[]
L_int_opt2018=[]
T018=[]
for i in range(len(FillNumber18)):
    #plotting results 
    plt.close("all")
    fig1,  ax1 = plt.subplots()
    print("######################################################################FILL",int(FillNumber18[i]),"#########################################################")
    text = str(int(FillNumber18[i])) #number of current fill
    #obtain the Times and luminosity evolution values for that fill
    f=open('ATLAS/ATLAS_fill_2018/{}_lumi_ATLAS.txt'.format(text),"r")
    lines=f.readlines()
    L_evolx=[]
    times=[]
    for x in lines:
        times.append(int(x.split(' ')[0]))  
        L_evolx.append(float(x.split(' ')[2]))
          
    f.close()
    Times = np.array(times)
    L_evol = np.array(L_evolx)
    
    T018.append(Times[0])
    
    #deleting the null values of the luminosity
    zero=np.where(L_evol<100)
    L_zero=np.delete(L_evol, zero)
    T_zero=np.delete(Times, zero)
        
    #check for enough points
    if len(L_zero)<10:
        zero=np.where(L_evol<5)
        L_zero=np.delete(L_evol, zero)
        T_zero=np.delete(Times, zero)

    #defining the derivative 
    dy = np.zeros(L_zero.shape)
    dy[0:-1] = np.diff(L_zero)/np.diff(T_zero)

 
     #start to slim down the fit interval       
    L_tofit=[]
    T_tofit=[]
    for idx in range(len(L_zero)):
        #cancelling too strong derivative points
        if dy[idx]<0 and dy[idx]>-1.5:
            L_tofit.append(L_zero[idx])
            T_tofit.append(T_zero[idx])
        if dy[idx]>0 or dy[idx]<-1.5:
            continue     
        
    #evaluating the differences between two subsequent points
    diff=np.diff(L_tofit)
        
    #deleting the discrepancies
    thr=np.max(abs(diff))*0.05
    idx_diff= np.where(abs(diff)>thr)[0]+1
        
    #new slim down of data
    L_tofit2=np.delete(L_tofit, idx_diff)
    T_tofit2=np.delete(T_tofit, idx_diff)
        
    #check for enough points
    if len(L_tofit2) < 30:
        L_tofit2=L_tofit
        T_tofit2=T_tofit
        
    L_fit=L_tofit2
    T_fit=T_tofit2 
         
    #normalization of the fit interval    
    norm_T_fit=[]
    norm_T_fit=np.array(norm_T_fit)
    for element in T_fit:
        z=(element-np.amin(T_fit))/(np.amax(T_fit)-np.amin(T_fit))
        norm_T_fit=np.append(norm_T_fit, z)
         
    #performing fit of last segments of data
    model.set_param_hint('b', value=0.2, min=0, max=100)
    model.set_param_hint('d', value=0.2, min=0, max=100)
    fit_result=model.fit(L_fit, x=norm_T_fit, a=1, b=0.2, c=1, d=0.2)
    
    ax1.plot(T_fit, L_fit, "b.", label='Smoothed data', markersize=4)
    ax1.plot(T_fit, fit_result.best_fit, 'r-', label='Double Exponential fit')
    ax1.set_xlabel('Times [s]')
    ax1.set_ylabel('Luminosity evolution [Hz/\u03BCb]')
    ax1.plot([], [], 'kx ', label='Reduced Chi-Square ={:.5f}'.format(fit_result.redchi))
    plt.legend(loc='best')
    L_i=integrate.simps(fit_result.best_fit, T_fit)
    L_intfit18.append(L_i)
          
    a_18.append(fit_result.params['a'].value)
    b_18.append(fit_result.params['b'].value)
    c_18.append(fit_result.params['c'].value)
    d_18.append(fit_result.params['d'].value)
            
    ax1.set_title('{}'.format(text)) 
    plt.savefig('FitModel/{}_fitModel.pdf'.format(text)) 
    ##plt.show()     

plt.close('all')
fig3, ax3=plt.subplots()
ax3.hist(a_18, bins=25, density=True)
ax3.set_title('amp_1 double exp 2018')
ax3.set_ylabel('Normalized Frequencies')
ax3.set_xlabel('Parameter Values')
plt.savefig('FitModel/a_2018.pdf') 
#plt.show()
plt.close()
fig3, ax3=plt.subplots()
ax3.hist(b_18, bins=25, density=True)
ax3.set_title('lambda_1 double exp 2018')
ax3.set_ylabel('Normalized Frequencies')
ax3.set_xlabel('Parameter Values')
plt.savefig('FitModel/b_2018.pdf') 
#plt.show()
plt.close()
fig3, ax3=plt.subplots()
ax3.hist(c_18, bins=23, density=True)
ax3.set_title('amp_2 double exp 2018')
ax3.set_ylabel('Normalized Frequencies')
ax3.set_xlabel('Parameter Values')
plt.savefig('FitModel/c_2018.pdf') 
##plt.show()
plt.close()
fig3, ax3=plt.subplots()
ax3.hist(d_18, bins=25, density=True)
ax3.set_title('lambda_2 double exp 2018')
ax3.set_ylabel('Normalized Frequencies')
ax3.set_xlabel('Parameter Values')
plt.savefig('FitModel/d_2018.pdf') 
##plt.show()
plt.close()


#Correlation between parameters
corr1=np.corrcoef(a_16, b_16)
corr2=np.corrcoef(c_16, d_16)
corr3=np.corrcoef(a_16, c_16)
corr4=np.corrcoef(a_16, d_16)
corr5=np.corrcoef(b_16, d_16)
corr6=np.corrcoef(c_16, b_16)

print(corr1[0,1], corr2[0,1], corr3[0,1], corr4[0,1], corr5[0,1], corr6[0,1])

corr1=np.corrcoef(a_17, b_17)
corr2=np.corrcoef(c_17, d_17)
corr3=np.corrcoef(a_17, c_17)
corr4=np.corrcoef(a_17, d_17)
corr5=np.corrcoef(b_17, d_17)
corr6=np.corrcoef(c_17, b_17)

print(corr1[0,1], corr2[0,1], corr3[0,1], corr4[0,1], corr5[0,1], corr6[0,1])

corr1=np.corrcoef(a_18, b_18)
corr2=np.corrcoef(c_18, d_18)
corr3=np.corrcoef(a_18, c_18)
corr4=np.corrcoef(a_18, d_18)
corr5=np.corrcoef(b_18, d_18)
corr6=np.corrcoef(c_18, b_18)

print(corr1[0,1], corr2[0,1], corr3[0,1], corr4[0,1], corr5[0,1], corr6[0,1])

plt.close('all')
fig4, ax4=plt.subplots()
ax4.plot(a_16, b_16, "b.")
ax4.set_title('a16/b16')
ax4.set_ylabel('b_16')
ax4.set_xlabel('a_16')
plt.savefig('FitModel/a16_b16.pdf') 
plt.close()
fig4, ax4=plt.subplots()
ax4.plot(a_16, d_16, "b.")
ax4.set_title('a16/d16')
ax4.set_ylabel('d_16')
ax4.set_xlabel('a_16')
plt.savefig('FitModel/a16_d16.pdf') 
plt.close()
fig4, ax4=plt.subplots()
ax4.plot(b_16, d_16, "b.")
ax4.set_title('b16/d16')
ax4.set_ylabel('d_16')
ax4.set_xlabel('b_16')
plt.savefig('FitModel/b16_d16.pdf') 
plt.close()
fig4, ax4=plt.subplots()
ax4.plot(c_16, b_16, "b.")
ax4.set_title('c16/b16')
ax4.set_ylabel('b_16')
ax4.set_xlabel('c_16')
plt.savefig('FitModel/c16_b16.pdf') 
plt.close()
fig4, ax4=plt.subplots()
ax4.plot(a_16, c_16, "b.")
ax4.set_title('a16/c16')
ax4.set_ylabel('c_16')
ax4.set_xlabel('a_16')
plt.savefig('FitModel/a16_c16.pdf') 
plt.close()
fig4, ax4=plt.subplots()
ax4.plot(c_16, d_16, "b.")
ax4.set_title('c16/d16')
ax4.set_ylabel('d_16')
ax4.set_xlabel('c_16')
plt.savefig('FitModel/c16_d16.pdf') 
plt.close()
fig4, ax4=plt.subplots()
ax4.plot(a_17, b_17, "b.")
ax4.set_title('a17/b17')
ax4.set_ylabel('d_17')
ax4.set_xlabel('a_17')
plt.savefig('FitModel/a17_b17.pdf') 
plt.close()
fig4, ax4=plt.subplots()
ax4.plot(a_17, d_17, "b.")
ax4.set_title('a17/d17')
ax4.set_ylabel('d_17')
ax4.set_xlabel('a_17')
plt.savefig('FitModel/a17_d17.pdf') 
plt.close()
fig4, ax4=plt.subplots()
ax4.plot(b_17, d_17, "b.")
ax4.set_title('b17/d17')
ax4.set_ylabel('d_17')
ax4.set_xlabel('b_17')
plt.savefig('FitModel/b17_d17.pdf') 
plt.close()
fig4, ax4=plt.subplots()
ax4.plot(a_17, c_17, "b.")
ax4.set_title('a17/c17')
ax4.set_ylabel('c_17')
ax4.set_xlabel('a_17')
plt.savefig('FitModel/a17_c17.pdf') 
plt.close()
fig4, ax4=plt.subplots()
ax4.plot(c_17, d_17, "b.")
ax4.set_title('c17/d17')
ax4.set_ylabel('d_17')
ax4.set_xlabel('c_17')
plt.savefig('FitModel/c17_d17.pdf') 
plt.close()
fig4, ax4=plt.subplots()
ax4.plot(b_17, c_17, "b.")
ax4.set_title('b17/c17')
ax4.set_ylabel('c_17')
ax4.set_xlabel('b_17')
plt.savefig('FitModel/b17_c17.pdf') 
plt.close()
fig4, ax4=plt.subplots()
ax4.plot(a_18, b_18, "b.")
ax4.set_title('a18/b18')
ax4.set_ylabel('b_18')
ax4.set_xlabel('a_18')
plt.savefig('FitModel/a18_b18.pdf') 
plt.close()
fig4, ax4=plt.subplots()
ax4.plot(a_18, d_18, "b.")
ax4.set_title('a18/d18')
ax4.set_ylabel('d_18')
ax4.set_xlabel('a_18')
plt.savefig('FitModel/a18_d18.pdf') 
plt.close()
fig4, ax4=plt.subplots()
ax4.plot(b_18, d_18, "b.")
ax4.set_title('b18/d18')
ax4.set_ylabel('d_18')
ax4.set_xlabel('b_18')
plt.savefig('FitModel/b18_d18.pdf') 
plt.close()
fig4, ax4=plt.subplots()
ax4.plot(a_18, c_18, "b.")
ax4.set_title('a18/c18')
ax4.set_ylabel('c_18')
ax4.set_xlabel('a_18')
plt.savefig('FitModel/a18_c18.pdf') 
plt.close()
fig4, ax4=plt.subplots()
ax4.plot(c_18, d_18, "b.")
ax4.set_title('c18/d18')
ax4.set_ylabel('d_18')
ax4.set_xlabel('c_18')
plt.savefig('FitModel/c18_d18.pdf') 
plt.close()
fig4, ax4=plt.subplots()
ax4.plot(b_18, c_18, "b.")
ax4.set_title('b18/c18')
ax4.set_ylabel('c_18')
ax4.set_xlabel('b_18')
plt.savefig('FitModel/b18_c18.pdf') 
plt.close()


