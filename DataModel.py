import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import LoadData as ld
import LuminosityOptimization as lo
import scipy.integrate as integrate
from lmfit import Model
import sympy as sp
#from Opt_MeasuredLumi_Extrapolation_average import L_int_opt

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

#A,B,C,D, tau = sp.symbols("a,b,c,d, tau", real=True)
#x,t=sp.symbols("x, t", real=True, positive=True)
#exp1=A*sp.exp(-B*x)+C*sp.exp(-D*x)-(tau+x)*((A/B)-(A/B)*sp.exp(-B*x)+(C/D)-(C/D)*sp.exp(-D*x))
#exp2=A*sp.exp(-B*x)-(tau+x)*((A/B)-(A/B)*sp.exp(-B*x))
#solution1 = sp.solve(exp1, x)
#solution2 = sp.solve(exp2, x)
#defining the fit model
model=Model(fit)


#2016
L_int_opt2016=[]
DoubleExp=[]
SingleExp=[]
L_intfit16=[]
a_16=[]
b_16=[]
c_16=[]
d_16=[]
a1_16=[]
b1_16=[]
for i in range(len(FillNumber16)):
    #n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps= lo.Parameters2016() #importing the theoretical luminosity model parameters
    #evaluate the optimal time for the current fill
    #t_opt = lo.t_opt_eval(N_i, n_c, Xi, S_int, data_ta16_sec[i])
    #t_sec=data_tf16_sec[-1]-data_tf16_sec[0]
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
    L_fit=[]
    T_fit=[]
    for idx in range(len(L_zero)):
        #cancelling too strong derivative points
        if dy[idx]<0 and dy[idx]>-1.5:
            L_fit.append(L_zero[idx])
            T_fit.append(T_zero[idx])
        if dy[idx]>0 or dy[idx]<-1.5:
            continue     
    #normalization of the fit interval    
    norm_T_fit=[]
    norm_T_fit=np.array(norm_T_fit)
    for element in T_fit:
        z=(element-np.amin(T_fit))/(np.amax(T_fit)-np.amin(T_fit))
        norm_T_fit=np.append(norm_T_fit, z)
         
    #performing fit of last segments of data
    fit_result=model.fit(L_fit, x=norm_T_fit, a=1, b=0.2, c=1, d=0.2)
    print(fit_result.params['a'].value, fit_result.params['b'].value, fit_result.params['c'].value, fit_result.params['d'].value)
    
    
    #if the double exponential grows changing the fit model --> one exponential
    if fit_result.params['b'].value <0 or fit_result.params['d'].value <0: 
        print('Single exponential!')
        def fit_lin(x, a2, b2):
            return a2*np.exp((-b2)*x)
        
        model2=Model(fit_lin)
        fit_result2=model2.fit(L_fit, x=norm_T_fit, a2=1, b2=1)
        print(fit_result2.params['a2'].value, fit_result2.params['b2'].value)
        L_i=integrate.simps(fit_result2.best_fit, T_fit)
        L_intfit16.append(L_i)
        a1_16.append(fit_result2.params['a2'].value)
        b1_16.append(fit_result2.params['b2'].value)
        ax1.plot(T_fit, L_fit, "b.", label='Fit interval', markersize=4)
        ax1.plot(T_fit, fit_result2.best_fit, 'r-' )
        ax1.plot([], [], 'kx ', label='Reduced Chi-Square ={:.5f}'.format(fit_result2.redchi))
        ax1.plot([],[], 'r>', label='Single exponential')
        plt.legend(loc='best')
        SingleExp.append(int(FillNumber16[i]))
        #t2=solution2.subs((A, fit_result2.params['a2'].value), (B, fit_result2.params['b2'].value), (x, t_sec))
        
    elif fit_result.params['b'].value >0 or fit_result.params['d'].value >0:
        ax1.plot(T_fit, L_fit, "b.", label='Fit interval', markersize=4)
        ax1.plot(T_fit, fit_result.best_fit, 'r-')
        ax1.plot([], [], 'kx ', label='Reduced Chi-Square ={:.5f}'.format(fit_result.redchi))
        ax1.plot([],[], 'r>', label='Double exponential')
        plt.legend(loc='best')
        L_i=integrate.simps(fit_result.best_fit, T_fit)
        L_intfit16.append(L_i)
        DoubleExp.append(int(FillNumber16[i]))
        a_16.append(fit_result.params['a'].value)
        b_16.append(fit_result.params['b'].value)
        c_16.append(fit_result.params['c'].value)
        d_16.append(fit_result.params['d'].value)
        #t2=solution1.subs((A,fit_result.params['a'].value ), (B, fit_result.params['b'].value), (C, fit_result.params['c'].value), (D, fit_result.params['d'].value), (x, t_sec))
    
    ax1.set_title('{}'.format(text)) 
    plt.savefig('FitModel/{}_fitModel.pdf'.format(text)) 
    ##plt.show() 
    #print('t ottimale modello - t ottimale dati')
    #print(t_opt, t2)        

print(len(FillNumber16))   
y1=[]
y2=[] 
for i in range(len(SingleExp)):
    y1.append(1)
for i in range(len(DoubleExp)):
    y2.append(2)

    
fig2, ax2=plt.subplots()
ax2.plot(SingleExp, y1, "b*", label='{}'.format(len(SingleExp)))   
ax2.plot(DoubleExp, y2, "r*", label='{}'.format(len(DoubleExp)))  
ax2.plot([], [], "k.", label='{}'.format(len(SingleExp)+len(DoubleExp)))
ax2.set_title('Models 2016')
plt.legend(loc='best')
plt.savefig('FitModel/2016_WhatModel.pdf') 
##plt.show() 

print(max(a_16), min(a_16))
print(max(b_16), min(b_16))
print(max(c_16), min(c_16))
print(max(d_16), min(d_16))
print(max(a1_16), min(a1_16))
print(max(b1_16), min(b1_16))

plt.close('all')
fig3, ax3=plt.subplots()
ax3.hist(a_16, bins=9, density=True)
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
fig3, ax3=plt.subplots()
ax3.hist(a1_16, bins=9, density=True)
ax3.set_title('amplitude single exp 2016')
ax3.set_ylabel('Normalized Frequencies')
ax3.set_xlabel('Parameter Values')
plt.savefig('FitModel/a1_2016.pdf') 
#plt.show()
plt.close()
fig3, ax3=plt.subplots()
ax3.hist(b1_16, bins=10, density=True)
ax3.set_title('lambda single exp 2016')
ax3.set_ylabel('Normalized Frequencies')
ax3.set_xlabel('Parameter Values')
plt.savefig('FitModel/b1_2016.pdf') 
#plt.show()
plt.close()

#2017
a_17=[]
b_17=[]
c_17=[]
d_17=[]
a1_17=[]
b1_17=[]
L_intfit17=[]
L_int_opt2017=[]
DoubleExp=[]
SingleExp=[]
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
    L_fit=[]
    T_fit=[]
    for idx in range(len(L_zero)):
        #cancelling too strong derivative points
        if dy[idx]<0 and dy[idx]>-1.5:
            L_fit.append(L_zero[idx])
            T_fit.append(T_zero[idx])
        if dy[idx]>0 or dy[idx]<-1.5:
            continue     
    #normalization of the fit interval    
    norm_T_fit=[]
    norm_T_fit=np.array(norm_T_fit)
    for element in T_fit:
        z=(element-np.amin(T_fit))/(np.amax(T_fit)-np.amin(T_fit))
        norm_T_fit=np.append(norm_T_fit, z)
         
    #performing fit of last segments of data
    fit_result=model.fit(L_fit, x=norm_T_fit, a=1, b=0.2, c=1, d=0.2)
    print(fit_result.params['a'].value, fit_result.params['b'].value, fit_result.params['c'].value, fit_result.params['d'].value)
    
    
    #if the double exponential grows changing the fit model --> one exponential
    if fit_result.params['b'].value <0 or fit_result.params['d'].value <0: 
        print('Single exponential!')
        def fit_lin(x, a2, b2):
            return a2*np.exp((-b2)*x)
        
        model2=Model(fit_lin)
        fit_result2=model2.fit(L_fit, x=norm_T_fit, a2=1, b2=1)
        print(fit_result2.params['a2'].value, fit_result2.params['b2'].value)
        L_i=integrate.simps(fit_result2.best_fit, T_fit)
        L_intfit17.append(L_i)
        a1_17.append(fit_result2.params['a2'].value)
        b1_17.append(fit_result2.params['b2'].value)
        ax1.plot(T_fit, L_fit, "b.", label='Fit interval', markersize=4)
        ax1.plot(T_fit, fit_result2.best_fit, 'r-' )
        ax1.plot([], [], 'kx ', label='Reduced Chi-Square ={:.5f}'.format(fit_result2.redchi))
        ax1.plot([],[], 'r>', label='Single exponential')
        plt.legend(loc='best')
        SingleExp.append(int(FillNumber17[i]))
        
    elif fit_result.params['b'].value >0 or fit_result.params['d'].value >0:
        a_17.append(fit_result.params['a'].value)
        b_17.append(fit_result.params['b'].value)
        c_17.append(fit_result.params['c'].value)
        d_17.append(fit_result.params['d'].value)
        L_i=integrate.simps(fit_result.best_fit, T_fit)
        L_intfit17.append(L_i)
        ax1.plot(T_fit, L_fit, "b.", label='Fit interval', markersize=4)
        ax1.plot(T_fit, fit_result.best_fit, 'r-')
        ax1.plot([], [], 'kx ', label='Reduced Chi-Square ={:.5f}'.format(fit_result.redchi))
        ax1.plot([],[], 'r>', label='Double exponential')
        plt.legend(loc='best')
        DoubleExp.append(int(FillNumber17[i]))
    
    ax1.set_title('{}'.format(text)) 
    plt.savefig('FitModel/{}_fitModel.pdf'.format(text)) 
    ##plt.show()     

print(len(FillNumber17))   
y1=[]
y2=[] 
for i in range(len(SingleExp)):
    y1.append(1)
for i in range(len(DoubleExp)):
    y2.append(2)

    
fig2, ax2=plt.subplots()
ax2.plot(SingleExp, y1, "b*", label='{}'.format(len(SingleExp)))   
ax2.plot(DoubleExp, y2, "r*", label='{}'.format(len(DoubleExp)))  
ax2.plot([], [], "k.", label='{}'.format(len(SingleExp)+len(DoubleExp)))
ax2.set_title('Models 2017')
plt.legend(loc='best')
plt.savefig('FitModel/2017_whatModel.pdf') 
##plt.show() 

print(max(a_17), min(a_17))
print(max(b_17), min(b_17))
print(max(c_17), min(c_17))
print(max(d_17), min(d_17))
print(max(a1_17), min(a1_17))
print(max(b1_17), min(b1_17))

plt.close('all')
fig3, ax3=plt.subplots()
n3, bins3, patches3 = ax3.hist(a_17, bins=20, density=True)
ax3.set_title('amp_1 double exp 2017')
plt.savefig('FitModel/a_2017.pdf') 
ax3.set_ylabel('Normalized Frequencies')
ax3.set_xlabel('Parameter Values')
#plt.show()
plt.close()
fig3, ax3=plt.subplots()
ax3.hist(b_17, bins=18, density=True)
ax3.set_title('lambda_1 double exp 2017')
plt.savefig('FitModel/b_2017.pdf') 
ax3.set_ylabel('Normalized Frequencies')
ax3.set_xlabel('Parameter Values')
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
plt.savefig('FitModel/d_2017.pdf') 
ax3.set_ylabel('Normalized Frequencies')
ax3.set_xlabel('Parameter Values')
#plt.show()
plt.close()
fig3, ax3=plt.subplots()
ax3.hist(a1_17, bins=10, density=True)
ax3.set_title('amplitude single exp 2017')
plt.savefig('FitModel/a1_2017.pdf') 
ax3.set_ylabel('Normalized Frequencies')
ax3.set_xlabel('Parameter Values')
#plt.show()
plt.close()
fig3, ax3=plt.subplots()
ax3.hist(b1_17, bins=10, density=True)
ax3.set_title('lambda single exp 2017')
plt.savefig('FitModel/b1_2017.pdf') 
ax3.set_ylabel('Normalized Frequencies')
ax3.set_xlabel('Parameter Values')
#plt.show()
plt.close()

#2018
a_18=[]
b_18=[]
c_18=[]
d_18=[]
L_intfit18=[]
a1_18=[]
b1_18=[]
L_int_opt2018=[]
DoubleExp=[]
SingleExp=[]
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
    L_fit=[]
    T_fit=[]
    for idx in range(len(L_zero)):
        #cancelling too strong derivative points
        if dy[idx]<0 and dy[idx]>-1:
            L_fit.append(L_zero[idx])
            T_fit.append(T_zero[idx])
        if dy[idx]>0 or dy[idx]<-1:
            continue     
    #normalization of the fit interval    
    norm_T_fit=[]
    norm_T_fit=np.array(norm_T_fit)
    for element in T_fit:
        z=(element-np.amin(T_fit))/(np.amax(T_fit)-np.amin(T_fit))
        norm_T_fit=np.append(norm_T_fit, z)
         
    #performing fit of last segments of data
    fit_result=model.fit(L_fit, x=norm_T_fit, a=1, b=0.2, c=1, d=0.2)
    print(fit_result.params['a'].value, fit_result.params['b'].value, fit_result.params['c'].value, fit_result.params['d'].value)
    
    
    #if the double exponential grows changing the fit model --> one exponential

    if fit_result.params['b'].value <0 or fit_result.params['d'].value <0: 
        print('Single exponential!')
        def fit_lin(x, a2, b2):
            return a2*np.exp((-b2)*x)
        
        model2=Model(fit_lin)
        fit_result2=model2.fit(L_fit, x=norm_T_fit, a2=1, b2=1)
        print(fit_result2.params['a2'].value, fit_result2.params['b2'].value)
        L_i=integrate.simps(fit_result2.best_fit, T_fit)
        L_intfit18.append(L_i)
        a1_18.append(fit_result2.params['a2'].value)
        b1_18.append(fit_result2.params['b2'].value)
        ax1.plot(T_fit, L_fit, "b.", label='Fit interval', markersize=4)
        ax1.plot(T_fit, fit_result2.best_fit, 'r-' )
        ax1.plot([], [], 'kx ', label='Reduced Chi-Square ={:.5f}'.format(fit_result2.redchi))
        ax1.plot([],[], 'r>', label='Single exponential')
        plt.legend(loc='best') 
        SingleExp.append(int(FillNumber18[i]))
        if fit_result.params['b'].value >0 or fit_result.params['d'].value >0:
            ax1.plot(T_fit, L_fit, "b.", label='Fit interval', markersize=4)
            ax1.plot(T_fit, fit_result.best_fit, 'r-')
            ax1.plot([], [], 'kx ', label='Reduced Chi-Square ={:.5f}'.format(fit_result.redchi))
            ax1.plot([],[], 'r>', label='Double exponential')
            plt.legend(loc='best')
            L_i=integrate.simps(fit_result.best_fit, T_fit)
            L_intfit18.append(L_i)
            DoubleExp.append(int(FillNumber18[i]))
            a_18.append(fit_result.params['a'].value)
            b_18.append(fit_result.params['b'].value)
            c_18.append(fit_result.params['c'].value)
            d_18.append(fit_result.params['d'].value)
            
            ax1.set_title('{}'.format(text)) 
            plt.savefig('FitModel/{}_fitModel.pdf'.format(text)) 
            ##plt.show()     

print(len(FillNumber18))   
y1=[]
y2=[] 
for i in range(len(SingleExp)):
    y1.append(1)
for i in range(len(DoubleExp)):
    y2.append(2)

    
fig2, ax2=plt.subplots()
ax2.plot(SingleExp, y1, "b*", label='{}'.format(len(SingleExp)))   
ax2.plot(DoubleExp, y2, "r*", label='{}'.format(len(DoubleExp)))  
ax2.plot([], [], "k.", label='{}'.format(len(SingleExp)+len(DoubleExp)))
ax2.set_title('Models 2018')
plt.legend(loc='best')
plt.savefig('FitModel/2018_WhatModel.pdf') 
##plt.show()

print(max(a_18), min(a_18))
print(max(b_18), min(b_18))
print(max(c_18), min(c_18))
print(max(d_18), min(d_18))
print(max(a1_18), min(a1_18))
print(max(b1_18), min(b1_18))

plt.close('all')
fig3, ax3=plt.subplots()
ax3.hist(a_18, bins=25, density=True)
ax3.set_title('amp_1 double exp 2018')
plt.savefig('FitModel/a_2018.pdf') 
ax3.set_ylabel('Normalized Frequencies')
ax3.set_xlabel('Parameter Values')
#plt.show()
plt.close()
fig3, ax3=plt.subplots()
ax3.hist(b_18, bins=25, density=True)
ax3.set_title('lambda_1 double exp 2018')
plt.savefig('FitModel/b_2018.pdf') 
ax3.set_ylabel('Normalized Frequencies')
ax3.set_xlabel('Parameter Values')
#plt.show()
plt.close()
fig3, ax3=plt.subplots()
ax3.hist(c_18, bins=23, density=True)
ax3.set_title('amp_2 double exp 2018')
plt.savefig('FitModel/c_2018.pdf') 
ax3.set_ylabel('Normalized Frequencies')
ax3.set_xlabel('Parameter Values')
##plt.show()
plt.close()
fig3, ax3=plt.subplots()
ax3.hist(d_18, bins=25, density=True)
ax3.set_title('lambda_2 double exp 2018')
plt.savefig('FitModel/d_2018.pdf') 
ax3.set_ylabel('Normalized Frequencies')
ax3.set_xlabel('Parameter Values')
##plt.show()
plt.close()

fig3, ax3=plt.subplots()
ax3.hist(a1_18, bins=25, density=True)
ax3.set_title('amplitude single exp 2018')
plt.savefig('FitModel/a1_2018.pdf') 
ax3.set_ylabel('Normalized Frequencies')
ax3.set_xlabel('Parameter Values')
#plt.show()
plt.close()
fig3, ax3=plt.subplots()
ax3.hist(b1_18, bins=23, density=True)
ax3.set_title('lambda single exp 2018')
plt.savefig('FitModel/b1_2018.pdf') 
ax3.set_ylabel('Normalized Frequencies')
ax3.set_xlabel('Parameter Values')
#plt.show()
plt.close()


