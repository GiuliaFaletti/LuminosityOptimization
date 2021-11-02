##############################################################################################################################################
# @Giulia Faletti
# Extrapolating the luminosity evolution for smaller or longer fill times
# based on the optimization model proposed
##############################################################################################################################################
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import LoadData as ld
import LuminosityOptimization as lo
import scipy.integrate as integrate
from lmfit import Model
import time 


#defining the start time of the code
start=time.time()

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

#defining the fit model
model=Model(fit)

#2016
L_int_opt2016=[]
for i in range(len(FillNumber16)):
     print(i)
     n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps= lo.Parameters2016() #importing the theoretical luminosity model parameters
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
     
     #evaluate the optimal time for the current fill
     t_opt = lo.t_opt_eval(N_i, n_c, Xi, S_int, data_ta16_sec[i])
     print("t_opt=", t_opt,"Last time=", Times[len(Times)-1], "ta=", data_ta16_sec[i])
     
     #compare the optimal time with the last value of the Times array
     t_opt = Times[0]+t_opt
     print("Times[0]+t_opt=", t_opt)
     control = t_opt>Times[len(Times)-1]   
     print("Times[0]+t_opt=", t_opt, "Extrapolation needed:", control)
     limit=Times[len(Times)-1] #last value of the Times array
     
     #plotting results 
     plt.close("all")
     fig,  ax = plt.subplots()
     
     #performing extrapolation when needed
     if int(t_opt)>limit:
         
        X = np.arange(limit+60, int(t_opt), 60) #definition of the extrapolation interval
        
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
            
        #defining the fit intervals
        j=0
        L_fit=np.array([])
        T_fit=np.array([])
        L_tofit2=np.array(L_tofit2)
        T_tofit2=np.array(T_tofit2)
        dy2 = np.zeros(L_tofit2.shape)
        dy2[0:-1] = np.diff(L_tofit2)/np.diff(T_tofit2)
        #filling the fit interval until it encounters a big step in the data given by the derivatives
        for a in range(len(L_tofit2)):
            if abs(dy2[len(dy2)-1-j])<0.6:
                L_fit=np.append(L_fit, L_tofit2[len(L_tofit2)-1-j])
                T_fit=np.append(T_fit, T_tofit2[len(L_tofit2)-1-j])
                j=j+1
            elif abs(dy2[len(dy2)-1-j])>1.5: 
                break         
           

        #normalization of the fit interval    
        norm_T_fit=[]
        norm_T_fit=np.array(norm_T_fit)
        norm_X=[]
        norm_X=np.array(norm_X)
        for element in T_fit:
            z=(element-np.amin(T_fit))/(np.amax(X)-np.amin(T_fit))
            norm_T_fit=np.append(norm_T_fit, z)
    
        for element in X:
            z=(element-np.amin(T_fit))/(np.amax(X)-np.amin(T_fit))
            norm_X=np.append(norm_X, z)   
             
        
        #performing fit of last segments of data
        fit_result=model.fit(L_fit, x=norm_T_fit, a=1, b=0.2, c=1, d=0.2)
        Y = fit(norm_X, fit_result.params['a'].value, fit_result.params['b'].value, fit_result.params['c'].value, fit_result.params['d'].value)
        print(fit_result.params['a'].value, fit_result.params['b'].value, fit_result.params['c'].value, fit_result.params['d'].value)
        #if the double exponential grows changing the fit model --> one exponential
        if fit_result.params['b'].value <0 or fit_result.params['d'].value <0: 
            print('Single exponential!')
            def fit_lin(x, a2, b2):
                return a2*np.exp((-b2)*x)

            model2=Model(fit_lin)
            fit_result2=model2.fit(L_fit, x=norm_T_fit, a2=1, b2=1)
            print(fit_result2.params['a2'].value, fit_result2.params['b2'].value)
            Y = fit_lin(norm_X, fit_result2.params['a2'].value, fit_result2.params['b2'].value)
        
        k1=0   
        L=np.append(L_zero, Y)
        T=np.append(T_zero, X)
        Lumi_evol=np.array([])
        Time=np.array([])
        for m in T:
            if int(t_opt)>=m:
                Lumi_evol=np.append(Lumi_evol, L[k1])
                Time=np.append(Time, T[k1])
                k1=k1+1
        
        ax.plot(T_zero, L_zero*1e30, "y.", label='Deleted zero', markersize=8)
        ax.plot(T_tofit2, L_tofit2*1e30, "g.", label='Deleted differences and derivatives', markersize=6)
        ax.plot(T_fit, L_fit*1e30, "b.", label='Fit interval', markersize=4)
        ax.plot(X,Y*1e30, "k.", label='Extrapolated Interval', markersize=4)
        ax.plot(Time, Lumi_evol*1e30, "r.", label='Luminosity Evolution', markersize=2)
        
        plt.legend(loc='best')
        ax.set_xlabel("Normalized Times")
        ax.set_ylabel("Luminosity [cm^-2 s^-1]")
        ax.set_title('Luminosity evolution of optimal fill {}'.format(text))
        plt.savefig('ATLAS/OptimalFillsLuminosityEvolution2016/{}_Opt_Lumi_Evol.pdf'.format(text)) 
        #plt.show()
         
     elif int(t_opt)<limit:       
        k2=0 
        Lumi_evol=[]
        Time=[]
        for m in Times:
            if int(t_opt)>=m:
                Lumi_evol.append(L_evol[k2])
                Time.append(Times[k2])
                k2=k2+1
            

        Lumi_evol = np.array(Lumi_evol).flatten('F')
        Time = np.array(Time).flatten('F')
        
        #normalizing the time interval 
        #norm_Time=[]
        #norm_Time=np.array(norm_Time)
        #for t in Time:
            #z=(t-np.amin(Time))/(np.amax(Time)-np.amin(Time))
            #norm_Time=np.append(norm_Time, z)
        #Time1=Time   
        #Time=norm_Time
        

        ax.plot(Time, Lumi_evol*1e30, "r.", label='Luminosity Evolution', markersize=2)
        #ax.legend(loc='best')
        ax.set_xlabel("Normalized Times")
        ax.set_ylabel("Luminosity [cm^-2 s^-1]")
        ax.set_title('Luminosity evolution of optimal fill {}'.format(text))
        plt.savefig('ATLAS/OptimalFillsLuminosityEvolution2016/{}_Opt_Lumi_Evol.pdf'.format(text)) 
        #plt.show()
     
     #evaluating the fill integral luminosity
     L_int_opt = integrate.simps(Lumi_evol*1e30, Time)
     L_int_opt2016.append(L_int_opt)

#2017   
L_int_opt2017=[]  
for i in range(len(FillNumber17)):
     print(i)
     n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps= lo.Parameters2017() #importing the theoretical luminosity model parameters
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
     
     #evaluate the optimal time for the current fill
     t_opt = lo.t_opt_eval(N_i, n_c, Xi, S_int, data_ta17_sec[i])
     print("t_opt=", t_opt,"Last time=", Times[len(Times)-1], "ta=", data_ta17_sec[i])
     
     #compare the optimal time with the last value of the Times array
     t_opt = Times[0]+t_opt
     print("Times[0]+t_opt=", t_opt)
     control = t_opt>Times[len(Times)-1]   
     print("Times[0]+t_opt=", t_opt, "Extrapolation needed:", control)
     limit=Times[len(Times)-1] #last value of the Times array
     
     
     #plotting results 
     plt.close("all")
     fig,  ax = plt.subplots()
     
     #performing extrapolation when needed
     if int(t_opt)>limit:
         
        X = np.arange(limit+60, int(t_opt), 60) #definition of the extrapolation interval
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
            
         #defining the fit intervals
        j1=0
        L_fit=np.array([])
        T_fit=np.array([])
        L_tofit2=np.array(L_tofit2)
        T_tofit2=np.array(T_tofit2)
        dy2 = np.zeros(L_tofit2.shape)
        dy2[0:-1] = np.diff(L_tofit2)/np.diff(T_tofit2)
        for a in range(len(L_tofit2)):
            if abs(dy2[len(dy2)-1-j1])<0.6: #and dy2[len(dy2)-1-j]>-2.5:
                L_fit=np.append(L_fit, L_tofit2[len(L_tofit2)-1-j1])
                T_fit=np.append(T_fit, T_tofit2[len(L_tofit2)-1-j1])
                j1=j1+1
            elif abs(dy2[len(dy2)-j1])>0.6: #and dy2[len(dy2)-1-j]>-2.5:
                break
        if len(L_fit)<3:
            L_fit=L_tofit2
            T_fit=T_tofit2

        #normalization of the fit interval    
        norm_T_fit=[]
        norm_T_fit=np.array(norm_T_fit)
        norm_X=[]
        norm_X=np.array(norm_X)
        for element in T_fit:
            z=(element-np.amin(T_fit))/(np.amax(X)-np.amin(T_fit))
            norm_T_fit=np.append(norm_T_fit, z)
    
        for element in X:
            z=(element-np.amin(T_fit))/(np.amax(X)-np.amin(T_fit))
            norm_X=np.append(norm_X, z)   
                    
        #performing fit of last segments of data
        fit_result=model.fit(L_fit, x=norm_T_fit, a=1, b=0.2, c=1, d=0.2)
        Y = fit(norm_X, fit_result.params['a'].value, fit_result.params['b'].value, fit_result.params['c'].value, fit_result.params['d'].value)
        print(fit_result.params['a'].value, fit_result.params['b'].value, fit_result.params['c'].value, fit_result.params['d'].value)
        #if the double exponential grows changing the fit model --> one exponential
        if fit_result.params['b'].value <0 or fit_result.params['d'].value <0: 
            print('Single exponential!')
            def fit_lin(x, a2, b2):
                return a2*np.exp((-b2)*x)

            model2=Model(fit_lin)
            fit_result2=model2.fit(L_fit, x=norm_T_fit, a2=1, b2=1)
            print(fit_result2.params['a2'].value, fit_result2.params['b2'].value)
            Y = fit_lin(norm_X, fit_result2.params['a2'].value, fit_result2.params['b2'].value)
        
        k1=0   
        L=np.append(L_zero, Y)
        T=np.append(T_zero, X)
        Lumi_evol=np.array([])
        Time=np.array([])
        for m in T:
            if int(t_opt)>=m:
                Lumi_evol=np.append(Lumi_evol, L[k1])
                Time=np.append(Time, T[k1])
                k1=k1+1
        
        ax.plot(T_zero, L_zero*1e30, "y.", label='Deleted zero', markersize=8)
        ax.plot(T_tofit2, L_tofit2*1e30, "g.", label='Deleted differences and derivatives', markersize=6)
        ax.plot(T_fit, L_fit*1e30, "b.", label='Fit interval', markersize=4)
        ax.plot(X,Y*1e30, "k.", label='Extrapolated Interval', markersize=4)
        ax.plot(Time, Lumi_evol*1e30, "r.", label='Luminosity Evolution', markersize=2)
        
        plt.legend(loc='best')
        ax.set_xlabel("Normalized Times")
        ax.set_ylabel("Luminosity [cm^-2 s^-1]")
        ax.set_title('Luminosity evolution of optimal fill {}'.format(text))
        plt.savefig('ATLAS/OptimalFillsLuminosityEvolution2017/{}_Opt_Lumi_Evol.pdf'.format(text)) 
        #plt.show()
         
     elif int(t_opt)<limit:       
        k2=0 
        Lumi_evol=[]
        Time=[]
        for m in Times:
            if int(t_opt)>=m:
                Lumi_evol.append(L_evol[k2])
                Time.append(Times[k2])
                k2=k2+1
            

        Lumi_evol = np.array(Lumi_evol).flatten('F')
        Time = np.array(Time).flatten('F')
        
        

        ax.plot(Time, Lumi_evol*1e30, "r.", label='Luminosity Evolution', markersize=2)
        #ax.legend(loc='best')
        ax.set_xlabel("Normalized Times")
        ax.set_ylabel("Luminosity [cm^-2 s^-1]")
        ax.set_title('Luminosity evolution of optimal fill {}'.format(text))
        plt.savefig('ATLAS/OptimalFillsLuminosityEvolution2017/{}_Opt_Lumi_Evol.pdf'.format(text)) 
        #plt.show()
     
     #evaluating the fill integral luminosity
     L_int_opt = integrate.simps(Lumi_evol*1e30, Time)
     L_int_opt2017.append(L_int_opt)
        
     
#2018
L_int_opt2018=[]
for i in range(len(FillNumber18)):
     print(i)
     n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps= lo.Parameters2018() #importing the theoretical luminosity model parameters
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
     
     #evaluate the optimal time for the current fill
     t_opt = lo.t_opt_eval(N_i, n_c, Xi, S_int, data_ta18_sec[i])
     print("t_opt=", t_opt,"Last time=", Times[len(Times)-1], "ta=", data_ta18_sec[i])
     
     #compare the optimal time with the last value of the Times array
     t_opt = Times[0]+t_opt
     print("Times[0]+t_opt=", t_opt)
     control = t_opt>Times[len(Times)-1]   
     print("Times[0]+t_opt=", t_opt, "Extrapolation needed:", control)
     limit=Times[len(Times)-1] #last value of the Times array
     
     #plotting results 
     plt.close("all")
     fig,  ax = plt.subplots()
     
     #performing extrapolation when needed
     if int(t_opt)>limit:
         
        X = np.arange(limit+60, int(t_opt), 60) #definition of the extrapolation interval
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
            
         #defining the fit intervals
        j2=0
        L_fit=np.array([])
        T_fit=np.array([])
        L_tofit2=np.array(L_tofit2)
        T_tofit2=np.array(T_tofit2)
        dy2 = np.zeros(L_tofit2.shape)
        dy2[0:-1] = np.diff(L_tofit2)/np.diff(T_tofit2)
        for a in range(len(L_tofit2)):
            if abs(dy2[len(dy2)-1-j2])<0.6: 
                L_fit=np.append(L_fit, L_tofit2[len(L_tofit2)-1-j2])
                T_fit=np.append(T_fit, T_tofit2[len(L_tofit2)-1-j2])
                j2=j2+1
            elif abs(dy2[len(dy2)-j2])>0.6: 
                break
        if len(L_fit)<3:
            L_fit=L_tofit2
            T_fit=T_tofit2

        #normalization of the fit interval    
        norm_T_fit=[]
        norm_T_fit=np.array(norm_T_fit)
        norm_X=[]
        norm_X=np.array(norm_X)
        for element in T_fit:
            z=(element-np.amin(T_fit))/(np.amax(X)-np.amin(T_fit))
            norm_T_fit=np.append(norm_T_fit, z)
    
        for element in X:
            z=(element-np.amin(T_fit))/(np.amax(X)-np.amin(T_fit))
            norm_X=np.append(norm_X, z)   
                    
        #performing fit of last segments of data
        fit_result=model.fit(L_fit, x=norm_T_fit, a=1, b=0.2, c=1, d=0.2)
        Y = fit(norm_X, fit_result.params['a'].value, fit_result.params['b'].value, fit_result.params['c'].value, fit_result.params['d'].value)
        print(fit_result.params['a'].value, fit_result.params['b'].value, fit_result.params['c'].value, fit_result.params['d'].value)
        #if the double exponential grows changing the fit model --> one exponential
        if fit_result.params['b'].value <0 or fit_result.params['d'].value <0: 
            print('Single exponential!')
            def fit_lin(x, a2, b2):
                return a2*np.exp((-b2)*x)

            model2=Model(fit_lin)
            fit_result2=model2.fit(L_fit, x=norm_T_fit, a2=1, b2=1)
            print(fit_result2.params['a2'].value, fit_result2.params['b2'].value)
            Y = fit_lin(norm_X, fit_result2.params['a2'].value, fit_result2.params['b2'].value)
        
        k1=0   
        L=np.append(L_zero, Y)
        T=np.append(T_zero, X)
        Lumi_evol=np.array([])
        Time=np.array([])
        for m in T:
            if int(t_opt)>=m:
                Lumi_evol=np.append(Lumi_evol, L[k1])
                Time=np.append(Time, T[k1])
                k1=k1+1
        
        ax.plot(T_zero, L_zero*1e30, "y.", label='Deleted zero', markersize=8)
        ax.plot(T_tofit2, L_tofit2*1e30, "g.", label='Deleted differences and derivatives', markersize=6)
        ax.plot(T_fit, L_fit*1e30, "b.", label='Fit interval', markersize=4)
        ax.plot(X,Y*1e30, "k.", label='Extrapolated Interval', markersize=4)
        ax.plot(Time, Lumi_evol*1e30, "r.", label='Luminosity Evolution', markersize=2)
        
        plt.legend(loc='best')
        ax.set_xlabel("Normalized Times")
        ax.set_ylabel("Luminosity [cm^-2 s^-1]")
        ax.set_title('Luminosity evolution of optimal fill {}'.format(text))
        plt.savefig('ATLAS/OptimalFillsLuminosityEvolution2018/{}_Opt_Lumi_Evol.pdf'.format(text)) 
        #plt.show()
         
     elif int(t_opt)<limit:       
        k2=0 
        Lumi_evol=[]
        Time=[]
        for m in Times:
            if int(t_opt)>=m:
                Lumi_evol.append(L_evol[k2])
                Time.append(Times[k2])
                k2=k2+1
            

        Lumi_evol = np.array(Lumi_evol).flatten('F')
        Time = np.array(Time).flatten('F')
        
        

        ax.plot(Time, Lumi_evol*1e30, "r.", label='Luminosity Evolution', markersize=2)
        #ax.legend(loc='best')
        ax.set_xlabel("Normalized Times")
        ax.set_ylabel("Luminosity [cm^-2 s^-1]")
        ax.set_title('Luminosity evolution of optimal fill {}'.format(text))
        plt.savefig('ATLAS/OptimalFillsLuminosityEvolution2018/{}_Opt_Lumi_Evol.pdf'.format(text)) 
        #plt.show()
     
     #evaluating the fill integral luminosity
     L_int_opt = integrate.simps(Lumi_evol*1e30, Time)
     L_int_opt2018.append(L_int_opt)
        
        
     
#defining the dataframes 
df16=pd.DataFrame(L_int_opt2016, columns=['Opt. Measured Int. Lum.']) 
df17=pd.DataFrame(L_int_opt2017, columns=['Opt. Measured Int. Lum.'])  
df18=pd.DataFrame(L_int_opt2018, columns=['Opt. Measured Int. Lum.']) 


#evaluating the total luminosities   
L_tot_2016=[np.sum(L_int_opt2016)]
L_tot_2017=[np.sum(L_int_opt2017)]
L_tot_2018=[np.sum(L_int_opt2018)]

df16tot=pd.DataFrame(L_tot_2016, columns=['Opt. Measured tot. Lum.']) 
df17tot=pd.DataFrame(L_tot_2017, columns=['Opt. Measured tot. Lum.'])  
df18tot=pd.DataFrame(L_tot_2018, columns=['Opt. Measured tot. Lum.']) 

#write dataframes on an excel file
with pd.ExcelWriter('Optimized Measured Luminosity.xlsx') as writer:
        df16.to_excel(writer, sheet_name='2016 Integrated Luminosity')
        df17.to_excel(writer, sheet_name='2017 Integrated Luminosity')
        df18.to_excel(writer, sheet_name='2018 Integrated Luminosity')  
        df16tot.to_excel(writer, sheet_name='2016 Total Luminosity')
        df17tot.to_excel(writer, sheet_name='2017 Total Luminosity')
        df18tot.to_excel(writer, sheet_name='2018 Total Luminosity')

#defining the stop time of the program      
stop=time.time()
#evaluating the run time 
runtime_seq=stop-start
print('The runtime is:', runtime_seq, '[s]')   