import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import LoadData as ld
import LuminosityOptimization as lo
import scipy.integrate as integrate
from lmfit import Model
import DataModel as dm

#loading data
data_16, data_17, data_18, array16, array17, array18 = ld.Data()
data_tot, dataTot, arrayTot = ld.TotalDataSet(data_16, data_17, data_18)
data_ta16, data_tf16, data_ta17, data_tf17, data_ta18, data_tf18 = ld.loadFill()     
FillNumber16, FillNumber17, FillNumber18 = ld.FillNumber()
data_16_sec, data_ta16_sec, data_tf16_sec, data_17_sec, data_ta17_sec,\
   data_tf17_sec, data_18_sec, data_ta18_sec, data_tf18_sec=ld.Data_sec(array16,\
      data_ta16, data_tf16, array17, data_ta17, data_tf17, array18, data_ta18, data_tf18)
   
L_mes16, L_mes17, L_mes18=ld.MeasuredLuminosity()   

L_int16=[]
for i in range(len(FillNumber16)):
     print(i)
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
     L_int_opt = integrate.simps(L_evol, Times)
     L_int16.append(L_int_opt)
     
L_int17=[]
for i in range(len(FillNumber17)):
     print(i)
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
     L_int_opt = integrate.simps(L_evol, Times)
     L_int17.append(L_int_opt)


L_int18=[]
for i in range(len(FillNumber18)):
     print(i)
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
     L_int_opt = integrate.simps(L_evol, Times)
     L_int18.append(L_int_opt)  
     
     
    
L_inst16, L_intmodel16, L_tot16=lo.Model_L16(data_tf16_sec, data_ta16_sec)  
L_inst17, L_intmodel17, L_tot17=lo.Model_L17(data_tf17_sec, data_ta17_sec)  
L_inst18, L_intmodel18, L_tot18=lo.Model_L18(data_tf18_sec, data_ta18_sec)   
print('******2016*******')        
print(L_mes16[0], L_int16[0]/1e9, L_intmodel16[0]/1e43, dm.L_intfit16[0]/1e9)
print(L_mes16[int(len(L_mes16)/2)], L_int16[int(len(L_int16)/2)]/1e9, L_intmodel16[int(len(L_intmodel16)/2)]/1e43, dm.L_intfit16[int(len(dm.L_intfit16)/2)]/1e9)
print(L_mes16[-1], L_int16[-1]/1e9, L_intmodel16[-1]/1e43, dm.L_intfit16[-1]/1e9 )
print(len(L_mes16), len(L_int16), len(L_intmodel16), len(dm.L_intfit16))
print('******2017*******')        
print(L_mes17[0], L_int17[0]/1e9, L_intmodel17[0]/1e43, dm.L_intfit17[0]/1e9)
print(L_mes17[int(len(L_mes17)/2)], L_int17[int(len(L_int17)/2)]/1e9, L_intmodel17[int(len(L_intmodel17)/2)]/1e43, dm.L_intfit17[int(len(dm.L_intfit17)/2)]/1e9)
print(L_mes17[-1], L_int17[-1]/1e9, L_intmodel17[-1]/1e43, dm.L_intfit17[-1]/1e9)
print(len(L_mes16), len(L_int16), len(L_intmodel17), len(dm.L_intfit17))
print('******2018*******')        
print(L_mes18[0], L_int18[0]/1e9, L_intmodel18[0]/1e43, dm.L_intfit18[0]/1e9)
print(L_mes18[int(len(L_mes18)/2)], L_int18[int(len(L_int18)/2)]/1e9, L_intmodel18[int(len(L_intmodel18)/2)]/1e43, dm.L_intfit18[int(len(dm.L_intfit18)/2)]/1e9)
print(L_mes18[-1], L_int18[-1]/1e9, L_intmodel18[-1]/1e43, dm.L_intfit18[-1]/1e9)
print(len(L_mes18), len(L_int18), len(L_intmodel18), len(dm.L_intfit18))