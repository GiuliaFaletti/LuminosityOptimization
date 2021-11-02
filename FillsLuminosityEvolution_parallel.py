##############################################################################################################################################
#  @Giulia Faletti
# Creating plots of the luminosity evolution of 
# physics fills at LHC Run 2
# ___PARALLELIZED CODE____
##############################################################################################################################################
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import LoadData as ld
import multiprocessing as mp
import time as t

#defining the start time of the program
start=t.time()
     
FillNumber16, FillNumber17, FillNumber18 = ld.FillNumber()

#2016
def func16(i):
     text = str(int(i))#number of current fill
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
     
     norm_time=[]
     norm_time=np.array(norm_time)
     for element in Times:
         z=(element-np.amin(Times))/(np.amax(Times)-np.amin(Times))
         norm_time=np.append(norm_time, z)
         
     Times=norm_time  
       
     plt.close("all")
     fig,  ax = plt.subplots()
     ax.plot(Times, L_evol*1e30, "r+", label='Luminosity Evolution')
     ax.set_xlabel("Normalized Times")
     ax.set_ylabel("Luminosity [cm^2 s^-1]")
     ax.set_title('Luminosity evolution of fill {}'.format(text))
     return plt.savefig('ATLAS/FillsLuminosityEvolution2016/{}_Lumi_Evol.pdf'.format(text))
          
#2017
def func17(i):
     text = str(int(i))#number of current fill
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
     
     #normalization of the time interval 
     norm_time=[]
     norm_time=np.array(norm_time)
     for element in Times:
         z=(element-np.amin(Times))/(np.amax(Times)-np.amin(Times))
         norm_time=np.append(norm_time, z)
         
     Times=norm_time  
        
     #plotting results  
     plt.close("all")
     fig,  ax = plt.subplots()
     ax.plot(Times, L_evol*1e30, "r+", label='Luminosity Evolution')
     ax.set_xlabel("Normalized Times")
     ax.set_ylabel("Luminosity [cm^2 s^-1]")
     ax.set_title('Luminosity evolution of fill {}'.format(text))
     plt.savefig('ATLAS/FillsLuminosityEvolution2017/{}_Lumi_Evol.pdf'.format(text))
     
#2018    
def func18(i):
     text = str(int(i))#number of current fill
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

     norm_time=[]
     norm_time=np.array(norm_time)
     for element in Times:
         z=(element-np.amin(Times))/(np.amax(Times)-np.amin(Times))
         norm_time=np.append(norm_time, z)
         
     Times=norm_time    
          
     plt.close("all")
     fig,  ax = plt.subplots()
     ax.plot(Times, L_evol*1e30, "r+", label='Luminosity Evolution')
     #ax.legend(loc='best')
     ax.set_xlabel("Normalized Times")
     ax.set_ylabel("Luminosity [cm^2 s^-1]")
     ax.set_title('Luminosity evolution of fill {}'.format(text))
     plt.savefig('ATLAS/FillsLuminosityEvolution2018/{}_Lumi_Evol.pdf'.format(text))         

#implement the main of the program including the parallelization functions
if __name__ == '__main__':
    pool=mp.Pool(mp.cpu_count()) #defining the Pool for parallelization 
    #parallelizing the plot production
    pool.map(func16, [i for i in FillNumber16]) #parallelization of 2016
    pool.map(func17, [i for i in FillNumber17]) #parallelization of 2017
    pool.map(func18, [i for i in FillNumber18]) #parallelization of 2018
    
    #defining the stop time of the program      
    stop=t.time()
    #evaluating the run time 
    runtime_parallel=stop-start
    print('The runtime of the parallelized code is:', runtime_parallel, '[s]')
    pool.close() #close the pool
    pool.join()


