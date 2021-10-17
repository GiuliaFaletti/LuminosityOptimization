import numpy as np
from numpy.ma import array
import LoadData as ld
import LuminosityOptimization as lo
import matplotlib.pyplot as plt
import matplotlib as mpl
from art import *
from scipy import integrate
from lmfit import Model
import CreatingVariableBins as cvb
import Models as mod
from scipy.stats import ks_2samp
import pandas as pd


mpl.style.use('default')
plt.close("all")

#Loading DataSets
data_16, data_17, data_18, array16, array17, array18 = ld.Data()
data_tot, dataTot, arrayTot = ld.TotalDataSet(data_16, data_17, data_18)
data_ta16, data_tf16, data_ta17, data_tf17, data_ta18, data_tf18 = ld.loadFill()
FillNumber16, FillNumber17, FillNumber18 = ld.FillNumber()
L_mes16, L_mes17, L_mes18 = ld.MeasuredLuminosity()

#Transforming data from hours into seconds - Total Data set
data_tot_sec = []
for i in data_tot:
    i = i*3600
    data_tot_sec.append(i)   

MaxT = np.max(data_tot_sec)
MinT = np.min(data_tot_sec)
data_tot_sec = np.array(data_tot_sec)

#Transforming data from hours into seconds
data_16_sec, data_ta16_sec, data_tf16_sec, data_17_sec, data_ta17_sec,\
   data_tf17_sec, data_18_sec, data_ta18_sec, data_tf18_sec=ld.Data_sec(array16,\
      data_ta16, data_tf16, array17, data_ta17, data_tf17, array18, data_ta18, data_tf18)

#welcome
tprint("Optimize",font="double")
print("________________________________________________________________")
print("|                     WELCOME To OPTIMIZE                      |")
print("| Optimize is a program which produces the result of a study   |")
print("| on luminosity optimization for LHC. The program is dvided    |")
print("| in different parts: the first presents the results obtained  |")
print("| analyzing the LHC Run 2 data; the second one allows the user |")
print("| to evaluate the optimal values of fill time, inserting the   |")
print("| correct parameters for the luminosity evolution model;       |")
print("| the third part shows the statistical studies on turnaround   |")
print("| times.                                                       |")
print("| @Optimizing the LHC luminosity - GiuliaFaletti               |")
print("________________________________________________________________")

#first choice: 1: LHC Run 2 Optimization studies; 2: New Run; 3: Turn Around statistical analysis
CONTROL1=True
while CONTROL1:
     control1 = input("###########################################\n"\
                      "-Enter 1: LHC Run 2 Optimization studies;\n-Enter 2: New Run;\n-Enter 3: Turn Around statistical analysis;\n-Enter q to Exit.\n")
     
     #LHC Run 2 Optimization studies
     if control1 == "1":
        #second choice: 1: 2016 Preliminary Analysis; 2: 2017 Preliminary Analysis; 3: 2018 Preliminary Analysis; 4: 2016 Vincolated Analysis;
        CONTROL2 = True
        while CONTROL2:
            control2 = input("#####################\n"\
                             "-Enter 1: 2016 Preliminary Analysis;\n-Enter 2: 2017 Preliminary Analysis;\n-Enter 3: 2018 Preliminary Analysis;\n-Enter 4: 2016 Vincolated Analysis;\n-Enter 5: 2017 Vincolated Analysis;\n-Enter 6: 2018 Vincolated Analysis;\nEnter b to go back to the previous menu.\n")
            
            #2016 Preliminary Analysis
            if control2 == "1":
               #Defining the turn around times in seconds
               Max = np.max(data_16_sec)
               Min  = np.min(data_16_sec)
               t_a = np.linspace(Min, Max, 100) #turn around time in seconds

               #Defining the turn around times in hours
               max = np.max(data_16)
               min  = np.min(data_16)
               tta = np.linspace(min, max, 100) #turn around time in hours  
               
               #defining models parameters
               n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps= lo.Parameters2016()
               
               #sorted data array
               ttta = lo.selection_sort(data_ta16_sec)
               
               #evaluating optimal quantities
               t_opt, L_tot_opt, L_int_opt, t_a16, t_opt16, L_int16, L_tot16 = lo.L_optimal_16(t_a)
               t_opt_real,  L_opt_real, L_int_Opt, t_a16, t_opt16, L_int16, L_tot16= lo.L_optimal_16(ttta)
               L_inst_real, L_int_real, L_tot_real= lo.Model_L16(data_tf16_sec, data_ta16_sec) 
               Inst_x, Integ_x, Tot_x= lo.Model_L16(t_opt16, t_a16)
               
               CONTROL3=True
               while CONTROL3:
                   control3 = input("##################################################################################################################################\n"\
                       "-Enter 1: Print the model parameters;\n"\
                       "-Enter 2: Evaluate the Peak Luminosity;\n"\
                       "-Enter 3: Comparison between Optimal Fill Times and Fill Time Data;\n"\
                       "-Enter 4: Comparison between Optimal Integrated Luminosity and Integrated Luminosity Data;\n"\
                       "-Enter 5: Distribution of Optimal Fill Times and Fill Time Data;\n"\
                       "-Enter 6: Distribution of Optimal Integrated Luminosity and Integrated Luminosity data;\n" 
                       "-Enter 7: Print the Average turn around time of the year and its respective optimal fill time and total luminosity;\n"    
                       "-Enter b to go back to the previous menu.\n"   )
                   
                   #Print the model parameters
                   if control3 == "1":
                        print("______________________________________________________________________________________________")
                        print("|_____________________________________Model Parameters________________________________________")
                        print("| Bunch Population                                         n_i =",n_i/1e11,"e+11")
                        print("| Value of the beta-function at the collision point [m]    B_* =",B_s)
                        print("| RMS normalized transverse emittance               [m]    E_* =",E_s)
                        print("| Relativistic beta factor                                 B_r =",B_r)
                        print("| Relativistic gamma factor                                G_r =", G_r)
                        print("| Half Crossing Angle                               [rad]  T_hc =", T_hc)
                        print("| Longitudinal RMS dimension                        [m]    S_z =", S_z)
                        print("| Tranverse RMS dimension                           [m]    S_s =", S_s)
                        print("| Interaction cross section                         [rad]  S_int =", S_int)
                        print("| Revolution Frequency                              [Hz]   f_rev =", f_rev)
                        print("| Reduction in volume overlap factor                [m]    Fe =", Fe," ")
                        print("| HL Collision points                                      n_c =", n_c)
                        print("| Colliding Bunches                                        k_b =", k_b)
                        print("|                                                          \u039E =", Xi)
                        print("|                                                          \u03B5 =", Eps)
                        print("| Total Physics Time                                [s]    T_ph =", T_ph)
                        print("_____________________________________________________________________________________________")
                   
                   #Evaluate the Peak Luminosity
                   elif control3 == "2":
                       t_ff = 0
                       L_inst, L_int, L_tot = lo.Model_L16(t_ff, t_a)
                       print("________________________________________________________________")
                       print("_______________________Peak Luminosity__________________________")
                       print("|                                                               ")
                       print("| L_model(0)) = ", L_inst*1e-4, "[cm^-2 s^-1]")
                       print("| L_LHC(0) = 1.4e+34 cm^-2 s^-1")
                       print("________________________________________________________________")
                   
                   #Comparison between Optimal Fill Times and Fill Time Data
                   elif control3 == "3":
                       plt.close("all")
                       fig1, ax1 = plt.subplots()
                       ax1.plot(tta, t_opt/3600, 'k-', label="Optimal Fill Times")
                       ax1.plot([2.3,23], [18.459584259140506, 18.459584259140506], 'r-', label="Optimal Fill Time for <ta>")
                       ax1.plot(data_ta16, data_tf16, 'b+', label= "2016 Data")
                       ax1.set_title('Comparison between Optimal Fill Times and 2016 Fill Time Data')
                       ax1.set_xlabel('Turn Around Times [h]')
                       ax1.set_ylabel('Fill Times [h]')
                       ax1.legend(loc='best')
                       plt.savefig('1_Immagini/2016/Optimal_tf-Data_tf.pdf')
                       plt.show()
                       
                   #Comparison between Optimal Integrated Luminosity and Integrated Luminosity Data
                   elif control3 == "4":
                       plt.close("all")
                       fig, ax = plt.subplots()
                       ax.plot(ttta/3600, L_int_Opt/1e43, 'k-', label="Optimal Luminosity")
                       ax.plot([2.3,23], [Integ_x/1e43, Integ_x/1e43], 'r-', label="Optimal Luminosity for <ta>")
                       ax.plot(data_ta16, L_mes16, 'b+', label="2016 Data")
                       ax.set_ylabel('Integrated Luminosity [fb^-1]')
                       ax.set_xlabel('Turn Around Times [h]')
                       ax.set_title('Comparison between Optimal Integrated Luminosity\nand 2016 Integrated Luminosity Data')
                       ax.legend(loc='upper right')
                       plt.savefig('1_Immagini/2016/Optimal_L_Int-Data_L_Int.pdf')
                       plt.show()
                       
                   #Distribution of Optimal Fill Times and Fill Time Data
                   elif control3 == "5":
                       c=N_i*n_c*Xi*S_int
                       a = 0.37082720
                       b = 2*3600
                       d = 0.1/3600
                       n = 0.73454371
                       f = lambda x: (a/(np.power((x-b), n)))*np.exp(-d*x)
                       val, err=integrate.quad(f, Min, Max)
                       g = lambda t_opt: ((a*2*c*t_opt)/(np.power((c*(np.power(t_opt, 2))-b), n)))*np.exp(-d*c*(np.power(t_opt, 2)))
                       a = 0.37082720/val 
                       
                       plt.close("all")
                       ###### HISTOGRAMS OF the distribution of t_opt_real and L_opt_real
                       fig5,  ax5 = plt.subplots(figsize=(5,4))
                       n5, bins5, patches5 = ax5.hist(t_opt_real/3600, bins=15, alpha=0.4, density=True, label="Optimal Fill Times")
                       ax5.plot(t_opt/3600, g(t_opt)*3600, "r-", label='Theorical Distribution of t_opt')
                       ax5.legend(loc='best')
                       ax5.set_xlabel("Optimal Fill Times [h]")
                       ax5.set_ylabel("Normalized Frequencies")
                       ax5.set_title("2016")
                       plt.savefig('1_Immagini/2016/t_opt_theory-distribution.pdf')
                       plt.show()
                     
                       T_ph_opt= np.sum(t_opt_real)
                       T_ph_opt = T_ph_opt/3600
                       T_ph_real= np.sum(data_tf16_sec)
                       T_ph_real= T_ph_real/3600
                       print("________________________________________________________________")
                       print("____________________Total Physics time__________________________")
                       print("|                                                               ")
                       print("| T_ph_opt = ", T_ph_opt/24, "days")
                       print("| T_ph_data =", T_ph_real/24, "days")
                       print("|")
                       print("________________________________________________________________")
                       fig4,  ax4 = plt.subplots(figsize=(4,4))
                       n4, bins4, patches4 = ax4.hist(t_opt_real/3600, bins=np.arange(1, 40 + 2.3, 2.3), histtype='step', alpha=0.2, fill=True, label="Optimal Fill Times")
                       n4, bins4, patches4 = ax4.hist(data_tf16_sec/3600, bins=np.arange(1, 40 + 2.3, 2.3), histtype='step', fill=False, label="2016 Data")
                       ax4.legend(loc='best')
                       ax4.set_xlabel("Fill Times [h]")
                       ax4.set_ylabel("Number of Fills")
                       plt.savefig('1_Immagini/2016/t_opt_Distributions.pdf')

                       plt.show()
                       
                   #Distribution of Optimal Integrated Luminosity and Integrated Luminosity data
                   elif control3 == "6":
                       
                       L_tot1 = np.sum(L_int_Opt)
                       L_tot2 = np.sum(L_int_real)
                       L_tot3 = np.sum(L_mes16)
                       testo1 = 'L_tot_opt={:.2f} [fb^-1]'.format(L_tot1/1e43)
                       testo2 = 'L_tot_eval={:.2f} [fb^-1]'.format(L_tot2/1e43)
                       testo3 = 'L_tot_data={:.2f} [fb^-1]'.format(L_tot3)
                       plt.close("all")
                       ###### HISTOGRAMS OF the distribution of t_opt_real and L_opt_real
                       fig6, ax6 = plt.subplots()
                       n6, bins6, patches6 = ax6.hist(L_int_Opt/1e43, bins=np.arange(0.1, 1 + 0.05, 0.05), facecolor='steelblue', density=True, alpha=0.4, label="Optimal Integrated Luminosity")
                       #n7, bins7, patches7 = ax6.hist(L_int_real/1e43, bins=np.arange(0.1, 1 + 0.05, 0.05), color='red', histtype='step', density=True, label="Integrated Luminosity evaluated\n from data")
                       n8, bins8, patches8 = ax6.hist(L_mes16, bins=np.arange(0.1, 1 + 0.05, 0.05), color='red', histtype='step', density=True, label="Mesured Integrated Luminosity")
                       ax6.text(0.6, 3.5, testo1, bbox=dict(facecolor='none', edgecolor='steelblue', boxstyle='round,pad=0.5'))
                       #ax6.text(0.6, 2.8, testo2, bbox=dict(facecolor='none', edgecolor='red', boxstyle='round,pad=0.5'))
                       ax6.text(0.6, 2.9, testo3, bbox=dict(facecolor='none', edgecolor='red', boxstyle='round,pad=0.5'))
                       ax6.set_xlabel("Integrated Luminosity [fb^-1]")
                       ax6.set_ylabel("Normalized Frequencies")
                       ax6.legend(loc='best')
                       ax6.set_title("2016")
                       plt.savefig('1_Immagini/2016/L_Int-Data_L_Int_Distributions.pdf') 
                       plt.show()

                       print("________________________________________________________________")
                       print("______Comparison between total optimal luminosity_______________")
                       print("________and total luminosity evaluated from data________________")
                       print("____________________and the measured one________________________")
                       print("|")
                       print("| L_tot_opt=", L_tot1/1e43, "[fb^-1]")
                       print("|")
                       print("| L_tot_eval=", L_tot2/1e43, "[fb^-1]")
                       print("|")
                       print("| L_tot_data=", L_tot3, "[fb^-1]")
                       print("|")
                       print("________________________________________________________________")
                                        
                   #Print the Average turn around time of the year and its respective optimal fill time and total luminosity
                   elif control3 == "7":
                       print("________________________________________________________________")
                       print("_____________________Average Turn Around Time___________________")
                       print("|")
                       print("| <ta>=", t_a16/3600, "[h]")
                       print("|")
                       print("________________________Optimal Fill Time_______________________")
                       print("|")
                       print("| t_opt=", t_opt16/3600, "[h]")
                       print("|")
                       print("_________________Optimal Integrated Luminosity__________________")
                       print("|")
                       print("| L_int=", L_int16/1e43, "[fb^-1]")
                       print("|")
                       print("________________________________________________________________")
                       print("_________________Optimal Total Luminosity_______________________")
                       print("|")
                       print("| L_opt=", L_tot16/1e43, "[fb^-1]")
                       print("|")
                       print("________________________________________________________________")
                       
                   elif control3 == "b"or"q":
                       CONTROL3=False
                   else:
                       print("You did not enter a valid command!") 
            
            #2017 Preliminary Analysis                                   
            elif control2 == "2":
               #Defining the turn around times in seconds
               Max = np.max(data_17_sec)
               Min  = np.min(data_17_sec)
               t_a = np.linspace(Min, Max, 100) #turn around time in seconds

               #Defining the turn around times in hours
               max = np.max(data_17)
               min  = np.min(data_17)
               tta = np.linspace(min, max, 100) #turn around time in hours 
               
               #model parameters
               n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps= lo.Parameters2017()
               
               #sorted data array
               ttta = lo.selection_sort(data_ta17_sec)
               
               #evaluating optimized quantities
               t_opt, L_tot_opt, L_int_opt, t_a17, t_opt17, L_int17, L_tot17 = lo.L_optimal_17(t_a)
               t_opt_real,  L_opt_real, L_int_Opt, t_a17, t_opt17, L_int17, L_tot17= lo.L_optimal_17(ttta)
               L_inst_real, L_int_real, L_tot_real= lo.Model_L17(data_tf17_sec, data_ta17_sec) 
               Inst_x, Integ_x, Tot_x= lo.Model_L17(t_opt17, t_a17)
                       
               CONTROL3=True
               while CONTROL3:
                   control3 = input("##################################################################################################################################\n"\
                       "-Enter 1: Print the model parameters;\n"\
                       "-Enter 2: Evaluate the Peak Luminosity;\n"\
                       "-Enter 3: Comparison between Optimal Fill Times and Fill Time Data;\n"\
                       "-Enter 4: Comparison between Optimal Integrated Luminosity and Integrated Luminosity Data;\n"\
                       "-Enter 5: Distribution of Optimal Fill Times and Fill Time Data;\n"\
                       "-Enter 6: Distribution of Optimal Integrated Luminosity and Integrated Luminosity data;\n" 
                       "-Enter 7: Print the Average turn around time of the year and its respective optimal fill time and total luminosity;\n"    
                       "-Enter b to go back to the previous menu.\n"   )
                   
                   #Print the model parameters
                   if control3 == "1":
                        print("______________________________________________________________________________________________")
                        print("|_____________________________________Model Parameters________________________________________")
                        print("| Bunch Population                                         n_i =",n_i/1e11,"e+11")
                        print("| Value of the beta-function at the collision point [m]    B_* =",B_s)
                        print("| RMS normalized transverse emittance               [m]    E_* =",E_s)
                        print("| Relativistic beta factor                                 B_r =",B_r)
                        print("| Relativistic gamma factor                                G_r =", G_r)
                        print("| Half Crossing Angle                               [rad]  T_hc =", T_hc)
                        print("| Longitudinal RMS dimension                        [m]    S_z =", S_z)
                        print("| Tranverse RMS dimension                           [m]    S_s =", S_s)
                        print("| Interaction cross section                         [rad]  S_int =", S_int)
                        print("| Revolution Frequency                              [Hz]   f_rev =", f_rev)
                        print("| Reduction in volume overlap factor                [m]    Fe =", Fe," ")
                        print("| HL Collision points                                      n_c =", n_c)
                        print("| Colliding Bunches                                        k_b =", k_b)
                        print("|                                                          \u039E =", Xi)
                        print("|                                                          \u03B5 =", Eps)
                        print("| Total Physics Time                                [s]    T_ph =", T_ph)
                        print("_____________________________________________________________________________________________")
                        
                   #Evaluate the Peak Luminosity
                   elif control3 == "2":
                       t_ff = 0
                       L_inst, L_int, L_tot = lo.Model_L17(t_ff, t_a)
                       print("________________________________________________________________")
                       print("_______________________Peak Luminosity__________________________")
                       print("|                                                               ")
                       print("| L_model(0)) = ", L_inst*1e-4, "[cm^-2 s^-1]")
                       print("| L_LHC(0) = 1.7e+34 cm^-2 s^-1")
                       print("________________________________________________________________")
                       
                   #Comparison between Optimal Fill Times and Fill Time Data
                   elif control3 == "3":
                       plt.close("all")
                       fig1, ax1 = plt.subplots()
                       ax1.plot(tta, t_opt/3600, 'k-', label="Optimal Fill Times")
                       ax1.plot([1.7,30], [17.191537228560026 , 17.191537228560026 ], 'r-', label="Optimal Fill Time for <ta>")
                       ax1.plot(data_ta17, data_tf17, 'b+', label= "2017 Data")
                       ax1.set_title('Comparison between Optimal Fill Times and 2017 Fill Time Data')
                       ax1.set_xlabel('Turn Around Times [h]')
                       ax1.set_ylabel('Fill Times [h]')
                       ax1.legend(loc='best')
                       plt.savefig('1_Immagini/2017/Optimal_tf-Data_tf.pdf')
                       plt.show()
                       
                   #Comparison between Optimal Integrated Luminosity and Integrated Luminosity Data
                   elif control3 == "4":
                       plt.close("all")
                       fig, ax = plt.subplots()
                       ax.plot(ttta/3600, L_int_Opt/1e43, 'k-', label="Optimal Luminosity")
                       ax.plot([2.3,23], [Integ_x/1e43, Integ_x/1e43], 'r-', label="Optimal Luminosity for <ta>")
                       ax.plot(data_ta17, L_mes17, 'b+', label="2017 Data")
                       ax.set_ylabel('Integrated Luminosity [fb^-1]')
                       ax.set_xlabel('Turn Around Times [h]')
                       ax.set_title('Comparison between Optimal Integrated Luminosity\nand 2017 Integrated Luminosity Data')
                       ax.legend(loc='lower right')
                       plt.savefig('1_Immagini/2017/Optimal_L_Int-Data_L_Int.pdf')
                       plt.show()
                       
                   #Distribution of Optimal Fill Times and Fill Time Data
                   elif control3 == "5":
                       c=N_i*n_c*Xi*S_int
                       a = 0.29687182
                       b = 2*3600
                       d = 0.10422466/3600
                       n = 0.50729741
                       f = lambda x: (a/(np.power((x-b), n)))*np.exp(-d*x)
                       val, err=integrate.quad(f, Min, Max)
                       g = lambda t_opt: ((a*2*c*t_opt)/(np.power((c*(np.power(t_opt, 2))-b), n)))*np.exp(-d*c*(np.power(t_opt, 2)))
                       a = 0.29687182/val 
                       
                       plt.close("all")
                       ###### HISTOGRAMS OF the distribution of t_opt_real and L_opt_real
                       fig5,  ax5 = plt.subplots(figsize=(5,4))
                       n5, bins5, patches5 = ax5.hist(t_opt_real/3600, bins=15, alpha=0.4, density=True, label="Optimal Fill Times")
                       ax5.plot(t_opt/3600, g(t_opt)*3600, "r-", label='Theorical Distribution of t_opt')
                       ax5.legend(loc='best')
                       ax5.set_xlabel("Optimal Fill Times [h]")
                       ax5.set_ylabel("Normalized Frequencies")
                       ax5.set_title("2017")
                       plt.savefig('1_Immagini/2017/t_opt_theory-distribution.pdf')
                       
                       T_ph_opt= np.sum(t_opt_real)
                       T_ph_opt = T_ph_opt/3600
                       T_ph_real= np.sum(data_tf17_sec)
                       T_ph_real=  T_ph_real/3600
                       print("________________________________________________________________")
                       print("____________________Total Physics time__________________________")
                       print("|                                                               ")
                       print("| T_ph_opt = ", T_ph_opt/24, "days")
                       print("| T_ph_data =", T_ph_real/24, "days")
                       print("|")
                       print("________________________________________________________________")
                       
                       fig4,  ax4 = plt.subplots(figsize=(4,4))
                       n4, bins4, patches4 = ax4.hist(t_opt_real/3600, bins=np.arange(1, 30 + 2, 2), histtype='step', alpha=0.2, fill=True, label="Optimal Fill Times")
                       n4, bins4, patches4 = ax4.hist(data_tf17_sec/3600, bins=np.arange(1, 30 + 2, 2), histtype='step', fill=False, label="2017 Data")
                       ax4.legend(loc='best')
                       ax4.set_xlabel("Fill Times [h]")
                       ax4.set_ylabel("Number of Fills")
                       plt.savefig('1_Immagini/2017/t_opt_Distributions.pdf')
                       plt.show()
                       
                   #Distribution of Optimal Integrated Luminosity and Integrated Luminosity data
                   elif control3 == "6":
                       L_tot1 = np.sum(L_int_Opt)
                       L_tot2 = np.sum(L_int_real)
                       L_tot3 = np.sum(L_mes17)
                       testo1 = 'L_tot_opt={:.2f} [fb^-1]'.format(L_tot1/1e43)
                       testo2 = 'L_tot_data={:.2f} [fb^-1]'.format(L_tot3)
                       plt.close("all")
                       ###### HISTOGRAMS OF the distribution of t_opt_real and L_opt_real
                       fig6, ax6 = plt.subplots()
                       n6, bins6, patches6 = ax6.hist(L_int_Opt/1e43, bins=np.arange(0.1, 1 + 0.05, 0.05), facecolor='steelblue', density=True, alpha=0.4, label="Optimal Integrated\nLuminosity")
                       n7, bins7, patches7 = ax6.hist(L_mes17, bins=np.arange(0.1, 1 + 0.05, 0.05), color='red', histtype='step', density=True, label="Mesured Integrated\nLuminosity")
                       ax6.text(0.1, 4, testo1, bbox=dict(facecolor='none', edgecolor='steelblue', boxstyle='round,pad=0.5'))
                       ax6.text(0.1, 3.5, testo2, bbox=dict(facecolor='none', edgecolor='red', boxstyle='round,pad=0.5'))
                       ax6.set_ylabel("Normalized Frequencies")
                       ax6.set_xlabel("Integrated Luminosity [fb^-1]")
                       ax6.set_title("2017")
                       ax6.legend(loc='upper right')
                       plt.savefig('1_Immagini/2017/L_Int-Data_L_Int_Distributions.pdf') 
                       plt.show()         
                    
                      
                       print("________________________________________________________________")
                       print("______Comparison between total optimal luminosity_______________")
                       print("________and total luminosity evaluated from data________________")
                       print("____________________and the measured one________________________")
                       print("|")
                       print("| L_tot_opt=", L_tot1/1e43, "[fb^-1]")
                       print("|")
                       print("| L_tot_eval=", L_tot2/1e43, "[fb^-1]")
                       print("|")
                       print("| L_tot_data=", L_tot3, "[fb^-1]")
                       print("|")
                       print("________________________________________________________________")  

                   #Print the Average turn around time of the year and its respective optimal fill time and total luminosity
                   elif control3 == "7":
                       print("________________________________________________________________")
                       print("_____________________Average Turn Around Time___________________")
                       print("|")
                       print("| <ta>=", t_a17/3600, "[h]")
                       print("|")
                       print("________________________Optimal Fill Time_______________________")
                       print("|")
                       print("| t_opt=", t_opt17/3600, "[h]")
                       print("|")
                       print("_________________Optimal Integrated Luminosity__________________")
                       print("|")
                       print("| L_int=", L_int17/1e43, "[fb^-1]")
                       print("|")
                       print("________________________________________________________________")
                       print("_________________Optimal Total Luminosity_______________________")
                       print("|")
                       print("| L_opt=", L_tot17/1e43, "[fb^-1]")
                       print("|")
                       print("________________________________________________________________")
                       
                   #Exit
                   elif control3 == "b"or"q":
                       CONTROL3=False
                   else:
                       print("You did not enter a valid command!") 
            
            #2018 Preliminary Analysis
            elif control2 == "3":
               #Defining the turn around times in seconds
               Max = np.max(data_18_sec)
               Min  = np.min(data_18_sec)
               t_a = np.linspace(Min, Max, 100) #turn around time in seconds

               #Defining the turn around times in hours
               max = np.max(data_18)
               min  = np.min(data_18)
               tta = np.linspace(min, max, 100) #turn around time in hours  
               
               #model parameters
               n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps= lo.Parameters2018()
               
               #sorted data array
               ttta = lo.selection_sort(data_ta18_sec)
               
               #evaluating optimized quantities
               t_opt, L_tot_opt, L_int_opt, t_a18, t_opt18, L_int18, L_tot18 = lo.L_optimal_18(t_a)
               t_opt_real,  L_opt_real, L_int_Opt, t_a18, t_opt18, L_int18, L_tot18 = lo.L_optimal_18(ttta)
               L_inst_real, L_int_real, L_tot_real= lo.Model_L18(data_tf18_sec, data_ta18_sec) 
               Inst_x, Integ_x, Tot_x= lo.Model_L18(t_opt18, t_a18)
               
               CONTROL3=True
               while CONTROL3:
                   control3 = input("##################################################################################################################################\n"\
                       "-Enter 1: Print the model parameters;\n"\
                       "-Enter 2: Evaluate the Peak Luminosity;\n"\
                       "-Enter 3: Comparison between Optimal Fill Times and Fill Time Data;\n"\
                       "-Enter 4: Comparison between Optimal Integrated Luminosity and Integrated Luminosity Data;\n"\
                       "-Enter 5: Distribution of Optimal Fill Times and Fill Time Data;\n"\
                       "-Enter 6: Distribution of Optimal Integrated Luminosity and Integrated Luminosity data;\n" 
                       "-Enter 7: Print the Average turn around time of the year and its respective optimal fill time and total luminosity;\n"    
                       "-Enter b to go back to the previous menu.\n"   )
                   
                   #Print the model parameters
                   if control3 == "1":
                        print("______________________________________________________________________________________________")
                        print("|_____________________________________Model Parameters________________________________________")
                        print("| Bunch Population                                         n_i =",n_i/1e11,"e+11")
                        print("| Value of the beta-function at the collision point [m]    B_* =",B_s)
                        print("| RMS normalized transverse emittance               [m]    E_* =",E_s)
                        print("| Relativistic beta factor                                 B_r =",B_r)
                        print("| Relativistic gamma factor                                G_r =", G_r)
                        print("| Half Crossing Angle                               [rad]  T_hc =", T_hc)
                        print("| Longitudinal RMS dimension                        [m]    S_z =", S_z)
                        print("| Tranverse RMS dimension                           [m]    S_s =", S_s)
                        print("| Interaction cross section                         [rad]  S_int =", S_int)
                        print("| Revolution Frequency                              [Hz]   f_rev =", f_rev)
                        print("| Reduction in volume overlap factor                [m]    Fe =", Fe," ")
                        print("| HL Collision points                                      n_c =", n_c)
                        print("| Colliding Bunches                                        k_b =", k_b)
                        print("|                                                          \u039E =", Xi)
                        print("|                                                          \u03B5 =", Eps)
                        print("| Total Physics Time                                [s]    T_ph =", T_ph)
                        print("_____________________________________________________________________________________________")
                        
                   #Evaluate the Peak Luminosity
                   elif control3 == "2":
                       t_ff = 0
                       L_inst, L_int, L_tot = lo.Model_L18(t_ff, t_a)
                       print("________________________________________________________________")
                       print("_______________________Peak Luminosity__________________________")
                       print("|                                                               ")
                       print("| L_model(0)) = ", L_inst*1e-4, "[cm^-2 s^-1]")
                       print("| L_LHC(0) = 2.1e+34 cm^-2 s^-1")
                       print("________________________________________________________________")
                       
                   #Comparison between Optimal Fill Times and Fill Time Data
                   elif control3 == "3":                       
                       plt.close("all")
                       fig1, ax1 = plt.subplots()
                       ax1.plot(tta, t_opt/3600, 'k-', label="Optimal Fill Times")
                       ax1.plot([1.7,20], [13.775100322919943 , 13.775100322919943 ], 'r-', label="Optimal Fill Time for <ta>")
                       ax1.plot(data_ta18, data_tf18, 'b+', label= "2018 Data")
                       ax1.set_title('Comparison between Optimal Fill Times and 2018 Fill Time Data')
                       ax1.set_xlabel('Turn Around Times [h]')
                       ax1.set_ylabel('Fill Times [h]')
                       ax1.legend(loc='best')
                       plt.savefig('1_Immagini/2018/Optimal_tf-Data_tf.pdf')
                       plt.show()
                       
                   #Comparison between Optimal Integrated Luminosity and Integrated Luminosity Data
                   elif control3 == "4":                       
                       plt.close("all")
                       fig, ax = plt.subplots()
                       ax.plot(ttta/3600, L_int_Opt/1e43, 'k-', label="Optimal Luminosity")
                       ax.plot([2.3,23], [Integ_x/1e43, Integ_x/1e43], 'r-', label="Optimal Luminosity for <ta>")
                       ax.plot(data_ta18, L_mes18, 'b+', label="2018 Data")
                       ax.set_ylabel('Integrated Luminosity [fb^-1]')
                       ax.set_xlabel('Turn Around Times [h]')
                       ax.set_title('Comparison between Optimal Integrated Luminosity\nand 2018 Integrated Luminosity Data')
                       ax.legend(loc='best')
                       plt.savefig('1_Immagini/2018/Optimal_L_Int-Data_L_Int.pdf')
                       plt.show()
                       
                   #Distribution of Optimal Fill Times and Fill Time Data
                   elif control3 == "5":
                       c=N_i*n_c*Xi*S_int
                       a = 0.33245703
                       b = 1.50000008*3600
                       d = 0.16421911/3600
                       n = 0.27151320
                       f = lambda x: (a/(np.power((x-b), n)))*np.exp(-d*x)
                       val, err=integrate.quad(f, Min, Max)
                       g = lambda t_opt: ((a*2*c*t_opt)/(np.power((c*(np.power(t_opt, 2))-b), n)))*np.exp(-d*c*(np.power(t_opt, 2)))
                       a = 0.33245703/val 
                       
                       plt.close("all")
                       ###### HISTOGRAMS OF the distribution of t_opt_real and L_opt_real
                       fig5,  ax5 = plt.subplots(figsize=(5,4))
                       n5, bins5, patches5 = ax5.hist(t_opt_real/3600, bins=15, alpha=0.4, density=True, label="Optimal Fill Times")
                       ax5.plot(t_opt/3600, g(t_opt)*3600, "r-", label='Theorical Distribution of t_opt')
                       ax5.legend(loc='best')
                       ax5.set_xlabel("Optimal Fill Times [h]")
                       ax5.set_ylabel("Normalized Frequencies")
                       ax5.set_title("2018")
                       plt.savefig('1_Immagini/2018/t_opt_theory-distribution.pdf')
                       
                       T_ph_opt= np.sum(t_opt_real)
                       T_ph_opt = T_ph_opt/3600
                       T_ph_real= np.sum(data_tf18_sec)
                       T_ph_real=  T_ph_real/3600
                       print("________________________________________________________________")
                       print("____________________Total Physics time__________________________")
                       print("|                                                               ")
                       print("| T_ph_opt = ", T_ph_opt/24, "days")
                       print("| T_ph_data =", T_ph_real/24, "days")
                       print("|")
                       print("________________________________________________________________")

                       fig4,  ax4 = plt.subplots(figsize=(4,4))
                       n4, bins4, patches4 = ax4.hist(t_opt_real/3600, bins=np.arange(1, 29 + 2, 2), histtype='step', alpha=0.2, fill=True, label="Optimal Fill Times")
                       n4, bins4, patches4 = ax4.hist(data_tf18_sec/3600, bins=np.arange(1, 29 + 2, 2), histtype='step', fill=False, label="2018 Data")
                       ax4.legend(loc='best')
                       ax4.set_xlabel("Fill Times [h]")
                       ax4.set_ylabel("Number of Fills")
                       plt.savefig('1_Immagini/2018/t_opt_Distributions.pdf')
                       plt.show()
                       
                   #Distribution of Optimal Integrated Luminosity and Integrated Luminosity data
                   elif control3 == "6":
                       L_tot1 = np.sum(L_int_Opt)
                       L_tot2 = np.sum(L_int_real)
                       L_tot3 = np.sum(L_mes18)
                       testo1 = 'L_tot_opt={:.2f} [fb^-1]'.format(L_tot1/1e43)
                       testo2 = 'L_tot_data={:.2f} [fb^-1]'.format(L_tot3)
                       plt.close("all")
                       ###### HISTOGRAMS OF the distribution of t_opt_real and L_opt_real
                       fig6, ax6 = plt.subplots()
                       n6, bins6, patches6 = ax6.hist(L_int_Opt/1e43, bins=np.arange(0.1, 1 + 0.05, 0.05), facecolor='steelblue', density=True, alpha=0.4, label="Optimal Integrated\nLuminosity")
                       n7, bins7, patches7 = ax6.hist(L_mes18, bins=np.arange(0.1, 1 + 0.05, 0.05), color='red', histtype='step', density=True, label="Measured Integrated\nLuminosity")
                       ax6.text(0.1, 4.5, testo1, bbox=dict(facecolor='none', edgecolor='steelblue', boxstyle='round,pad=0.5'))
                       ax6.text(0.1, 4, testo2, bbox=dict(facecolor='none', edgecolor='red', boxstyle='round,pad=0.5'))
                       ax6.set_ylabel("Normalized Frequencies")
                       ax6.set_xlabel("Integrated Luminosity [fb^-1]")
                       ax6.set_title("2018")
                       ax6.legend(loc='upper right')
                       plt.savefig('1_Immagini/2018/L_Int-Data_L_Int_Distributions.pdf') 
                       plt.show() 
                       

                       print("________________________________________________________________")
                       print("______Comparison between total optimal luminosity_______________")
                       print("________and total luminosity evaluated from data________________")
                       print("____________________and the measured one________________________")
                       print("|")
                       print("| L_tot_opt=", L_tot1/1e43, "[fb^-1]")
                       print("|")
                       print("| L_tot_eval=", L_tot2/1e43, "[fb^-1]")
                       print("|")
                       print("| L_tot_data=", L_tot3, "[fb^-1]")
                       print("|")
                       print("________________________________________________________________")                  

                   #Print the Average turn around time of the year and its respective optimal fill time and total luminosity
                   elif control3 == "7":
                       print("________________________________________________________________")
                       print("_____________________Average Turn Around Time___________________")
                       print("|")
                       print("| <ta>=", t_a18/3600, "[h]")
                       print("|")
                       print("________________________Optimal Fill Time_______________________")
                       print("|")
                       print("| t_opt=", t_opt18/3600, "[h]")
                       print("|")
                       print("_________________Optimal Integrated Luminosity__________________")
                       print("|")
                       print("| L_int=", L_int18/1e43, "[fb^-1]")
                       print("|")
                       print("________________________________________________________________")
                       print("_________________Optimal Total Luminosity_______________________")
                       print("|")
                       print("| L_opt=", L_tot18/1e43, "[fb^-1]")
                       print("|")
                       print("________________________________________________________________")
                       
                   #Exit
                   elif control3 == "b"or"q":
                       CONTROL3=False
                   else:
                       print("You did not enter a valid command!") 
            
            #2016 Vincolated Analysis           
            elif control2 == '4':
               
               #model parameters 
               n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps= lo.Parameters2016()
               
               #sorted data array
               ttta = lo.selection_sort(data_ta16_sec)
               
               #evaluating optimized quantities 
               t_opt_real,  L_opt_real, L_int_Opt, t_a16, t_opt16, L_int16, L_tot16= lo.L_optimal_16(ttta)
               L_inst_real, L_int_real, L_tot_real= lo.Model_L16(data_tf16_sec, data_ta16_sec) 
               Inst_x, Integ_x, Tot_x= lo.Model_L16(t_opt16, t_a16)
               
               T_ph_data=np.sum(data_ta16_sec)+np.sum(data_tf16_sec)
               N = T_ph_data/(t_a16+t_opt16)
               T_ph_opt=N*(t_opt16+t_a16)
               L_opt_tot = N*Integ_x
               L_real_tot = np.sum(L_int_real)
               print("________________________________________________________________")
               print("______________Total physics Time and Total Luminosity___________")
               print("_____________________statistical evaluation_____________________")
               print("|                                                               ")             
               print("|   T_ph_data=",T_ph_data/(3600*24), "T_ph_opt=", T_ph_opt/(3600*24) )
               print("|   N_data=", len(data_ta16_sec), "N_opt=", N)
               print("|")
               print("________________________________________________________________")
               print("________________Total luminosity evaluated with models__________")
               print("|")
               print("|   L_real_tot_eval=", L_real_tot/1e43, "[fb^-1]", "L_opt_tot_eval=", L_opt_tot/1e43, "[fb^-1]")
               print("________________________________________________________________")
               print("________________Total luminosity from data extrapolation________")
               print("|")
               print("| L_tot_data_mes=", np.sum(L_mes16), "[fb^-1]")
               print("________________________________________________________________")
               new_t_opt=[]
               new_ta=[]
               new_t_opt=np.array(new_t_opt)
               new_ta=np.array(new_ta)
               for i in data_ta16_sec:
                   if np.sum(new_ta)+np.sum(new_t_opt)<T_ph_data:
                       new_t_opt=np.append(new_t_opt, t_opt16)
                       new_ta=np.append(new_ta, i)
                
               L_tot_real_ta=len(new_t_opt)*L_int16 
               T_ph_check=np.sum(new_ta)+np.sum(new_t_opt)
               print("________________________________________________________________")
               print("______________Total physics Time and Total Luminosity___________")
               print("_______________________with occurred tta________________________")
               print("|")
               print("|  T_ph=", T_ph_check/(3600*24), "N=", len(new_ta))
               print("|   L_tot_real_ta=", L_tot_real_ta/1e43, "[fb^-1]")
               print("________________________________________________________________")
               
               print(np.sum)    
                       
            #2017 Vincolated Analysis
            elif control2 == '5':
               n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps= lo.Parameters2017()
               ttta = lo.selection_sort(data_ta17_sec)
               t_opt_real,  L_opt_real, L_int_Opt, t_a17, t_opt17, L_int17, L_tot17= lo.L_optimal_17(ttta)
               L_inst_real, L_int_real, L_tot_real= lo.Model_L17(data_tf17_sec, data_ta17_sec) 
               Inst_x, Integ_x, Tot_x= lo.Model_L17(t_opt17, t_a17)
               
               T_ph_data=np.sum(data_ta17_sec)+np.sum(data_tf17_sec)
               N = T_ph_data/(t_a17+t_opt17)
               T_ph_opt=N*(t_opt17+t_a17)
               L_opt_tot = N*L_int17
               L_real_tot = np.sum(L_int_real)
               print("________________________________________________________________")
               print("______________Total physics Time and Total Luminosity___________")
               print("_____________________statistical evaluation_____________________")
               print("|                                                               ")            
               print("|   T_ph_data=",T_ph_data/(3600*24), "T_ph_opt=", T_ph_opt/(3600*24) )
               print("|   N_data=", len(data_ta17_sec), "N_opt=", N)
               print("|")
               print("________________________________________________________________")
               print("________________Total luminosity evaluated with models__________")
               print("|")
               print("|   L_real_tot_eval=", L_real_tot/1e43, "[fb^-1]", "L_opt_tot_eval=", L_opt_tot/1e43, "[fb^-1]")
               print("________________________________________________________________")
               print("________________Total luminosity from data extrapolation________")
               print("|")
               print("|   L_tot_data_mes=", np.sum(L_mes17), "[fb^-1]")
               print("________________________________________________________________")
               new_t_opt=[]
               new_ta=[]
               new_t_opt=np.array(new_t_opt)
               new_ta=np.array(new_ta)
               for i in data_ta17_sec:
                   if np.sum(new_ta)+np.sum(new_t_opt)<T_ph_data:
                       new_t_opt=np.append(new_t_opt, t_opt16)
                       new_ta=np.append(new_ta, i)
                
               L_tot_real_ta=len(new_t_opt)*L_int17
               T_ph_check=np.sum(new_ta)+np.sum(new_t_opt)
               print("________________________________________________________________")
               print("______________Total physics Time and Total Luminosity___________")
               print("_______________________with occurred tta________________________")
               print("|")
               print("|  T_ph=", T_ph_check/(3600*24), "N=", len(new_ta))
               print("|   L_tot_real_ta=", L_tot_real_ta/1e43, "[fb^-1]")
               print("________________________________________________________________")
            
            #2018 Vincolated Analysis
            elif control2 == '6':    
               n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps= lo.Parameters2018()
               ttta = lo.selection_sort(data_ta18_sec)
               t_opt_real,  L_opt_real, L_int_Opt, t_a18, t_opt18, L_int18, L_tot18= lo.L_optimal_18(ttta)
               L_inst_real, L_int_real, L_tot_real= lo.Model_L18(data_tf18_sec, data_ta18_sec) 
               Inst_x, Integ_x, Tot_x= lo.Model_L18(t_opt18, t_a18)
               
               T_ph_data=np.sum(data_ta18_sec)+np.sum(data_tf18_sec)
               N = T_ph_data/(t_a18+t_opt18)
               T_ph_opt=N*(t_opt18+t_a18)
               L_opt_tot = N*L_int18
               L_real_tot = np.sum(L_int_real)
               print("________________________________________________________________")
               print("______________Total physics Time and Total Luminosity___________")
               print("_____________________statistical evaluation_____________________")
               print("|                                                               ")             
               print("|   T_ph_data=",T_ph_data/(3600*24), "T_ph_opt=", T_ph_opt/(3600*24) )
               print("|   N_data=", len(data_ta18_sec), "N_opt=", N)
               print("|")
               print("________________________________________________________________")
               print("________________Total luminosity evaluated with models__________")
               print("|")
               print("|   L_real_tot_eval=", L_real_tot/1e43, "[fb^-1]", "L_opt_tot_eval=", L_opt_tot/1e43, "[fb^-1]")
               print("________________________________________________________________")
               print("________________Total luminosity from data extrapolation________")
               print("|")
               print("|   L_tot_data_mes=", np.sum(L_mes18), "[fb^-1]")
               print("________________________________________________________________")
               new_t_opt=[]
               new_ta=[]
               new_t_opt=np.array(new_t_opt)
               new_ta=np.array(new_ta)
               for i in data_ta18_sec:
                   if np.sum(new_ta)+np.sum(new_t_opt)<T_ph_data:
                       new_t_opt=np.append(new_t_opt, t_opt18)
                       new_ta=np.append(new_ta, i)
                
               L_tot_real_ta=len(new_t_opt)*L_int18 
               T_ph_check=np.sum(new_ta)+np.sum(new_t_opt)
               print("________________________________________________________________")
               print("______________Total physics Time and Total Luminosity___________")
               print("_______________________with occurred tta________________________")
               print("|")
               print("|  T_ph=", T_ph_check/(3600*24), "N=", len(new_ta))
               print("|   L_tot_real_ta=", L_tot_real_ta/1e43, "[fb^-1]")
               print("________________________________________________________________")             
            
            #Exit
            elif control2 == "b"or"q":
                CONTROL2=False
            else:
               print("You did not enter a valid command!")
     
     #New Run                         
     elif control1 == "2":
         
        #function to determine if a generic input is a number 
        def isanumber(a):
            """Function that evaluate if a generic input is a number or not."""
            try:
                float(a)
                bool_a = True
            except:
                bool_a = False

            return bool_a

        #loading measured data
        data_ta16, data_tf16, data_ta17, data_tf17, data_ta18, data_tf18 =ld.loadFill()
        alpha=0.05 #95% confidence level for the kolmogorov smirnov test

        #function to determine the optimal fill time
        def t_Opt(t_a):
            """Function that evaluates the optimal fill time.

            Args:
                t_a (any): (average) turn around time 

            Returns:
                any: [description]
            """
            return np.sqrt(t_a)/(np.sqrt(N_i*n_c*Xi*S_int)) 

        CONTROLz=True
        while CONTROLz:
            #choose the previous year
            previous_year=input('Insert the  previous year: ')

            #Evaluating the old sample model parameters for the evaluation of the optimal fill time 
            #and the old sample
            if previous_year=='2016':
                old_sample = np.array(data_ta16)
                old_sample = old_sample.flatten('F')
            elif previous_year=='end' or previous_year=='q':
                 break  
            elif previous_year=='2017':    
                old_sample = np.array(data_ta17)
                old_sample = old_sample.flatten('F')
            elif previous_year=='2018': 
                    old_sample = np.array(data_ta18)
                    old_sample = old_sample.flatten('F')
            elif previous_year!='2016' and previous_year!='2017' and previous_year!='2018': 
                    print('Remember to compile the excel file present in the program directory\ncalled Previous_Year_Data.xlsx where in a column called ta_stat\nthere must be the statistics sample for the old turnarund times\nand in another column, average_ta, the average of the turnaround times.\nHere should also be present the machine parameters for the current year.')
                        
                    file = pd.read_excel(r'Previous_Year_Data.xlsx', sheet_name='Parameters')
                    #n_i
                    df1 = pd.DataFrame(file, columns=['n_i']).dropna()
                    list1 = df1.values.tolist()
                    n_i=float(list1[0][0])
                    #B_s
                    df3 = pd.DataFrame(file, columns=['B_s']).dropna()
                    list3 = df3.values.tolist()
                    B_s=float(list3[0][0])
                    #E_s
                    df4 = pd.DataFrame(file, columns=['E_s']).dropna()
                    list4 = df4.values.tolist()
                    E_s=float(list4[0][0])
                    #B_r
                    df5 = pd.DataFrame(file, columns=['B_r']).dropna()
                    list5 = df5.values.tolist()
                    B_r=float(list5[0][0])
                    #G_r
                    df6 = pd.DataFrame(file, columns=['G_r']).dropna()
                    list6 = df6.values.tolist()
                    G_r=float(list6[0][0])
                    #T_hc
                    df7 = pd.DataFrame(file, columns=['T_hc']).dropna()
                    list7 = df7.values.tolist()
                    T_hc=float(list7[0][0])
                    #S_z
                    df8 = pd.DataFrame(file, columns=['S_z']).dropna()
                    list8 = df8.values.tolist()
                    S_z=float(list8[0][0])
                    #S_int
                    df9 = pd.DataFrame(file, columns=['S_int']).dropna()
                    list9 = df9.values.tolist()
                    S_int=float(list9[0][0])
                    #f_rev
                    df10 = pd.DataFrame(file, columns=['f_rev']).dropna()
                    list10 = df10.values.tolist()
                    f_rev=float(list10[0][0])
                    #n_c
                    df11 = pd.DataFrame(file, columns=['n_c']).dropna()
                    list11 = df11.values.tolist()
                    n_c=float(list11[0][0])
                    #k_b
                    df12 = pd.DataFrame(file, columns=['k_b']).dropna()
                    list12 = df12.values.tolist()
                    k_b=float(list12[0][0])
                    #T_ph
                    df13 = pd.DataFrame(file, columns=['T_ph']).dropna()
                    list13 = df13.values.tolist()
                    T_ph=float(list13[0][0])
                        
                    #intensity
                    N_i = k_b*n_i
            
                    #Definition of sigma* - Tranverse RMS dimension
                    S_s = np.sqrt((B_s*E_s)/(B_r*G_r))

                    #Definition of F(theta_c, sigma_z, sigma*) -  that accounts for the reduction in volume overlap
                    # between the colliding bunches due to the presence of a crossing angle
                    Fe = 1/(np.sqrt(1+(((T_hc*S_z)/(S_s))**2)))

                    #Definition of Xi
                    Xi = ((G_r*f_rev)/(4*np.pi*E_s*B_s*k_b))*Fe
            
                    # Epsilon Definition
                    Eps = (S_int*n_c*Xi)/f_rev
                    
                    
            #choosing the current year        
            current_year=input('Insert the  current year: ')

            #defining the output excel file 
            writer=pd.ExcelWriter('Outputs_file.xlsx')

            #evaluating the model parameters for the current year
            if current_year=='2016':
                n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps = lo.Parameters2016()
            elif current_year=='end' or current_year=='q':
                CONTROLz=False      
            elif current_year=='2017'or current_year=='linespace' or current_year=='single':    
                n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps = lo.Parameters2017()
            elif current_year=='2018': 
                n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps = lo.Parameters2018()
            elif current_year!='2016' and previous_year!='2017' and previous_year!='2018': 
                    #Importing datas from the excel file
                    data = pd.read_excel(r'Current_Year_Data.xlsx')
                    df1 = pd.DataFrame(data, columns=['ta_stat'])
                    old_sample_list = df1.values.tolist()
                    old_sample = df1.to_numpy()
                    old_sample= old_sample.flatten('F')
                    av_tta = pd.DataFrame(data, columns=['average_ta']).values.tolist()
                    
                    #importing the machine parameters
                    param=input("Insert the current year machine parameters with a file or one by one.\n\
                                 Enter 'file' in the first case and 'tip' in the latter.\n")
                    if param == 'file':
                        print('Remember to compile the excel file present in the program directory\ncalled Current_Year_Data.xlsx where in a column called ta_stat\nthere must be the statistics sample for the old turnarund times\n\
                                and in another column, average_ta, the average of the turnaround times.\n\
                                Here should also be present the machine parameters for the current year.')
                        
                        file2 = pd.read_excel(r'Current_Year_Data.xlsx', sheet_name='Parameters')
                        #n_i
                        df1 = pd.DataFrame(file2, columns=['n_i']).dropna()
                        list1 = df1.values.tolist()
                        n_i=float(list1[0][0])
                        #B_s
                        df3 = pd.DataFrame(file2, columns=['B_s']).dropna()
                        list3 = df3.values.tolist()
                        B_s=float(list3[0][0])
                        #E_s
                        df4 = pd.DataFrame(file2, columns=['E_s']).dropna()
                        list4 = df4.values.tolist()
                        E_s=float(list4[0][0])
                        #B_r
                        df5 = pd.DataFrame(file2, columns=['B_r']).dropna()
                        list5 = df5.values.tolist()
                        B_r=float(list5[0][0])
                        #G_r
                        df6 = pd.DataFrame(file2, columns=['G_r']).dropna()
                        list6 = df6.values.tolist()
                        G_r=float(list6[0][0])
                        #T_hc
                        df7 = pd.DataFrame(file2, columns=['T_hc']).dropna()
                        list7 = df7.values.tolist()
                        print(list7)
                        T_hc=float(list7[0][0])
                        #S_z
                        df8 = pd.DataFrame(file2, columns=['S_z']).dropna()
                        list8 = df8.values.tolist()
                        S_z=float(list8[0][0])
                        #S_int
                        df9 = pd.DataFrame(file2, columns=['S_int']).dropna()
                        list9 = df9.values.tolist()
                        S_int=float(list9[0][0])
                        #f_rev
                        df10 = pd.DataFrame(file2, columns=['f_rev']).dropna()
                        list10 = df10.values.tolist()
                        f_rev=float(list10[0][0])
                        #n_c
                        df11 = pd.DataFrame(file2, columns=['n_c']).dropna()
                        list11 = df11.values.tolist()
                        n_c=float(list11[0][0])
                        #k_b
                        df12 = pd.DataFrame(file2, columns=['k_b']).dropna()
                        list12 = df12.values.tolist()
                        k_b=float(list12[0][0])
                        #T_ph
                        df13 = pd.DataFrame(file2, columns=['T_ph']).dropna()
                        list13 = df13.values.tolist()
                        T_ph=float(list13[0][0])
                        
                        
                        #intensity
                        N_i = k_b*n_i
            
                        #Definition of sigma* - Tranverse RMS dimension
                        S_s = np.sqrt((B_s*E_s)/(B_r*G_r))

                        #Definition of F(theta_c, sigma_z, sigma*) -  that accounts for the reduction in volume overlap
                        # between the colliding bunches due to the presence of a crossing angle
                        Fe = 1/(np.sqrt(1+(((T_hc*S_z)/(S_s))**2)))

                        #Definition of Xi
                        Xi = ((G_r*f_rev)/(4*np.pi*E_s*B_s*k_b))*Fe
            
                        # Epsilon Definition
                        Eps = (S_int*n_c*Xi)/f_rev
                    
                    #insert manually the machine parameters    
                    elif param == 'tip':
                        print('Insert the current year machine parameters: ')
                        c = 299792458 #[m/s]
                        C_LHC = 26658.883 #[m]
                        m = 0.938 #[GeV]
                        n_c = 2 #Number of collision points
                        control=True
                        while control:
                            try:
                                #insert parameters
                                E =float(input("Total Energy[GeV] E= "))
                                G_r = E/m #gamma_r
                                B_r = np.sqrt(1-(1/(np.power(G_r, 2)))) #beta_r
                                S_z = float(input("Longitudinal RMS Dimension [m] \u03C3_z= "))#[m] sigma_z longitudinal RMS dimension 
                                S_int = float(input("Interaction cross section [m^2] \u03C3_int= ")) #[m^2]interaction cross section
                                n_i = float(input("Intensity of the beam n_i= ")) #Intensity of the beam (only time dependent beam parameter)
                                B_s = float(input("value of the beta-function at the collision point [m] \u03B2*= ")) #[m] beta* -  value of the beta-function at the collision point
                                E_s = float(input("RMS normalized transverse emittance [m] \u03B5*= ")) #[m] epsilon* - RMS normalized transverse emittance
                                T_hc = float(input("Half crossing angle [rad] \u03B8_hc= ")) #[rad] theta_hc - half crossing angle 
                                k_b = float(input("Number of colliding bunches k_b= "))  #number of colliding bunches
                                N_i = k_b*n_i #Intensity of the beam (only time dependent beam parameter)
                                T_ph = int(input("Total Physics Time [days] T_ph= "))*24*3600 #[s] total physics time 
                                t_a = float(input("Turn Around Time [s] t_ta= "))
                                control = False
                            except ValueError:
                                print('You insert an invalid quantity!')
                                continue
                            
                        #Definition of sigma* - Tranverse RMS dimension
                        S_s = np.sqrt((B_s*E_s)/(B_r*G_r))

                        #Definition of F(theta_c, sigma_z, sigma*) -  that accounts for the reduction in volume overlap
                        # between the colliding bunches due to the presence of a crossing angle
                        Fe = 1/(np.sqrt(1+(((T_hc*S_z)/(S_s))**2)))


                        #Definition of the revolution frequency
                        f_rev = c/(C_LHC) #Hz

                        #Definition of Xi
                        Xi = ((G_r*f_rev)/(4*np.pi*E_s*B_s*k_b))*Fe
                
                        #Epsilon Definition
                        Eps = (S_int*n_c*Xi)/f_rev
                    



            #printing the current year machine parameters
            print("___________New Sample Model:", current_year, "_______________________________________________")
            print("|_____________________________________Model Parameters________________________________________")
            print("| Bunch Population                                         n_i =",n_i/1e11,"e+11")
            print("| Value of the beta-function at the collision point [m]    B_* =",B_s)
            print("| RMS normalized transverse emittance               [m]    E_* =",E_s)
            print("| Relativistic beta factor                                 B_r =",B_r)
            print("| Relativistic gamma factor                                G_r =", G_r)
            print("| Half Crossing Angle                               [rad]  T_hc =", T_hc)
            print("| Longitudinal RMS dimension                        [m]    S_z =", S_z)
            print("| Tranverse RMS dimension                           [m]    S_s =", S_s)
            print("| Interaction cross section                         [rad]  S_int =", S_int)
            print("| Revolution Frequency                              [Hz]   f_rev =", f_rev)
            print("| Reduction in volume overlap factor                [m]    Fe =", Fe)
            print("| HL Collision points                                      n_c =", n_c)
            print("| Colliding Bunches                                        k_b =", k_b)
            print("|                                                          \u039E =", Xi)
            print("|                                                          \u03B5 =", Eps)
            print("| Total Physics Time                                [s]    T_ph =", T_ph)
            print("_____________________________________________________________________________________________") 




            #Defining the new sample of turnaround times

            #2016 data
            if current_year=='2016':
                    new_sample=np.array(data_ta16)
                    len_new_sample=len(new_sample) #number of fills
                    new_sample=new_sample.flatten('F')
                    new_Sample = []
                    average_Ta=[]
                    optimal_fill=[]
                    weights = []
                    new_Sample=np.array(new_Sample)
                    Weights=np.array(weights)
                    Opt_Fill=np.array(optimal_fill)
                    average_Ta=np.array(average_Ta)
                    
                    #K-S test on the first 30 turnaround 
                    for i in new_sample:
                        new_Sample = np.append(new_Sample, i)
                        if len(new_Sample)<= 30:
                            test = ks_2samp(new_Sample, old_sample)
                            p_value=test[1]
                            Weights = np.append(Weights, p_value)
                            #evaluation of the everage turn around time and the corrispondent optimal fill time
                            average_ta= float((np.sum(old_sample)+ np.sum(np.multiply(Weights, new_Sample)))/(1*len(old_sample)+np.sum(Weights)))
                            average_Ta=np.append(average_Ta, average_ta)
                            t_opt = t_Opt(average_ta*3600)
                            Opt_Fill=np.append(Opt_Fill, t_opt)
                        elif len(new_Sample)> 30:
                            #evaluation of the everage turn around time and the corrispondent optimal fill time
                            average_ta= float((np.sum(new_sample))/(len(new_sample)))
                            average_Ta=np.append(average_Ta, average_ta)
                            t_opt = t_Opt(average_ta*3600)
                            Opt_Fill=np.append(Opt_Fill, t_opt)
                        
                        #defining the output dataframe 
                        df_Opt_Fill=pd.DataFrame(Opt_Fill/3600, columns=['Opt t_fill'])
                        df_old_sample=pd.DataFrame(old_sample, columns=['{} sample'.format(previous_year)])  
                        df_new_sample=pd.DataFrame(new_Sample, columns=['{} sample'.format(current_year)]) 
                        
                        #writing on the output file
                        df_Opt_Fill.to_excel(writer, sheet_name='Optimal fill times')
                        df_new_sample.to_excel(writer, sheet_name='New sample')
                        df_old_sample.to_excel(writer, sheet_name='Old Sample')  
                                    
                        if (np.sum(Opt_Fill/3600)+np.sum(new_Sample))*3600>= T_ph:
                            print('You reach the nominal physics time!')
                            break 
            
            #Exit
            elif current_year=='end' or current_year=='q':
                 break
                 
            #2017 data                
            elif current_year=='2017':    
                new_sample=np.array(data_ta17)
                len_new_sample=len(new_sample) #number of fills
                new_sample=new_sample.flatten('F')
                new_Sample = []
                average_Ta=[]
                optimal_fill=[]
                weights = []
                new_Sample=np.array(new_Sample)
                Weights=np.array(weights)
                Opt_Fill=np.array(optimal_fill)
                average_Ta=np.array(average_Ta)
                
                #K-S test on the first 30 turnaround 
                for i in new_sample:
                    new_Sample = np.append(new_Sample, i)
                    if len(new_Sample)<= 30:
                        test = ks_2samp(new_Sample, old_sample)
                        p_value=test[1]
                        Weights = np.append(Weights, p_value)
                        #evaluation of the everage turn around time and the corrispondent optimal fill time
                        average_ta= float((np.sum(old_sample)+ np.sum(np.multiply(Weights, new_Sample)))/(1*len(old_sample)+np.sum(Weights)))
                        average_Ta=np.append(average_Ta, average_ta)
                        t_opt = t_Opt(average_ta*3600)
                        Opt_Fill=np.append(Opt_Fill, t_opt)
                    elif len(new_Sample)> 30:
                        #evaluation of the everage turn around time and the corrispondent optimal fill time
                        average_ta= float((np.sum(new_sample))/(len(new_sample)))
                        average_Ta=np.append(average_Ta, average_ta)
                        t_opt = t_Opt(average_ta*3600)
                        Opt_Fill=np.append(Opt_Fill, t_opt)
                    
                    #defining the output dataframe 
                    df_Opt_Fill=pd.DataFrame(Opt_Fill/3600, columns=['Opt t_fill'])
                    df_old_sample=pd.DataFrame(old_sample, columns=['{} sample'.format(previous_year)])  
                    df_new_sample=pd.DataFrame(new_Sample, columns=['{} sample'.format(current_year)]) 
                        
                    #writing on the output file
                    df_Opt_Fill.to_excel(writer, sheet_name='Optimal fill times')
                    df_new_sample.to_excel(writer, sheet_name='New sample')
                    df_old_sample.to_excel(writer, sheet_name='Old Sample')   
                                    
                                
                    if (np.sum(Opt_Fill/3600)+np.sum(new_Sample))*3600>= T_ph:
                        print('You reach the nominal physics time!')
                        break 
            
            # generated data         
            elif current_year=='generated2016':  
                max = np.max(data_ta16)
                min  = np.min(data_ta16)
                tta = np.linspace(min, max, 1000)    
                new_sample=tta
                len_new_sample=len(new_sample) #number of fills
                new_Sample = []
                average_Ta=[]
                optimal_fill=[]
                weights = []
                new_Sample=np.array(new_Sample)
                Weights=np.array(weights)
                Opt_Fill=np.array(optimal_fill)
                average_Ta=np.array(average_Ta)
                
                #K-S test on the first 30 turnaround 
                for i in new_sample:
                    new_Sample = np.append(new_Sample, i)
                    if len(new_Sample)<= 30:
                        test = ks_2samp(new_Sample, old_sample)
                        p_value=test[1]
                        Weights = np.append(Weights, p_value)
                        #evaluation of the everage turn around time and the corrispondent optimal fill time
                        average_ta= float((np.sum(old_sample)+ np.sum(np.multiply(Weights, new_Sample)))/(1*len(old_sample)+np.sum(Weights)))
                        average_Ta=np.append(average_Ta, average_ta)
                        t_opt = t_Opt(average_ta*3600)
                        Opt_Fill=np.append(Opt_Fill, t_opt)
                    elif len(new_Sample)> 30:
                        #evaluation of the everage turn around time and the corrispondent optimal fill time
                        average_ta= float((np.sum(new_sample))/(len(new_sample)))
                        average_Ta=np.append(average_Ta, average_ta)
                        t_opt = t_Opt(average_ta*3600)
                        Opt_Fill=np.append(Opt_Fill, t_opt)
                        
                    #defining the output dataframe 
                    df_Opt_Fill=pd.DataFrame(Opt_Fill/3600, columns=['Opt t_fill'])
                    df_old_sample=pd.DataFrame(old_sample, columns=['{} sample'.format(previous_year)])  
                    df_new_sample=pd.DataFrame(new_Sample, columns=['{} sample'.format(current_year)]) 
                        
                    #writing on the output file
                    df_Opt_Fill.to_excel(writer, sheet_name='Optimal fill times')
                    df_new_sample.to_excel(writer, sheet_name='New sample')
                    df_old_sample.to_excel(writer, sheet_name='Old Sample') 
                                
                    if (np.sum(Opt_Fill/3600)+np.sum(new_Sample))*3600>= T_ph:
                        print('You reach the nominal physics time!')
                        break 

            elif current_year=='generated2017':  
                max = np.max(data_ta17)
                min  = np.min(data_ta17)
                tta = np.linspace(min, max, 1000)    
                new_sample=tta
                len_new_sample=len(new_sample) #number of fills
                new_Sample = []
                average_Ta=[]
                optimal_fill=[]
                weights = []
                new_Sample=np.array(new_Sample)
                Weights=np.array(weights)
                Opt_Fill=np.array(optimal_fill)
                average_Ta=np.array(average_Ta)
                
                #K-S test on the first 30 turnaround 
                for i in new_sample:
                    new_Sample = np.append(new_Sample, i)
                    if len(new_Sample)<= 30:
                        test = ks_2samp(new_Sample, old_sample)
                        p_value=test[1]
                        Weights = np.append(Weights, p_value)
                        #evaluation of the everage turn around time and the corrispondent optimal fill time
                        average_ta= float((np.sum(old_sample)+ np.sum(np.multiply(Weights, new_Sample)))/(1*len(old_sample)+np.sum(Weights)))
                        average_Ta=np.append(average_Ta, average_ta)
                        t_opt = t_Opt(average_ta*3600)
                        Opt_Fill=np.append(Opt_Fill, t_opt)
                    elif len(new_Sample)> 30:
                        #evaluation of the everage turn around time and the corrispondent optimal fill time
                        average_ta= float((np.sum(new_sample))/(len(new_sample)))
                        average_Ta=np.append(average_Ta, average_ta)
                        t_opt = t_Opt(average_ta*3600)
                        Opt_Fill=np.append(Opt_Fill, t_opt)
                        
                    #defining the output dataframe 
                    df_Opt_Fill=pd.DataFrame(Opt_Fill/3600, columns=['Opt t_fill'])
                    df_old_sample=pd.DataFrame(old_sample, columns=['{} sample'.format(previous_year)])  
                    df_new_sample=pd.DataFrame(new_Sample, columns=['{} sample'.format(current_year)]) 
                        
                    #writing on the output file
                    df_Opt_Fill.to_excel(writer, sheet_name='Optimal fill times')
                    df_new_sample.to_excel(writer, sheet_name='New sample')
                    df_old_sample.to_excel(writer, sheet_name='Old Sample') 
                                
                    if (np.sum(Opt_Fill/3600)+np.sum(new_Sample))*3600>= T_ph:
                        print('You reach the nominal physics time!')
                        break 

            elif current_year=='generated2018':  
                max = np.max(data_ta18)
                min  = np.min(data_ta18)
                tta = np.linspace(min, max, 1000)    
                new_sample=tta
                len_new_sample=len(new_sample) #number of fills
                new_Sample = []
                average_Ta=[]
                optimal_fill=[]
                weights = []
                new_Sample=np.array(new_Sample)
                Weights=np.array(weights)
                Opt_Fill=np.array(optimal_fill)
                average_Ta=np.array(average_Ta)
                
                #K-S test on the first 30 turnaround 
                for i in new_sample:
                    new_Sample = np.append(new_Sample, i)
                    if len(new_Sample)<= 30:
                        test = ks_2samp(new_Sample, old_sample)
                        p_value=test[1]
                        Weights = np.append(Weights, p_value)
                        #evaluation of the everage turn around time and the corrispondent optimal fill time
                        average_ta= float((np.sum(old_sample)+ np.sum(np.multiply(Weights, new_Sample)))/(1*len(old_sample)+np.sum(Weights)))
                        average_Ta=np.append(average_Ta, average_ta)
                        t_opt = t_Opt(average_ta*3600)
                        Opt_Fill=np.append(Opt_Fill, t_opt)
                    elif len(new_Sample)> 30:
                        #evaluation of the everage turn around time and the corrispondent optimal fill time
                        average_ta= float((np.sum(new_sample))/(len(new_sample)))
                        average_Ta=np.append(average_Ta, average_ta)
                        t_opt = t_Opt(average_ta*3600)
                        Opt_Fill=np.append(Opt_Fill, t_opt)
                        
                    #defining the output dataframe 
                    df_Opt_Fill=pd.DataFrame(Opt_Fill/3600, columns=['Opt t_fill'])
                    df_old_sample=pd.DataFrame(old_sample, columns=['{} sample'.format(previous_year)])  
                    df_new_sample=pd.DataFrame(new_Sample, columns=['{} sample'.format(current_year)]) 
                        
                    #writing on the output file
                    df_Opt_Fill.to_excel(writer, sheet_name='Optimal fill times')
                    df_new_sample.to_excel(writer, sheet_name='New sample')
                    df_old_sample.to_excel(writer, sheet_name='Old Sample') 
                                
                    if (np.sum(Opt_Fill/3600)+np.sum(new_Sample))*3600>= T_ph:
                        print('You reach the nominal physics time!')
                        break 
                    
                            
            #2018 data        
            elif current_year=='2018': 
                new_sample=np.array(data_ta18)
                len_new_sample=len(new_sample) #number of fills
                new_sample=new_sample.flatten('F')
                new_Sample = []
                average_Ta=[]
                optimal_fill=[]
                weights = []
                new_Sample=np.array(new_Sample)
                Weights=np.array(weights)
                Opt_Fill=np.array(optimal_fill)
                average_Ta=np.array(average_Ta)
                
                #K-S test on the first 30 turnaround 
                for i in new_sample:
                    new_Sample = np.append(new_Sample, i)
                    if len(new_Sample)<= 30:
                        test = ks_2samp(new_Sample, old_sample)
                        p_value=test[1]
                        Weights = np.append(Weights, p_value)
                        #evaluation of the everage turn around time and the corrispondent optimal fill time
                        average_ta= float((np.sum(old_sample)+ np.sum(np.multiply(Weights, new_Sample)))/(1*len(old_sample)+np.sum(Weights)))
                        average_Ta=np.append(average_Ta, average_ta)
                        t_opt = t_Opt(average_ta*3600)
                        Opt_Fill=np.append(Opt_Fill, t_opt)
                    elif len(new_Sample)> 30:
                        #evaluation of the everage turn around time and the corrispondent optimal fill time
                        average_ta= float((np.sum(new_sample))/(len(new_sample)))
                        average_Ta=np.append(average_Ta, average_ta)
                        t_opt = t_Opt(average_ta*3600)
                        Opt_Fill=np.append(Opt_Fill, t_opt)
                        
                    #defining the output dataframe 
                    df_Opt_Fill=pd.DataFrame(Opt_Fill/3600, columns=['Opt t_fill'])
                    df_old_sample=pd.DataFrame(old_sample, columns=['{} sample'.format(previous_year)])  
                    df_new_sample=pd.DataFrame(new_Sample, columns=['{} sample'.format(current_year)]) 
                        
                    #writing on the output file
                    df_Opt_Fill.to_excel(writer, sheet_name='Optimal fill times')
                    df_new_sample.to_excel(writer, sheet_name='New sample')
                    df_old_sample.to_excel(writer, sheet_name='Old Sample')  
                                        
                                
                    if (np.sum(Opt_Fill/3600)+np.sum(new_Sample))*3600>= T_ph:
                        print('You reach the nominal physics time!')
                        break 

            #new year data - online data entry            
            elif current_year!='2016' and previous_year!='2017' and previous_year!='2018': 
                CONTROL=True
                new_Sample = []
                average_Ta=[]
                optimal_fill=[]
                weights = []
                new_Sample=np.array(new_Sample)
                Weights=np.array(weights)
                Opt_Fill=np.array(optimal_fill)
                average_Ta=np.array(average_Ta)
                
                while CONTROL:
                    new_ta = input("Enter the new value of the turn around time,\n\
                                    when the year ends or when you want to exit print end.\nNew Turnaround time= ")
                    if new_ta == "end" or new_ta == "End" or new_ta=='q':
                        CONTROL=False
                    elif (new_ta != "end" or new_ta != "End" or new_ta!='q') and (isanumber(new_ta) == True):
                        new_Sample = np.append(new_Sample, float(new_ta))
                        #evaluation of the everage turn around time and the corrispondent optimal fill time
                        if len(new_Sample)<= 30:
                            test = ks_2samp(new_Sample, old_sample)
                            p_value=test[1]
                            Weights = np.append(Weights, p_value)
                            
                            # case in which the samples are compatible
                            if p_value>alpha:
                                print('The new sample seems to be distributed according to the previous year distribution!')
                                print('The new sample is:', new_Sample)
                                choice= input('Do you want to discharge the last value? yes/no\n')    
                                if choice == 'yes':
                                    new_Sample=np.delete(new_Sample, len(new_Sample)-1)
                                    Weights = np.delete(Weights, len(Weights)-1)
                                    print('The new sample is:', new_Sample )
                                    average_ta= float((np.sum(old_sample)+ np.sum(np.multiply(Weights, new_Sample)))/(1*len(old_sample)+np.sum(Weights))) 
                                    t_opt = t_Opt(average_ta*3600)            
                                    Opt_Fill=np.append(Opt_Fill, t_opt)
                                    print("The average t_a is:", average_ta)
                                    print("The optimal fill time t_opt is:", t_opt/3600) 
                                    
                                    #defining the output dataframe 
                                    df_Opt_Fill=pd.DataFrame(Opt_Fill/3600, columns=['Opt t_fill'])
                                    df_old_sample=pd.DataFrame(old_sample, columns=['{} sample'.format(previous_year)])  
                                    df_new_sample=pd.DataFrame(new_Sample, columns=['{} sample'.format(current_year)]) 
                        
                                    #writing on the output file
                                    df_Opt_Fill.to_excel(writer, sheet_name='Optimal fill times')
                                    df_new_sample.to_excel(writer, sheet_name='New sample')
                                    df_old_sample.to_excel(writer, sheet_name='Old Sample') 
                                
                                    if (np.sum(Opt_Fill/3600)+np.sum(new_Sample))*3600>= T_ph:
                                        print('You reach the nominal physics time!')
                                        break 
                                    
                                if choice == 'no':
                                    print('The new sample is still:', new_Sample )
                                    average_ta= float((np.sum(old_sample)+ np.sum(np.multiply(Weights, new_Sample)))/(1*len(old_sample)+np.sum(Weights))) 
                                    t_opt = t_Opt(average_ta*3600)
                                    Opt_Fill=np.append(Opt_Fill, t_opt)
                                    print("The average t_a is:", average_ta)
                                    print("The optimal fill time t_opt is:", t_opt/3600) 
                                    
                                    #defining the output dataframe 
                                    df_Opt_Fill=pd.DataFrame(Opt_Fill/3600, columns=['Opt t_fill'])
                                    df_old_sample=pd.DataFrame(old_sample, columns=['{} sample'.format(previous_year)])  
                                    df_new_sample=pd.DataFrame(new_Sample, columns=['{} sample'.format(current_year)]) 
                        
                                    #writing on the output file
                                    df_Opt_Fill.to_excel(writer, sheet_name='Optimal fill times')
                                    df_new_sample.to_excel(writer, sheet_name='New sample')
                                    df_old_sample.to_excel(writer, sheet_name='Old Sample') 
                                    
                                    if (np.sum(Opt_Fill/3600)+np.sum(new_Sample))*3600>= T_ph:
                                        print('You reach the nominal physics time!')
                                        break
                                
                                    continue
                            # case in which the samples are not compatible
                            elif p_value<alpha:
                                print('It seems that the new sample does not follow the old sample distribution!')
                                print('The new sample is:', new_Sample )
                                choice= input('Do you want to discharge the last value? yes/no\n')    
                                if choice == 'yes' or choice=='q' or choice=='end':
                                    new_Sample=np.delete(new_Sample, len(new_Sample)-1)
                                    Weights = np.delete(Weights, len(Weights)-1)
                                    print('The new sample is:', new_Sample )
                                    average_ta=float(np.sum(new_Sample)/len(new_Sample))
                                    t_opt = t_Opt(average_ta*3600)
                                    Opt_Fill=np.append(Opt_Fill, t_opt)
                                    print("The average t_a is:", average_ta)
                                    print("The optimal fill time t_opt is:", t_opt/3600)
                                    
                                    #defining the output dataframe 
                                    df_Opt_Fill=pd.DataFrame(Opt_Fill/3600, columns=['Opt t_fill'])
                                    df_old_sample=pd.DataFrame(old_sample, columns=['{} sample'.format(previous_year)])  
                                    df_new_sample=pd.DataFrame(new_Sample, columns=['{} sample'.format(current_year)]) 
                        
                                    #writing on the output file
                                    df_Opt_Fill.to_excel(writer, sheet_name='Optimal fill times')
                                    df_new_sample.to_excel(writer, sheet_name='New sample')
                                    df_old_sample.to_excel(writer, sheet_name='Old Sample') 
                                    
                                    if (np.sum(Opt_Fill/3600)+np.sum(new_Sample))*3600>= T_ph:
                                        print('You reach the nominal physics time!')
                                        break   
                                                    
                                if choice == 'no':
                                    print('The new sample is still:', new_Sample )
                                    average_ta=float(np.sum(new_Sample)/len(new_Sample) )
                                    t_opt = t_Opt(average_ta*3600)
                                    Opt_Fill=np.append(Opt_Fill, t_opt)
                                    print("The average t_a is:", average_ta)
                                    print("The optimal fill time t_opt is:", t_opt)
                                    
                                    #defining the output dataframe 
                                    df_Opt_Fill=pd.DataFrame(Opt_Fill/3600, columns=['Opt t_fill'])
                                    df_old_sample=pd.DataFrame(old_sample, columns=['{} sample'.format(previous_year)])  
                                    df_new_sample=pd.DataFrame(new_Sample, columns=['{} sample'.format(current_year)]) 
                        
                                    #writing on the output file
                                    df_Opt_Fill.to_excel(writer, sheet_name='Optimal fill times')
                                    df_new_sample.to_excel(writer, sheet_name='New sample')
                                    df_old_sample.to_excel(writer, sheet_name='Old Sample') 
                                    
                                    if (np.sum(Opt_Fill/3600)+np.sum(new_Sample))*3600>= T_ph:
                                        print('You reach the nominal physics time!')
                                        break
                                                    
                                    continue  
                        # for enough bigger sample evaluating its own average value of turn around time and 
                        #its matching fill time 
                        elif len(new_Sample)>30:  
                            average_ta= float((np.sum(new_Sample))/(len(new_Sample)))
                            average_Ta=np.append(average_Ta, average_ta)
                            t_opt = t_Opt(average_ta*3600)
                            Opt_Fill=np.append(Opt_Fill, t_opt)
                            print("The average t_a is:", average_ta)
                            print("The optimal fill time t_opt is:", t_opt/3600)      
            
                        
                    elif (new_ta != "end" or new_ta != "End" or new_ta!='q') and isanumber(new_ta) == False:   
                        print('You not enter a valid input! Here what you enter: ', new_ta)


            #Verify the Total Physics time
            print('Nominal Physics time', T_ph/3600, '[h]')
            print('The total physics time is now:', np.sum(Opt_Fill/3600)+np.sum(new_Sample), '[h]')
            writer.save()
            
            end=input("Do you want to go back to the previous menu? yes/no\n")
            if end=='yes' or end=='end' or end=='q':
                CONTROLz=False
            elif end=='no':
                continue    
    
     #Turn Around statistical analysis    
     elif control1 == "3":
        CONTROL2 = True
        while CONTROL2:
         
         #load data and bins   
         data_tot, dataTot, array_tot = ld.TotalDataSet(data_16, data_17, data_18)   
         data_tot_A, data_tot_B, data_tot_C, array_totA, array_totB, array_totC = ld.PartialDataSets(data_16, data_17, data_18)
         bi16, bi17, bi18, biT, biTA, biTB, biTC = cvb.CreateBins(array16, array17, array18, array_tot, array_totA, array_totB, array_totC)
         
         #load fitting models
         NEl_mod = mod.N_ExpLaw_Model() 
         NTpl_mod = mod.N_TruncLaw_Model()
         control2 = input("###########################################\n"\
                      "-Enter 1: 2016;\n-Enter 2: 2017;\n-Enter 3: 2018;\n-Enter 4: Total Run 2;\n-Enter 5: Kolmogorov-Smirnov Test;\n-Enter b to go back to the previous menu.\n")
         
         #2016
         if control2=="1":
             #fit parameters
              A = np.min(array16)
              B = np.max(array16)
              
              #fitting models
              def HighNormalized_PowerLaw (x, Off, n):
                  return ((1-n)/(((B-Off)**(1-n))-(A-Off)**(n-1)))*(1/((x-Off)**n))
              HnPl_mod = Model(HighNormalized_PowerLaw) 
              CONTROL3 = True
              while CONTROL3==True:
                control3=input("###########################################\n\
                    -Enter 1: Fits with 3 different fitting functions;\n\
                    -Enter 2: Truncated Power Law Fit;\n\
                    -Enter 3: Average value of the turn around time;\n\
                    -Enter b to go back to the previous menu.\n")
                
                #Fits with 3 different fitting functions
                if control3=="1":
                    
                      #Fits
                      fig1, ax1 = plt.subplots()
                      n1, bins1, patches1 = ax1.hist(data_16, bins=bi16, facecolor='darkseagreen', alpha=0.5, density=True, label="2016 Data")
                      binscenters = np.array([0.5 * (bins1[i] + bins1[i+1]) for i in range(len(bins1)-1)])
                      x = binscenters
                      y = n1
                      hnpl_result = HnPl_mod.fit(y, x=x, Off=1.8, n=0.5) 
                      ax1.plot(binscenters, hnpl_result.best_fit, 'r--', label="Power law")
                      ax1.plot([], [], 'rx ', label='Reduced Chi-Square ={:.5f}'.format(hnpl_result.redchi))
                      NTpl_mod.set_param_hint('off', value=1.7, min=1.5, max=2.0)
                      NTpl_mod.set_param_hint('lam', value=0.2, min=0.1, max=0.4)
                      ntpl_result = NTpl_mod.fit(y, x=x, off=1.7, amp=1, lam=0.2, n=0.9) 
                      ax1.plot(binscenters, ntpl_result.best_fit, 'b--', label="Truncated Power law")
                      ax1.plot([], [], 'bx ', label='Reduced Chi-Square ={:.5f}'.format(ntpl_result.redchi))
                      nel_result = NEl_mod.fit(y, x=x, off=2, lam=0) 
                      ax1.plot(binscenters, nel_result.best_fit, 'k:', label="Exponential Law")
                      ax1.plot([], [], 'kx ', label='Reduced Chi-Square ={:.5f}'.format(nel_result.redchi))
                      plt.title('Turn Around Times 2016')
                      plt.xlabel('Turn Around Times [h]')
                      plt.ylabel('Normalized Frequency')
                      plt.legend(loc='best')
                      
                      #Residuals
                      fig2, (ax2A, ax2B, ax2C) = plt.subplots(3,1, sharex=True, sharey=True, figsize=(10, 5))
                      hnpl_result.plot_residuals(ax2A, datafmt='r*', yerr=None, data_kws=None, fit_kws=None, ax_kws=None, parse_complex='abs')
                      ax2A.set_ylabel("Residuals")
                      ax2A.set_title("Residuals for the Power Law - 2016")
                      ntpl_result.plot_residuals(ax2B, datafmt='b*', yerr=None, data_kws=None, fit_kws=None, ax_kws=None, parse_complex='abs')
                      ax2B.set_ylabel("Residuals")
                      ax2B.set_title("Residuals for the Truncated Power Law")
                      nel_result.plot_residuals(ax2C, datafmt='k*', yerr=None, data_kws=None, fit_kws=None, ax_kws=None, parse_complex='abs')
                      plt.xlabel("Bincenters")
                      ax2C.set_ylabel("Residuals")
                      ax2C.set_title("Residuals for the Exponential Law")
                      
                      #printing results
                      print('____________________2016_______________________________')
                      print('RESULTs OF THE POWER LAW FIT')
                      print('-------------------------------------------------------')
                      print(hnpl_result.fit_report())
                      print('_______________________________________________________')
                      print('_______________________________________________________')
                      print('_______________________________________________________')
                      print('RESULTs OF THE TRUNCATED POWER LAW FIT')
                      print('-------------------------------------------------------')
                      print(ntpl_result.fit_report())
                      print('_______________________________________________________')
                      print('_______________________________________________________')
                      print('_______________________________________________________')
                      print('RESULTs OF THE EXPONENTIAL FIT')
                      print('-------------------------------------------------------')
                      print(nel_result.fit_report())
                      print('_______________________________________________________')
                      print('_______________________________________________________')
                      plt.show()
                
                #Truncated Power Law Fit
                elif control3 == "2":
                    
                      #Truncated power law fit
                      plt.close("all")
                      fig1, ax1 = plt.subplots()
                      n1, bins1, patches1 = ax1.hist(data_16, bins=bi16, facecolor='darkseagreen', alpha=0.5, density=True, label="2016 Data")
                      binscenters = np.array([0.5 * (bins1[i] + bins1[i+1]) for i in range(len(bins1)-1)])
                      x = binscenters
                      y = n1
                      NTpl_mod.set_param_hint('off', value=1.7, min=1.5, max=2.0)
                      NTpl_mod.set_param_hint('lam', value=0.2, min=0.1, max=0.4)
                      ntpl_result = NTpl_mod.fit(y, x=x, off=1.7, amp=1, lam=0.2, n=0.9) 
                      ax1.plot(binscenters, ntpl_result.best_fit, 'b--', label="Truncated Power law")
                      ax1.plot([], [], 'b. ', label='Reduced Chi-Square ={:.5f}'.format(ntpl_result.redchi))
                      ax1.plot([], [], 'bx', label='Offset = {:.3f} +/- {:.3f}'.format(ntpl_result.params['off'].value, ntpl_result.params['off'].stderr ))
                      ax1.plot([], [], 'bx', label='Exponent = {:.3f} +/- {:.3f}'.format(ntpl_result.params['n'].value, ntpl_result.params['n'].stderr ))
                      ax1.plot([], [], 'bx', label='Amplitude = {:.3f} +/- {:.3f}'.format(ntpl_result.params['amp'].value, ntpl_result.params['amp'].stderr ))
                      ax1.plot([], [], 'bx', label='\u03BB = {:.3f} +/- {:.3f}'.format(ntpl_result.params['lam'].value, ntpl_result.params['lam'].stderr ))
                      plt.title('Turn Around Times 2016')
                      plt.xlabel('Turn Around Times [h]')
                      plt.ylabel('Normalized Frequency')
                      plt.legend(loc='best')
                      plt.show()
                
                #Average value of the turn around time
                elif control3 == "3":
                    
                      #### Test for the interval of the integrand - 2016
                      x1 = np.linspace(2, 150, 1000)
                      a1 = 0.37082720*x1*(np.power((x1-2),(-0.73454371)))*(np.exp(-0.10000000*x1))
                      plt.close("all")
                      fig1, ax1 = plt.subplots()
                      ax1.plot(x1, a1, "b-")
                      ax1.set_title("f_trunc16(t)")
                      #plt.show()
                      
                      #####2016 Epectation Value
                      def f_trunc16(t):
                         a = 0.37082720
                         b = 2
                         d = 0.10000000
                         n = 0.73454371
                         return t*((a/(np.power(t-b, n)))*np.exp(-d*t))
                      val16, err16 = integrate.quad(f_trunc16, 2, 250)
                      
                      #èrinting results
                      print("_____________________________________________________")
                      print("2016: E[t_ta]=",val16,"+-", err16)
                      print("_____________________________________________________")                     
                
                #Exit
                elif control3 == "b"or"q":
                     CONTROL3=False
         
         #2017                    
         elif control2=="2":
              #fit parameters
              A = np.min(array17)
              B = np.max(array17)
              
              #fitting models
              def HighNormalized_PowerLaw (x, Off, n):
                 return ((1-n)/(((B-Off)**(1-n))-(A-Off)**(n-1)))*(1/((x-Off)**n))
              HnPl_mod = Model(HighNormalized_PowerLaw) 
              CONTROL3 = True
              while CONTROL3==True:
                control3=input("###########################################\n\
                    -Enter 1: Fits with 3 different fitting functions;\n\
                    -Enter 2: Truncated Power Law Fit;\n\
                    -Enter 3: Average value of the turn around time;\n\
                    -Enter b to go back to the previous menu.\n")
                
                #Fits with 3 different fitting functions
                if control3=="1":
                  #Fits  
                  plt.close("all")
                  fig3, ax3 = plt.subplots()
                  n3, bins3, patches3 = ax3.hist(data_17, bins=bi17, facecolor='steelblue', alpha=0.4, density=True, label="2017 Data")
                  binscenters = np.array([0.5 * (bins3[i] + bins3[i+1]) for i in range(len(bins3)-1)])
                  x = binscenters
                  y = n3
                  hnpl_result = HnPl_mod.fit(y, x=x, Off=1.8, n=0.5) 
                  ax3.plot(binscenters, hnpl_result.best_fit, 'r--', label="Power law")
                  ax3.plot([], [], 'rx ', label='Reduced Chi-Square ={:.5f}'.format(hnpl_result.redchi))
                  NTpl_mod.set_param_hint('off', value=1.7, min=1.5, max=2.0)
                  NTpl_mod.set_param_hint('lam', value=0.2, min=0.1, max=0.4)
                  ntpl_result = NTpl_mod.fit(y, x=x, off=1.7, amp=0.1, lam=0.2, n=1) 
                  ax3.plot(binscenters, ntpl_result.best_fit, 'b--', label="Truncated Power law")
                  ax3.plot([], [], 'bx ', label='Reduced Chi-Square ={:.5f}'.format(ntpl_result.redchi))
                  nel_result = NEl_mod.fit(y, x=x, off=2, lam=0) 
                  plt.plot(binscenters, nel_result.best_fit, 'k:', label="Exponential Power law")
                  ax3.plot([], [], 'kx ', label='Reduced Chi-Square ={:.5f}'.format(nel_result.redchi))
                  plt.title('Turn Around Times 2017')
                  plt.xlabel('Turn Around Times [h]')
                  plt.ylabel('Normalized Frequency')
                  plt.legend(loc='best')
                  
                  #residuals
                  fig4, (ax4A, ax4B, ax4C) = plt.subplots(3,1, sharex=True, sharey=True, figsize=(10, 5))
                  hnpl_result.plot_residuals(ax4A, datafmt='r*', yerr=None, data_kws=None, fit_kws=None, ax_kws=None, parse_complex='abs')
                  ax4A.set_ylabel("Residuals")
                  ax4A.set_title("Residuals for the Power Law - 2017")
                  ntpl_result.plot_residuals(ax4B, datafmt='b*', yerr=None, data_kws=None, fit_kws=None, ax_kws=None, parse_complex='abs')
                  ax4B.set_ylabel("Residuals")
                  ax4B.set_title("Residuals for the Truncated Power Law")
                  nel_result.plot_residuals(ax4C, datafmt='k*', yerr=None, data_kws=None, fit_kws=None, ax_kws=None, parse_complex='abs')
                  plt.xlabel("Bincenters")
                  ax4C.set_ylabel("Residuals")
                  ax4C.set_title("Residuals for the Exponential Law")
                  
                  #printing results
                  print('__________________2017_________________________________')
                  print('RESULTs OF THE POWER LAW FIT')
                  print('-------------------------------------------------------')
                  print(hnpl_result.fit_report())
                  print('_______________________________________________________')
                  print('_______________________________________________________')
                  print('_______________________________________________________')
                  print('RESULTs OF THE TRUNCATED POWER LAW FIT')
                  print('-------------------------------------------------------')
                  print(ntpl_result.fit_report())
                  print('_______________________________________________________')
                  print('_______________________________________________________')
                  print('_______________________________________________________')
                  print('RESULTs OF THE EXPONENTIAL FIT')
                  print('-------------------------------------------------------')
                  print(nel_result.fit_report())
                  print('_______________________________________________________')
                  print('_______________________________________________________')
                  plt.show()

                #Truncated Power Law Fit
                elif control3 == "2":
                    
                  plt.close("all")
                  fig5, ax5 = plt.subplots()
                  n3, bins3, patches3 = ax5.hist(data_17, bins=bi17, facecolor='steelblue', alpha=0.4, density=True, label="2017 Data")
                  binscenters = np.array([0.5 * (bins3[i] + bins3[i+1]) for i in range(len(bins3)-1)])
                  x = binscenters
                  y = n3
                  NTpl_mod.set_param_hint('off', value=1.7, min=1.5, max=2.0)
                  NTpl_mod.set_param_hint('lam', value=0.2, min=0.1, max=0.4)
                  ntpl_result = NTpl_mod.fit(y, x=x, off=1.7, amp=0.1, lam=0.2, n=1) 
                  ax5.plot(binscenters, ntpl_result.best_fit, 'b--', label="Truncated Power law")
                  ax5.plot([], [], 'b. ', label='Reduced Chi-Square ={:.5f}'.format(ntpl_result.redchi))
                  ax5.plot([], [], 'bx', label='Offset = {:.3f} +/- {:.3f}'.format(ntpl_result.params['off'].value, ntpl_result.params['off'].stderr ))
                  ax5.plot([], [], 'bx', label='Exponent = {:.3f} +/- {:.3f}'.format(ntpl_result.params['n'].value, ntpl_result.params['n'].stderr ))
                  ax5.plot([], [], 'bx', label='Amplitude = {:.3f} +/- {:.3f}'.format(ntpl_result.params['amp'].value, ntpl_result.params['amp'].stderr ))
                  ax5.plot([], [], 'bx', label='\u03BB = {:.3f} +/- {:.3f}'.format(ntpl_result.params['lam'].value, ntpl_result.params['lam'].stderr ))
                  ax5.set_title('Turn Around Times 2017')
                  ax5.set_xlabel('Turn Around Times [h]')
                  ax5.set_ylabel('Normalized Frequency')
                  plt.legend(loc='best')
                  plt.show()
                  
                #Average value of the turn around time
                elif control3 == "3":
                      #### Test for the interval of the integrand - 2017
                      x2 = np.linspace(2, 150, 1000)
                      a2 = 0.29687182*x2*(np.power((x2-2),(-0.50729741)))*(np.exp(-0.10422466*x2))
                      fig2, ax2 = plt.subplots()
                      ax2.plot(x2, a2, "b-")
                      ax2.set_title("f_trunc17(t)")
                      #plt.show()
                      
                      #####2017 Epectation Value
                      def f_trunc17(t):
                         a = 0.29687182
                         b = 2
                         d = 0.10422466
                         n = 0.50729741
                         return t*((a/(np.power(t-b, n)))*np.exp(-d*t))
                      val17, err17 = integrate.quad(f_trunc17, 2, 250)
                      print("_____________________________________________________")
                      print("2017: E[t_ta]=",val17,"+-", err17)
                      print("_____________________________________________________")                    
                
                #Exit
                elif control3 == "b"or"q":
                     CONTROL3=False
         
         #2018            
         elif control2=="3":
             #fitting parameters
              A = np.min(array18)
              B = np.max(array18)
              
              #fitting models
              def HighNormalized_PowerLaw (x, Off, n):
                    return ((1-n)/(((B-Off)**(1-n))-(A-Off)**(n-1)))*(1/((x-Off)**n))
              HnPl_mod = Model(HighNormalized_PowerLaw)   
                          
              CONTROL3 = True
              while CONTROL3==True:
                control3=input("###########################################\n\
                    -Enter 1: Fits with 3 different fitting functions;\n\
                    -Enter 2: Truncated Power Law Fit;\n\
                    -Enter 3: Average value of the turn around time;\n\
                    -Enter b to go back to the previous menu.\n")
                
                #Fits with 3 different fitting functions
                if control3=="1":
                  #Fits  
                  plt.close("all")
                  fig6, ax6 = plt.subplots()
                  n3, bins3, patches3 = ax6.hist(data_18, bins=bi18, facecolor='r', alpha=0.2, density=True, label="2018 Data")
                  binscenters = np.array([0.5 * (bins3[i] + bins3[i+1]) for i in range(len(bins3)-1)])
                  x = binscenters
                  y = n3
                  hnpl_result = HnPl_mod.fit(y, x=x, Off=1, n=0.2) 
                  ax6.plot(binscenters, hnpl_result.best_fit, 'r--', label="Power law")
                  ax6.plot([], [], 'rx ', label='Reduced Chi-Square ={:.5f}'.format(hnpl_result.redchi))
                  #NTpl_mod.set_param_hint('off', value=1.0, min=0, max=2.0)
                  #ntpl_result = NTpl_mod.fit(y, x=x, off=1, amp=1, lam=0.2, n=0.9) 
                  NTpl_mod.set_param_hint('off', value=1.9, min=1.5, max=2.0)
                  NTpl_mod.set_param_hint('lam', value=0.2, min=0.1, max=0.4)
                  ntpl_result = NTpl_mod.fit(y, x=x, off=1.9, amp=0.5, lam=0.2, n=0.5) 
                  ax6.plot(binscenters, ntpl_result.best_fit, 'b--', label="Truncated Power law")
                  ax6.plot([], [], 'bx ', label='Reduced Chi-Square ={:.5f}'.format(ntpl_result.redchi))
                  nel_result = NEl_mod.fit(y, x=x, off=2, lam=0) 
                  ax6.plot(binscenters, nel_result.best_fit, 'k:', label="Exponential Power law")
                  ax6.plot([], [], 'kx ', label='Reduced Chi-Square ={:.5f}'.format(nel_result.redchi))
                  plt.title('Turn Around Times 2018')
                  plt.xlabel('Turn Around Times [h]')
                  plt.ylabel('Normalized Frequency')
                  plt.legend(loc='best')

                  #Residuals
                  fig7, (ax7A, ax7B, ax7C) = plt.subplots(3,1, sharex=True, sharey=True, figsize=(10, 5))
                  hnpl_result.plot_residuals(ax7A, datafmt='r*', yerr=None, data_kws=None, fit_kws=None, ax_kws=None, parse_complex='abs')
                  ax7A.set_ylabel("Residuals")
                  ax7A.set_title("Residuals for the Power Law - 2018")
                  ntpl_result.plot_residuals(ax7B, datafmt='b*', yerr=None, data_kws=None, fit_kws=None, ax_kws=None, parse_complex='abs')
                  ax7B.set_ylabel("Residuals")
                  ax7B.set_title("Residuals for the Truncated Power Law")
                  nel_result.plot_residuals(ax7C, datafmt='k*', yerr=None, data_kws=None, fit_kws=None, ax_kws=None, parse_complex='abs')
                  plt.xlabel("Bincenters")
                  ax7C.set_ylabel("Residuals")
                  ax7C.set_title("Residuals for the Exponential Law")

                  #plotting results
                  print('__________________2018_________________________________')
                  print('RESULTs OF THE POWER LAW FIT')
                  print('-------------------------------------------------------')
                  print(hnpl_result.fit_report())
                  print('_______________________________________________________')
                  print('_______________________________________________________')
                  print('_______________________________________________________')
                  print('RESULTs OF THE TRUNCATED POWER LAW FIT')
                  print('-------------------------------------------------------')
                  print(ntpl_result.fit_report())
                  print('_______________________________________________________')
                  print('_______________________________________________________')
                  print('_______________________________________________________')
                  print('RESULTs OF THE EXPONENTIAL FIT')
                  print('-------------------------------------------------------')
                  print(nel_result.fit_report())
                  print('_______________________________________________________')
                  print('_______________________________________________________')
                  plt.show()
                
                #Truncated Power Law Fit
                elif control3 == "2":
                    
                      plt.close("all")
                      fig8, ax8 = plt.subplots()
                      n3, bins3, patches3 = ax8.hist(data_18, bins=bi18, facecolor='r', alpha=0.2, density=True, label="2018 Data")
                      binscenters = np.array([0.5 * (bins3[i] + bins3[i+1]) for i in range(len(bins3)-1)])
                      x = binscenters
                      y = n3
                      NTpl_mod.set_param_hint('off', value=1.9, min=1.5, max=2.0)
                      NTpl_mod.set_param_hint('lam', value=0.2, min=0.1, max=0.4)
                      ntpl_result = NTpl_mod.fit(y, x=x, off=1.9, amp=0.5, lam=0.2, n=0.5) 
                      ax8.plot(binscenters, ntpl_result.best_fit, 'b--', label="Truncated Power law")
                      ax8.plot([], [], 'b. ', label='Reduced Chi-Square ={:.5f}'.format(ntpl_result.redchi))
                      ax8.plot([], [], 'bx', label='Offset = {:.3f} +/- {:.3f}'.format(ntpl_result.params['off'].value, ntpl_result.params['off'].stderr ))
                      ax8.plot([], [], 'bx', label='Exponent = {:.3f} +/- {:.3f}'.format(ntpl_result.params['n'].value, ntpl_result.params['n'].stderr ))
                      ax8.plot([], [], 'bx', label='Amplitude = {:.3f} +/- {:.3f}'.format(ntpl_result.params['amp'].value, ntpl_result.params['amp'].stderr ))
                      ax8.plot([], [], 'bx', label='\u03BB = {:.3f} +/- {:.3f}'.format(ntpl_result.params['lam'].value, ntpl_result.params['lam'].stderr ))
                      plt.title('Turn Around Times 2018')
                      plt.xlabel('Turn Around Times [h]')
                      plt.ylabel('Normalized Frequency')
                      plt.legend(loc='best')
                      plt.show()
                
                #Average value of the turn around time
                elif control3 == "3":
                      #### Test for the interval of the integrand - 2018
                      x3 = np.linspace(1.50000008, 150, 1000)
                      a3 = 0.316*x3*(np.power((x3-1.901),(-0.964)))*(np.exp(-0.017*x3))
                      fig3, ax3 = plt.subplots()
                      ax3.plot(x3, a3, "b-")
                      ax3.set_title("f_trunc18(t)")
                      # plt.show()
                      
                      #####2018 Epectation Value
                      def f_trunc18(t):
                         a = 0.33245703
                         b = 1.50000008
                         d = 0.16421911
                         n = 0.27151320
                         return t*((a/(np.power(t-b, n)))*np.exp(-d*t))
                      val18, err18 = integrate.quad(f_trunc18, 1.50000008, 250)
                      print("_____________________________________________________")
                      print("2018: E[t_ta]=",val18,"+-", err18)
                      print("_____________________________________________________")                     
         
         #Total Run             
         elif control2=="4":
              #fitting parameters
              A = np.min(array_tot)
              B = np.max(array_tot)
              
              #fitting models 
              def HighNormalized_PowerLaw (x, Off, n):
                 return ((1-n)/(((B-Off)**(1-n))-(A-Off)**(n-1)))*(1/((x-Off)**n))
              HnPl_mod = Model(HighNormalized_PowerLaw)
                            
              CONTROL3 = True
              while CONTROL3==True:
                control3=input("##############################################\n\
                    -Enter 1: Fits with 3 different fitting functions;\n\
                    -Enter 2: Truncated Power Law Fit;\n\
                    -Enter b to go back to the previous menu.\n")
                
                #Fits with 3 different fitting functions
                if control3=="1":
                      
                      #Fits
                      plt.close("all")
                      fig9, ax9 = plt.subplots()
                      n4, bins4, patches4 = ax9.hist(dataTot, bins=biT, facecolor='gold', alpha=0.5, density=True, label="Run 2 Data")
                      binscenters = np.array([0.5 * (bins4[i] + bins4[i+1]) for i in range(len(bins4)-1)])
                      x = binscenters
                      y = n4
                      hnpl_result = HnPl_mod.fit(y, x=x, Off=1.8, n=0.5) 
                      ax9.plot(binscenters, hnpl_result.best_fit, 'r--', label="Power law")
                      ax9.plot([], [], 'rx ', label='Reduced Chi-Square ={:.5f}'.format(hnpl_result.redchi))
                      NTpl_mod.set_param_hint('off', value=1.9, min=1.5, max=2.0)
                      ntpl_result = NTpl_mod.fit(y, x=x, off=1.9, amp=1, lam=0, n=0.9) 
                      ax9.plot(binscenters, ntpl_result.best_fit, 'b--', label="Truncated Power law")
                      ax9.plot([], [], 'bx ', label='Reduced Chi-Square ={:.5f}'.format(ntpl_result.redchi))
                      nel_result = NEl_mod.fit(y, x=x, off=2, lam=0) 
                      ax9.plot(binscenters, nel_result.best_fit, 'k:', label="Exponential law")
                      ax9.plot([], [], 'kx ', label='Reduced Chi-Square ={:.5f}'.format(nel_result.redchi))
                      plt.title('Turn Around Times Total')
                      plt.xlabel('Turn Around Times [h]')
                      plt.ylabel('Normalized Frequency')
                      plt.legend(loc='best')

                      #Residuals
                      fig10, (ax10A, ax10B, ax10C) = plt.subplots(3,1, sharex=True, sharey=True, figsize=(10, 5))
                      hnpl_result.plot_residuals(ax10A, datafmt='r*', yerr=None, data_kws=None, fit_kws=None, ax_kws=None, parse_complex='abs')
                      ax10A.set_ylabel("Residuals")
                      ax10A.set_title("Residuals for the Power Law - RUN 2")
                      ntpl_result.plot_residuals(ax10B, datafmt='b*', yerr=None, data_kws=None, fit_kws=None, ax_kws=None, parse_complex='abs')
                      ax10B.set_ylabel("Residuals")
                      ax10B.set_title("Residuals for the Truncated Power Law")
                      nel_result.plot_residuals(ax10C, datafmt='k*', yerr=None, data_kws=None, fit_kws=None, ax_kws=None, parse_complex='abs')
                      plt.xlabel("Bincenters")
                      ax10C.set_ylabel("Residuals")
                      ax10C.set_title("Residuals for the Exponential Law")

                      #printing results
                      print('__________________RUN 2_________________________________')
                      print('RESULTs OF THE POWER LAW FIT')
                      print('-------------------------------------------------------')
                      print(hnpl_result.fit_report())
                      print('_______________________________________________________')
                      print('_______________________________________________________')
                      print('_______________________________________________________')
                      print('RESULTs OF THE TRUNCATED POWER LAW FIT')
                      print('-------------------------------------------------------')
                      print(ntpl_result.fit_report())
                      print('_______________________________________________________')
                      print('_______________________________________________________')
                      print('_______________________________________________________')
                      print('RESULTs OF THE EXPONENTIAL FIT')
                      print('-------------------------------------------------------')
                      print(nel_result.fit_report())
                      print('_______________________________________________________')
                      print('_______________________________________________________')
                      plt.show()
                      
                #Truncated Power Law Fit      
                elif control3 == "2":
                    
                      plt.close("all")
                      fig11, ax11 = plt.subplots()
                      n4, bins4, patches4 = ax11.hist(dataTot, bins=biT, facecolor='gold', alpha=0.5, density=True, label="Run 2 Data")
                      binscenters = np.array([0.5 * (bins4[i] + bins4[i+1]) for i in range(len(bins4)-1)])
                      x = binscenters
                      y = n4
                      NTpl_mod.set_param_hint('off', value=1.9, min=1.5, max=2.0)
                      ntpl_result = NTpl_mod.fit(y, x=x, off=1.9, amp=1, lam=0, n=0.9) 
                      ax11.plot(binscenters, ntpl_result.best_fit, 'b--', label="Truncated Power law")
                      ax11.plot([], [], 'bx ', label='Reduced Chi-Square ={:.5f}'.format(ntpl_result.redchi))
                      ax11.plot([], [], 'bx', label='Offset = {:.3f} +/- {:.3f}'.format(ntpl_result.params['off'].value, ntpl_result.params['off'].stderr ))
                      ax11.plot([], [], 'bx', label='Exponent = {:.3f} +/- {:.3f}'.format(ntpl_result.params['n'].value, ntpl_result.params['n'].stderr ))
                      ax11.plot([], [], 'bx', label='Amplitude = {:.3f} +/- {:.3f}'.format(ntpl_result.params['amp'].value, ntpl_result.params['amp'].stderr ))
                      ax11.plot([], [], 'bx', label='\u03BB = {:.3f} +/- {:.3f}'.format(ntpl_result.params['lam'].value, ntpl_result.params['lam'].stderr ))
                      plt.title('Turn Around Times Total')
                      plt.xlabel('Turn Around Times [h]')
                      plt.ylabel('Normalized Frequency')
                      plt.legend(loc='best')
                      plt.show()                   
 
                #Exit
                elif control3 == "b"or"q":
                     CONTROL3=False
         
         #Kolmogorov Smirnov Test                
         elif control2=="5":
             
            #Performing the tests
            test1 = ks_2samp(data_16, data_17)
            test2 = ks_2samp(data_16, data_18)
            test3 = ks_2samp(data_17, data_18)

            print(test1)
            print(test2)
            print(test3)
            print('Lokking at the K-S results t is possible to say that\n\
                data coming from 2016 and 2017 seems to be distributed according\n\
                to the same distribution, as data coming from 2017 and 2018.\n\
                However, data for 2016 and 2018 samples it is not possible to accept\n\
                the hyothesis of equal distributions. This last conslusion can be linked\n\
                to the consistent difference in the number of objects in the two samples.\n\
                Furthermore, taking into account the results as a whole, it is possible to state,\n\
                albeit with a certain limit of uncertainty, that the turn-around times of different\n\
                years are distributed following the same distribution.')  
            
         #Exit                
         elif control2 == "b"or"q":
                CONTROL2=False        
     
     #Exit
     elif control1 == "q":
        print("      Thank you for choosing OPTIMIZE!    ")
        tprint("Have a Good Day!",font="white_square")
        break
     else:
        print("You did not enter a valid command!")  
        