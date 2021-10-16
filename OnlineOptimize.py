#################################################################################
#
# @Giulia Faletti
#Strategy to chose online the optimal fill time 
#
#################################################################################
import numpy as np
from numpy.ma import array
import pandas as pd
from scipy.stats import ks_2samp
import LuminosityOptimization as lo
import LoadData as ld



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

#choose the previous year
previous_year=input('Insert the  previous year:')

#Evaluating the old sample model parameters for the evaluation of the optimal fill time 
#and the old sample
if previous_year=='2016':
    old_sample = np.array(data_ta16)
    old_sample = old_sample.flatten('F')
elif previous_year=='2017':    
    old_sample = np.array(data_ta17)
    old_sample = old_sample.flatten('F')
elif previous_year=='2018': 
        old_sample = np.array(data_ta18)
        old_sample = old_sample.flatten('F')
elif previous_year!='2016' and previous_year!='2017' and previous_year!='2018': 
        print('Remember to compile the excel file present in the program directory\n\
            called Previous_Year_Data.xlsx where in a column called ta_stat\n\
            there must be the statistics sample for the old turnarund times\n\
                and in another column, average_ta, the average of the turnaround times.\n\
                    Here should also be present the machine parameters for the current year.')
            
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
current_year=input('Insert the  current year:')

#defining the output excel file 
writer=pd.ExcelWriter('Outputs_file.xlsx')

#evaluating the model parameters for the current year
if current_year=='2016':
    n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps = lo.Parameters2016()
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
        param=input("Insert the current year machine parameters with a file or one by one. Enter 'file' in the first case and 'tip' in the latter.")
        if param == 'file':
            print('Remember to compile the excel file present in the program directory\n\
            called Current_Year_Data.xlsx where in a column called ta_stat\n\
            there must be the statistics sample for the old turnarund times\n\
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
            print('Insert the current year machine parameters:')
            c = 299792458 #[m/s]
            C_LHC = 26658.883 #[m]
            m = 0.938 #[GeV]
            n_c = 2 #Number of collision points
            control=True
            while control:
                try:
                    #insert parameters
                    E =float(input("Total Energy[GeV] E="))
                    G_r = E/m #gamma_r
                    B_r = np.sqrt(1-(1/(np.power(G_r, 2)))) #beta_r
                    S_z = float(input("Longitudinal RMS Dimension [m] \u03C3_z="))#[m] sigma_z longitudinal RMS dimension 
                    S_int = float(input("Interaction cross section [m^2] \u03C3_int=")) #[m^2]interaction cross section
                    n_i = float(input("Intensity of the beam n_i= ")) #Intensity of the beam (only time dependent beam parameter)
                    B_s = float(input("value of the beta-function at the collision point [m] \u03B2*=")) #[m] beta* -  value of the beta-function at the collision point
                    E_s = float(input("RMS normalized transverse emittance [m] \u03B5*=")) #[m] epsilon* - RMS normalized transverse emittance
                    T_hc = float(input("Half crossing angle [rad] \u03B8_hc=")) #[rad] theta_hc - half crossing angle 
                    k_b = float(input("Number of colliding bunches k_b="))  #number of colliding bunches
                    N_i = k_b*n_i #Intensity of the beam (only time dependent beam parameter)
                    T_ph = int(input("Total Physics Time [days] T_ph="))*24*3600 #[s] total physics time 
                    t_a = float(input("Turn Around Time [s] t_ta="))
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
        new_ta = input("Enter the new value of the turn around time, when the year ends or you want to exit print end.\nNew Turnaround time= ")
        if new_ta == "end" or new_ta == "End":
            CONTROL=False
        elif (new_ta != "end" or new_ta != "End") and (isanumber(new_ta) == True):
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
                    choice= input('You want to discharge the last value? yes/no\n')    
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
                    choice= input('You want to discharge the last value? yes/no\n')    
                    if choice == 'yes':
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
  
            
        elif (new_ta != "end" or new_ta != "End") and isanumber(new_ta) == False:   
            print('You not enter a valid input! Here what you enter:', new_ta)



print(new_sample, new_Sample, Opt_Fill/3600)
print(len(new_sample), len(new_Sample))

#Verify the Total Physics time
print('Nominal Physics time', T_ph/3600, '[h]')
print('The total physics time is now:', np.sum(Opt_Fill/3600)+np.sum(new_Sample), '[h]')
writer.save()






        
        

