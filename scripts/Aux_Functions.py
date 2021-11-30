
import numpy as np

def time_string(t):
    if t>3600*24*365*1000000:
        t_string = 't = %.2f Myrs'%(t/3600/24/365/1000/1000)
    elif t>3600*24*365*1000:
        t_string = 't = %.2f kyrs'%(t/3600/24/365/1000)
    elif t>3600*24*365:
        t_string = 't = %.1f yrs'%(t/3600/24/365)
    elif t>3600*24:
        t_string = 't = %.1f days'%(t/3600/24)
    elif t>3600:
        t_string = 't = %.1f hrs'%(t/3600)
    elif t>60:
        t_string = 't = %.1f mins'%(t/60)
    elif t>=1:
        t_string = 't = %.1f secs'%(t)
    else:
        t_string = 't = %.1f msecs'%(t*1000)
    return t_string



def palette(min_color,max_color):
    from matplotlib.colors import LinearSegmentedColormap, to_hex

    cm=LinearSegmentedColormap.from_list('plt',[min_color,max_color])

    idx=np.linspace(0,255,num=256,dtype=int)
    plt=[to_hex(cm(i)) for i in idx]
    return plt



def Equil_P(T_x,C=True):
    
    if C:
        T_x+=273.15
    
    
    A0 = -1.94138504464560e05
    A1 = +3.31018213397926e03   
    A2 = -2.25540264493806e01
    A3 = +7.67559117787059e-2  
    A4 = -1.30465829788791e-4
    A5 = +8.86065316687571e-8

    B0 = -4.38921173434628e01
    B1 = +7.76302133739303e-1   
    B2 = -7.27291427030502e-3
    B3 = +3.85413985900724e-5   
    B4 = -1.03669656828834e-7
    B5 = +1.09882180475307e-10
    
    C0 = +9.652117566301735e-1
    C1 = +5.563942679876470e-2  
    C2 = +2.934835672207024e-2
    C3 = +7.696735279663661e-3  
    C4 = -6.147609081030884e-3
    C5 = -1.931115655395969e-3  
    C6 = +6.350841470341581e-4
    C7 = +1.682282887657391e-4

    
    if T_x > 3.254e2:
        Hydrate_Equilibrium_Pressure = 1.0e10
        
    elif (T_x > 2.75e2 and T_x <= 3.254e2):
        Hydrate_Equilibrium_Pressure = 1.0e6*np.exp( A0+T_x*
                                          (A1+T_x*
                                          (A2+T_x*
                                          (A3+T_x*
                                          (A4+T_x*
                                           A5)))) 
                                         )

                                                  
    elif (T_x > 2.745e2 and T_x <= 2.75e2):
        F1 = 1.0e6*np.exp( A0+T_x*(A1+T_x*(A2+T_x*(A3+T_x*(A4+T_x*A5)))))
        
        T_c = T_x - 273.15e0

        F2 = 1.0e6*np.exp(C0+T_c * (C1+T_c * (C2+T_c * (C3+T_c * (C4+T_c * (C5+T_c * (C6+T_c*C7)))))))
        
        Hydrate_Equilibrium_Pressure = F1 + 2.0e0*(2.75e2 - T_x)*(F2 - F1)

    elif (T_x > 2.72e2 and T_x <= 2.745e2):
        T_c = T_x - 273.15e0
        Hydrate_Equilibrium_Pressure = 1.0e6*np.exp( C0+T_c*(C1+T_c*(C2+T_c*(C3+T_c*(C4+T_c*(C5+T_c*(C6+T_c*C7)))))))

    elif (T_x > 2.715e2 and T_x <= 2.72e2):
        
        F1 = 1.0e6*np.exp( B0+T_x*  
                        (B1+T_x*  
                        (B2+T_x*  
                        (B3+T_x*  
                        (B4+T_x*  
                         B5))))   
                       )

        T_c = T_x-273.15e0

        F2 = 1.0e6*np.exp( C0+T_c*  
                        (C1+T_c*  
                        (C2+T_c*  
                        (C3+T_c*  
                        (C4+T_c*  
                        (C5+T_c*  
                        (C6+T_c*  
                         C7)))))) 
                       )

        Hydrate_Equilibrium_Pressure = F2 + 2.0e0*(2.72e2 - T_x)*(F1 - F2)

    elif (T_x <= 2.715e2 and T_x >= 1.488e2):
        Hydrate_Equilibrium_Pressure = 1.0e6*np.exp( B0+T_x* 
                                          (B1+T_x* 
                                          (B2+T_x* 
                                          (B3+T_x* 
                                          (B4+T_x* 
                                           B5))))  
                                         )
                                             

    else: 
        print('T is out of range')
    


    return Hydrate_Equilibrium_Pressure





v_Equil_P=np.vectorize(Equil_P)


def Tshift_NaCl(X_iA,Max_Tshift=2e0,Xmol_iA_atMax_Tshift=1.335e-2,InhibitorMW=5.8448e1):
    #For NaCl
    
    # Xmol_iA_atMax_Tshift = 1.335E-2
    # Max_Tshift = 2.0E0
    # InhibitorMW = 5.8448E1
    # 
    MW_H2O = 18.015268e+00
    
    moles_Nacl=(X_iA)/InhibitorMW
    moles_H2O=(1-X_iA)/MW_H2O
    
    Xmol_iA=moles_Nacl/(moles_Nacl+moles_H2O)
    
    Tshift=Max_Tshift*np.log(1-Xmol_iA)/np.log(1-Xmol_iA_atMax_Tshift)
    
    return Tshift


v_Tshift_NaCl=np.vectorize(Tshift_NaCl)


#%% Executable
def RunTOUGH(file,th_path=r'path_to_TOUGH_executable'):
    
    if th_path == r'path_to_TOUGH_executable':
        print('Run function as RunTOUGH(file, th_path = <correct path to executable>)')

    from subprocess import Popen,PIPE
    import os
    
    wd=os.getcwd()
    fd=os.path.dirname(file)
    fname=os.path.basename(file)
    
    nd=os.path.join(wd,fd)
    
    old_wd = os.getcwd()

    os.chdir(nd)
    
    op_files=os.listdir()
    up_file=r'Parameter_Update_File'
    
    
    
    f_list = [fname,up_file,'.ipynb_checkpoints']    
    for f in op_files:
        if not f in f_list:
            os.remove(f)
    
    
    
    file_in=fname
    file_out=fname[:-2]+r'out'
    cmd=th_path+r' <'+file_in+r'> '+file_out 
    
    print('Start simulation at '+os.path.basename(os.getcwd()))
    # proc=Popen(cmd, stdout=PIPE, shell=True)
    proc=Popen(cmd, stdout=PIPE, shell=True, preexec_fn=os.setpgrp)

    os.chdir(old_wd)

    return proc