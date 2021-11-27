
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


def minmax(pmt, data):
    from bokeh.models import LinearColorMapper
    from bokeh.palettes import Viridis256



    minv=data[pmt].min()
    maxv=data[pmt].max()
    fmt='0.0'
    
    
    #Parameters that are unit fractions or multipliers [0,1]
    fr_pmts=['S_hyd', 'S_aqu', 'S_gas', 'S_icd', 'k_rg', 'k_rw', 'k_adj_F']
    
    param_name={'P':'pressure [MPa]',
                'T':'temperature [degC]',
                'S_hyd':'hydrate sat. [V/V]',
                'S_aqu':'water sat. [V/V]',
                'S_gas':'gas sat. [V/V]',
                'S_icd':'ice sat. [V/V]',
                'X_inh':'salinity',
                'k_rg':'gas rel. perm. [frac]',
                'k_rw':'water rel. perm. [frac]',
                'krg_h':'gas rel. perm. [frac]',
                'krw_h':'water rel. perm. [frac]',
                'k_adj_F':'perm adj. factor [frac]',
                'perm_abs':'absolute perm. [mD]',
                'k_eff':'eff. perm. [mD]',
                'porosity':'porosity [frac]',
                'phi_eff':'eff. porosity [frac]',
                'P_cap':'capillay pres. [Pa]',
                'P_aqu':'press. aqueous ph. [Pa]',
                'ZONE':'zonation index'}
    
    try:
        pmt_name=param_name[pmt]
    except:
        pmt_name='unnamed'
            
    if pmt in fr_pmts or minv-maxv==0:
        minv=0
        maxv=1
        fmt='0.0'
    
    elif minv==0:
        minv=0
        
        if np.log10(maxv-minv)<-16:
            minv=0
            maxv=1

        
    elif maxv==0:
        
        if minv==0:
            maxv=0.5
        
        else:
            maxv=0
    
    
    else:
        if maxv<0:
            q_max=np.floor(np.log10(-maxv))-1
        
        elif maxv==0:
            q_max=0
        
        else:
            q_max=np.floor(np.log10(maxv))-1
            
        

        if minv<0:
            q_min=np.floor(np.log10(-minv))-1
        elif minv==0:
            q_min=0
        else:
            q_min=np.floor(np.log10(minv))-1
        
        
        maxv=10**q_max * np.ceil(maxv/(10**q_max))

        minv=10**q_min * np.floor(minv/10**q_min)

        
            
    if pmt=='T' and minv>0:
        minv=0
        
    # if pmt in ['perm_abs','k_eff']:

    #     print('logarithmic colorbar')
    #     cm=LogColorMapper(palette=Viridis256, low=minv, high=maxv)

        
    # else:
    cm=LinearColorMapper(palette=Viridis256, low=minv, high=maxv)


    return minv,maxv,pmt_name,cm


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