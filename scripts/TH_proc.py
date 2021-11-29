#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 15:37:24 2020

@author: gpb
"""

#%% Init functions
def init_HYDRATE():
    Hyd_dict={'H_1':{'NCom':1},
              'H_2':{'nameG':"'CH4'",
                     'hydrN':6,
                     'moleF':1},
              'H_3':{'N_ThC':1},
              'H_4':{'ThC':[4.5e-1]},
              'H_5':{'N_SpH':1},
              'H_6':{'SpH':[2.1e3]},
              'H_7':{'N_Rho':1},
              'H_8':{'Rho':[9.2e2]},
              'H_9':{'inhibitor_flag':False,
                     'Max_TShift':2,
                     'Y_atMax_TShift':1.335e-2,
                     'InhibitorMW':58.448,
                     'InhibitorDens':2.6e3,
                     'InhibitorEnthSol':6.6479e4,
                     'InhibitorCpCoeff':[41.293E0,3.3607e-2,-1.3927e-5]},
              'H_10':{'EquationOption':0},
              'H_11':{'Reaction_Type':"'EQUILIBRIUM'"},
              'H_12':{'ActivationEnergy':8.1e4,
                      'IntrinsicRateConstant':4.70848E+05,
                      'Area_Factor':1}
              }
    return Hyd_dict


def init_MEMORY():
    mem_dict={'M_2':{'EOS_Name':"'HYDRATE-EQUILIBRIUM'"},
            'M_3':{'NumCom':2,
                    'NumEq':3,
                    'NumPhases':4,
                    'binary_diffusion':False},
            'M_4':{'coordinate_system':"'Cartesian'",
                    'Max_NumElem':15,
                    'Max_NumConx':30,
                    'ElemNameLength':5,
                    'active_conx_only':False,
                    'boundaries_in_matrix':False},
            'M_5':{'Max_NumSS':2},
            'M_6':{'Max_NumMedia':2},          
            'M_7':{'element_by_element_properties':False,
                    'porosity_perm_dependence':False,
                    'scaled_capillary_pressure':False,
                    'Option_tortuosity_CompuMethod':"'Saturation'"},
            'M_8':{'coupled_geochemistry':False,
                    'property_update':"'Continuous'"},
            'M_9':{'coupled_geomechanics':False,
                    'geomechanical_code_name':"' '",
                    'property_update':"'Continuous'" ,
                    'num_geomech_param':0}

            }

    return mem_dict


def init_INCON():
    import numpy as np
    import pandas as pd

    dtypes = np.dtype([
          ('ElName', 'S5'),
          ('NSEQ', int),
          ('NADD', int),
          ('porosity', float),
          ('StateIndex', 'S3'),
          ('comment','S36'),
          ('PermX',float),
          ('PermY',float),
          ('PermZ',float),          
          ('phi',float),          
          ('var1',float),
          ('var2',float),
          ('var3',float),
          ('var4',float),
          ('var5',float),
          ('var6',float)
          ])
    data = np.empty(0, dtype=dtypes)
    empty_INCON = pd.DataFrame(data)

    return empty_INCON

def init_SUBDOMAIN(*names):
    import pandas as pd
    
    n_sd=1
    
    subdom=dict()
    
    reg_columns=['definition_mode',
                'number_of_elements',
                'format_to_read_data',
                'first_element_number',
                'first_element_name',
                'element_sequence_stride',
                'region_shape',
                'elements',
                'Xmin',
                'Xmax',
                'Ymin',
                'Ymax',
                'Zmin',
                'Zmax',
                'Rmin',
                'Rmax',
                'top_cylinder_center_xyz',
                'bot_cylinder_center_xyz',
                'sphere_center_xyz']
    
    # if isinstance(names, list):
    #     n_sd=len(names)
    #     subdom['number_of_subdomains'] = n_sd
    n_sd=0
    subdom['number_of_subdomains'] = n_sd
    for name in names:
        subdom[name]={'subdomain_name':name,
                      'number_of_regions':0,
                      'regions':pd.DataFrame(columns=reg_columns)}
        n_sd+=1
        
    subdom['number_of_subdomains'] = n_sd
    
    if n_sd == 0:
       print("No subdomain defined. Default subdomain structure is created")

       subdom={'number_of_subdomains':n_sd,
               'subdomain_1':{'subdomain_name':'',
                              'number_of_regions':0,
                              'regions':pd.DataFrame(columns=reg_columns)}}
        
    return subdom


def init_RPCAP():
    import pandas as pd

    rpcap = pd.DataFrame(columns=['RelPermEquationNum','RelPermParam1','RelPermParam2',
                                  'RelPermParam3','RelPermParam4','RelPermParam5',
                                  'RelPermParam6','RelPermParam7',
                                  'PcapEquationNum','PcapParam1','PcapParam2','PcapParam3',
                                  'PcapParam4','PcapParam5','PcapParam6','PcapParam7'],
                        index = range(2) )
                    
    return rpcap

#%%Reading Functions
def read_HYDRATE(hydrate_str):
    import re
    hydrate=init_HYDRATE()

    if type(hydrate_str)==list:
        hydrate_str='\n'.join(hydrate_str)

    hydrate_str=hydrate_str.strip()
    hydrate_str=hydrate_str.split('\n')

    lines = iter(range(1,len(hydrate_str)+1))

    #Process HYDRATE_1
    line = next(lines)
    proc_line=hydrate_str[line].split('!')[0].strip()
    hydrate['H_1']['NCom']=int(proc_line)

    #Process HYDRATE_2
    line = next(lines)
    proc_line=hydrate_str[line].split('!')[0].strip().split()

    #remove FORTRAN formats
    proc_line=proc_line[0]+'  '+re.sub('d','e','  '.join(proc_line[1:]))
    proc_line=proc_line.split()


    hydrate['H_2']['nameG']=proc_line[0]
    hydrate['H_2']['hydrN']=float(proc_line[1])
    hydrate['H_2']['moleF']=float(proc_line[2])

    #Process HYDRATE_3
    line = next(lines)
    proc_line=hydrate_str[line].split('!')[0].strip()
    hydrate['H_3']['N_ThC']=int(proc_line)

    #Process HYDRATE_4
    line = next(lines)
    proc_line=hydrate_str[line].split('!')[0].strip().split()
    hydrate['H_4']['ThC']=list(map(float,proc_line))

    #Process HYDRATE_5
    line = next(lines)
    proc_line=hydrate_str[line].split('!')[0].strip()
    hydrate['H_5']['N_SpH']=int(proc_line)

    #Process HYDRATE_6
    line = next(lines)
    proc_line=hydrate_str[line].split('!')[0].strip().split()
    hydrate['H_6']['SpH']=list(map(float,proc_line))    

    #Process HYDRATE_7
    line = next(lines)
    proc_line=hydrate_str[line].split('!')[0].strip()
    hydrate['H_7']['N_Rho']=int(proc_line)

    #Process HYDRATE_8
    if hydrate['H_7']['N_Rho'] > 0:
        line = next(lines)
        proc_line=hydrate_str[8].split('!')[0].strip().split()
        hydrate['H_8']['Rho']=list(map(float,proc_line))
    else:
        hydrate['H_8']['Rho']=None
    #Process HYDRATE_9
    line = next(lines)
    proc_line=hydrate_str[line].split('!')[0].strip().split()

    #remove FORTRAN formats
    proc_line=proc_line[0].strip('.').capitalize() + '  '+re.sub('d','e','  '.join(proc_line[1:]))
    proc_line=proc_line.split()
    hydrate['H_9']['inhibitor_flag']=eval(proc_line[0])   

    if hydrate['H_9']['inhibitor_flag']:
        hydrate['H_9']['Max_TShift']=float(proc_line[1])
        hydrate['H_9']['Y_atMax_TShift']=float(proc_line[2])
        hydrate['H_9']['InhibitorMW']=float(proc_line[3])
        hydrate['H_9']['InhibitorDens']=float(proc_line[4])
        hydrate['H_9']['InhibitorEnthSol']=float(proc_line[5])
        hydrate['H_9']['InhibitorCpCoeff']=list(map(float,proc_line[6:])) 

    else:
        hydrate['H_9']['Max_TShift']=None
        hydrate['H_9']['Y_atMax_TShift']=None
        hydrate['H_9']['InhibitorMW']=None
        hydrate['H_9']['InhibitorDens']=None
        hydrate['H_9']['InhibitorEnthSol']=None
        hydrate['H_9']['InhibitorCpCoeff']=None

    #Process HYDRATE_10
    line = next(lines)
    proc_line=hydrate_str[line].split('!')[0].strip()
    hydrate['H_10']['EquationOption']=int(proc_line)

    #Process HYDRATE_11
    line = next(lines)
    proc_line=hydrate_str[line].split('!')[0].strip()
    hydrate['H_11']['Reaction_Type']=proc_line

    #Process HYDRATE_12
    if 'KINETIC' in hydrate['H_11']['Reaction_Type']:
        line = next(lines)
        proc_line=hydrate_str[line].split('!')[0].strip()
        proc_line=re.sub('d','e',proc_line)
        proc_line=proc_line.split()

        hydrate['H_12']['ActivationEnergy']=float(proc_line[0])
        hydrate['H_12']['IntrinsicRateConstant']=float(proc_line[1])
        hydrate['H_12']['Area_Factor']=float(proc_line[2])

    else:
        hydrate['H_12']['ActivationEnergy']=None
        hydrate['H_12']['IntrinsicRateConstant']=None
        hydrate['H_12']['Area_Factor']=None

    return hydrate

def read_MEMORY(memory_str):
    import re

    memory=init_MEMORY()

    if type(memory_str) == list:
        memory_str = '\n'.join(memory_str)

    memory_str = memory_str.strip()
    memory_str = memory_str.split('\n')

    #Process MEMORY_2
    proc_line = memory_str[1].split('!')[0].strip()
    memory['M_2']['EOS_Name']=proc_line

    #Process MEMORY_3
    proc_line = memory_str[2].split('!')[0].strip().split()
    memory['M_3']['NumCom']=int(proc_line[0])
    memory['M_3']['NumEq']=int(proc_line[1])
    memory['M_3']['NumPhases']=int(proc_line[2])
    memory['M_3']['binary_diffusion']=eval(proc_line[3].strip('.').capitalize())

    #Process MEMORY_4
    proc_line = memory_str[3].split('!')[0].strip().split()
    memory['M_4']['coordinate_system']=proc_line[0]
    memory['M_4']['Max_NumElem']=int(proc_line[1])
    memory['M_4']['Max_NumConx']=int(proc_line[2])
    memory['M_4']['ElemNameLength']=int(proc_line[3])
    memory['M_4']['active_conx_only']=eval(proc_line[4].strip('.').capitalize())
    memory['M_4']['boundaries_in_matrix']=eval(proc_line[5].strip('.').capitalize())

    #Process MEMORY_5
    proc_line = memory_str[4].split('!')[0].strip().split()
    memory['M_5']['Max_NumSS']=int(proc_line[0])

    #Process MEMORY_6
    proc_line = memory_str[5].split('!')[0].strip().split()
    memory['M_6']['Max_NumMedia']=int(proc_line[0])

    #Process MEMORY_7
    proc_line = memory_str[6].split('!')[0].strip().split()
    memory['M_7']['element_by_element_properties'] = eval(proc_line[0].strip('.').capitalize())
    memory['M_7']['porosity_perm_dependence'] = eval(proc_line[1].strip('.').capitalize())
    memory['M_7']['scaled_capillary_pressure'] = eval(proc_line[2].strip('.').capitalize())
    memory['M_7']['Option_tortuosity_CompuMethod']=proc_line[3]

    #Process MEMORY_8
    proc_line = memory_str[7].split('!')[0].strip().split()
    memory['M_8']['coupled_geochemistry'] = eval(proc_line[0].strip('.').capitalize())
    memory['M_8']['property_update']=proc_line[1]

    #Process MEMORY_9
    proc_line = memory_str[8].replace("' '", "'_'")
    proc_line = proc_line.split('!')[0].strip().split()
    memory['M_9']['coupled_geomechanics'] = eval(proc_line[0].strip('.').capitalize())
    memory['M_9']['geomechanical_code_name']=proc_line[1].replace("'_'","' '")
    memory['M_9']['property_update']=proc_line[2]
    memory['M_9']['num_geomech_param']=int(proc_line[3])

    return memory


def read_ROCKS(rocks_str):
    
    #Read ROCKS section to DataFrame
    
    import pandas as pd
    from io import StringIO
    import os
    
    # if os.path.isfile(rocks_str):
    #     with open(rocks_str) as r_file:
    #         rocks_str=r_file.read()
    
    if type(rocks_str)==list:
        rocks_str='\n'.join(rocks_str)
        
    
    rocks_str=rocks_str.strip()
    rocks_str=rocks_str.split('\n')
    rocks_str=rocks_str[1:]
    
    rtable=pd.DataFrame(columns=['Name','NAD','DensG','Poros','Perm1','Perm2','Perm3','KThrW','SpcHt','PoMedRGrain',
                             'Compr','Expan','KThrD','Tortu','Klink','OrgCF','CritSat','PermExpon','Beta','Gama','PhiZeroStress',
                             'RelPermEquationNum','RelPermParam1','RelPermParam2','RelPermParam3','RelPermParam4','RelPermParam5','RelPermParam6','RelPermParam7',
                             'PcapEquationNum','PcapParam1','PcapParam2','PcapParam3','PcapParam4','PcapParam5','PcapParam6','PcapParam7',
                             'PhiPolyOrder','PhiCoeff0','PhiCoeff1','PhiCoeff2','PhiCoeff3','PhiCoeff4','PhiCoeff5','PhiCoeff6',
                             'LoComp','SatAtLoComp','HiComp','SatAtHiComp','DeltaSat'])
    
    
    def ROCKS1(line,df=rtable):
        r1_header=df.columns[:10]
        r1_widths=[5]*2+[10]*8

        r1_dtype=dict()
        for i,header in enumerate(r1_header):
            if i==0:
                r1_dtype[header]=str
            elif i==1:
                r1_dtype[header]=int
            else:
                r1_dtype[header]=float

        r1=pd.read_fwf(StringIO(line),header=None,names=r1_header,widths=r1_widths, dtype=r1_dtype)
        return r1.loc[0]
    
    def ROCKS2(line,df=rtable):
        r2_header=df.columns[10:21]
        r2_widths=[10]*11

        r2_dtype=float
        try:
            r2=pd.read_fwf(StringIO(line),header=None,names=r2_header,widths=r2_widths, dtype=r2_dtype)
        except:
            r2=pd.read_fwf(StringIO(line),header=None,names=r2_header,widths=r2_widths)
        return r2.loc[0]
        
    def ROCKS3(line,df=rtable):
        r3_header=df.columns[21:29]
        r3_widths=[10]*8

        r3_dtype=dict(zip(r3_header,[int]+7*[float]))

        r3=pd.read_fwf(StringIO(line),header=None,names=r3_header,widths=r3_widths, dtype=r3_dtype)
        return r3.loc[0]

    def ROCKS4(line,df=rtable):
        r4_header=df.columns[29:37]
        r4_widths=[10]*8

        r4_dtype=dict(zip(r4_header,[int]+7*[float]))

        r4=pd.read_fwf(StringIO(line),header=None,names=r4_header,widths=r4_widths, dtype=r4_dtype)
        return r4.loc[0]
    
    r_iter=iter(rocks_str)
    while True:     
        try:
            line=next(r_iter)
            line1=ROCKS1(line)
            # print(line)
            nad=line1['NAD']
            rock=line1
            
            if nad>=1:
                line=next(r_iter)

                line2=ROCKS2(line)
                rock=rock.append(line2)
                    

                
                if nad>=2:
                    line=next(r_iter)

                    line3=ROCKS3(line)
                    
                    line=next(r_iter)
                    line4=ROCKS4(line)
                    rock=rock.append(line3)
                    rock=rock.append(line4)
            
            rtable=rtable.append(rock,ignore_index=True)
            

        except:
            break
    
    rtable['RelPermEquationNum']= rtable.RelPermEquationNum.fillna(0).astype(int)
    rtable['PcapEquationNum']= rtable.PcapEquationNum.fillna(0).astype(int)
    return rtable

def read_RPCAP(rpcap_str):
    import pandas as pd
    from io import StringIO

    if type(rpcap_str) == list:
        rpcap_str = '\n'.join(rpcap_str)

    rpcap_str = rpcap_str.strip()
    rpcap_str = rpcap_str.split('\n')
    rpcap_str = rpcap_str[1:]
    rpcap_str = '\n'.join(rpcap_str)
    
    rpcap_header=['EquationNum','Param1','Param2',
                    'Param3','Param4','Param5',
                    'Param6','Param7']

    rpcap_widths=[10]*8

    rpcap_dtype=dict(zip(rpcap_header,[int]+7*[float]))

    rpcap=pd.read_fwf(StringIO(rpcap_str),header=None,names=rpcap_header,
                    widths=rpcap_widths, dtype=rpcap_dtype)

    rpcap = rpcap.set_index(pd.Index(['RelPerm', 'CapPres']))

    return rpcap

            
        

def read_INCON(INCON_data):
    import pandas as pd
    import numpy as np
    
    
    #Clean up header and tail
    if type(INCON_data)==list:
        INCON_data='\n'.join(INCON_data)
        

    INCON_data = INCON_data.strip()
    INCON_data = INCON_data.split(':::')
    INCON_data = INCON_data[0].strip()
    INCON_data = INCON_data.split('\n')
    INCON_data = INCON_data[1:]
    

    if len(INCON_data) == 0:
        print("No INCON data processed")
        return []
    
    
    
        
    """
    INCON_1
    Format (A5, 2I5, E15.8, 2x, A3, 36x, 3(E15.8))
    ElName5C, NSEQ, NADD, porosity, StateIndex, (perm(i), i=1,3)
    """
    INCON_1=pd.DataFrame(data=INCON_data[::2],columns=['raw1'])
    INCON_1['ElName']=INCON_1['raw1'].str[:5]
    INCON_1['NSEQ']=INCON_1['raw1'].str[5:10]
    INCON_1['NADD']=INCON_1['raw1'].str[10:15]
    INCON_1['porosity']=INCON_1['raw1'].str[15:30]
    INCON_1['StateIndex']=INCON_1['raw1'].str[32:35]
    INCON_1['comment']=INCON_1['raw1'].str[35:71]
    INCON_1['PermX']=INCON_1['raw1'].str[71:86]
    INCON_1['PermY']=INCON_1['raw1'].str[86:101]
    INCON_1['PermZ']=INCON_1['raw1'].str[101:116]
    INCON_1['phi']=INCON_1['raw1'].str[116:131]
    INCON_1=INCON_1.set_index('ElName',drop=False)
    
    cols_1=list(INCON_1.columns)
    

    """
    INCON_2
    Format (6E20.13)
    X(i), i = 1,NumCom+1
    """
    # INCON_2_raw=pd.Series(data=INCON_data[1::2])
    INCON_2=pd.DataFrame(data=INCON_data[1::2],columns=['raw2'])
    col_start=0
    col_length=20
    for col_idx in range(6):
        new_col='var{:d}'.format(col_idx+1)
        
        # print('processing', new_col)
        col_end=col_start+col_length
        try:
            INCON_2[new_col]=INCON_2['raw2'].str[col_start:col_end].astype(float)
        except:
            INCON_2[new_col]=np.nan
        col_start=col_end


    INCON_2['ElName']=INCON_1.index

    INCON_2=INCON_2.set_index('ElName')

    cols_2=list(INCON_2.columns)

    cols=cols_1[1:]+cols_2[:-1]+cols_1[:1]+cols_2[-1:]
    
    INCON=pd.concat([INCON_1,INCON_2],axis=1,sort=False)
    
    INCON=INCON[cols]

    INCON.loc[:,'NSEQ']=pd.to_numeric(INCON.loc[:,'NSEQ'],'coerce')
    INCON.loc[:,'NADD']=pd.to_numeric(INCON.loc[:,'NADD'],'coerce')
    INCON.loc[:,'porosity']=pd.to_numeric(INCON.loc[:,'porosity'],'coerce')
    INCON.loc[:,'PermX']=pd.to_numeric(INCON.loc[:,'PermX'],'coerce')
    INCON.loc[:,'PermY']=pd.to_numeric(INCON.loc[:,'PermY'],'coerce')
    INCON.loc[:,'PermZ']=pd.to_numeric(INCON.loc[:,'PermZ'],'coerce')
    INCON.loc[:,'phi']=pd.to_numeric(INCON.loc[:,'phi'],'coerce')


    return INCON


def read_CONNE(CONNE_data):
     
    import pandas as pd
    
    if type(CONNE_data)==list:
        CONNE_data='\n'.join(CONNE_data)
        

    CONNE_data=CONNE_data.strip()
    CONNE_data=CONNE_data.split('\n')
    conne=CONNE_data[1:]
    
    
    conne=pd.DataFrame(data=conne,columns=['raw'])
    conne.loc[:,'raw']=conne['raw'].str.lstrip()
    ftr=conne['raw']!=''
    conne=conne[ftr]
    
    """
    Format (A5, A5, 4I5, 5E10.4)
    ConxName1, ConxName2, NSEQ, NAD1, NAD2, ConxKi,
    ConxD1, ConxD2, ConxArea, ConxBeta, emissivity
    """
    
    conne['ConxName1']=conne['raw'].str[:5]     #A5
    conne['ConxName2']=conne['raw'].str[5:10]   #A5
    conne['NSEQ']=conne['raw'].str[10:15]       #I5
    conne['NAD1']=conne['raw'].str[15:20]       #I5
    conne['NAD2']=conne['raw'].str[20:25]       #I5
    conne['ConxKi']=conne['raw'].str[25:30]     #I5
    conne['ConxD1']=conne['raw'].str[30:40]     #E10.4
    conne['ConxD2']=conne['raw'].str[40:50]     #E10.4
    conne['ConxArea']=conne['raw'].str[50:60]   #E10.4
    conne['ConxBeta']=conne['raw'].str[60:70]   #E10.4
    conne['emissivity']=conne['raw'].str[70:80] #E10.4
    
    
    cols=list(conne.columns)
    cols=cols[1:]+cols[:1]
    conne=conne[cols]
    
    return conne


def read_ELEME(ELEME_data):

    import pandas as pd

    #Clean up string data header and tail
    if type(ELEME_data)==list:
        ELEME_data='\n'.join(ELEME_data)
        

    ELEME_data=ELEME_data.strip()
    ELEME_data=ELEME_data.split('\n')
    ELEME_data=ELEME_data[1:]
   
    #Create empty dataframel
    elem_data = pd.DataFrame()
    elem_data['raw'] = ELEME_data
    
    #Retrieve indexes for lines with fix boundaries
    ina_query=elem_data['raw'].str.contains('ina')
    try:
        INA_idx=elem_data[ina_query].index[0]
        
        #fix boundary elements
        fixb=elem_data.loc[INA_idx+1:,'raw']
        
        length_query=fixb.str.len()>=81
        
        #truncate long lines
        fixb[length_query]=fixb[length_query].str[:81]+'I  '
        
        #complete short lines
        fixb[~length_query]=fixb[~length_query].str.pad(width=81,side='right')
        fixb[~length_query]=fixb[~length_query].astype(str)+'I  '
        
        #Delete ina line
        elem_data=elem_data.drop(INA_idx)
        #Reset index
        elem_data=elem_data.reset_index(drop=True)
    except:
        pass
    

    
    #Parse elem_data    
    elem_data['ElName']=    elem_data['raw'].str[:5]
    elem_data['NSEQ']=      elem_data['raw'].str[5:10]
    elem_data['NADD']=      elem_data['raw'].str[10:15]
    elem_data['MA12']=      elem_data['raw'].str[15:20]
    elem_data['elem_vol']=  pd.to_numeric(elem_data['raw'].str[20:30],errors='coerce')
    elem_data['elem_aht']=  pd.to_numeric(elem_data['raw'].str[30:40],errors='coerce')
    elem_data['elem_pm']=   pd.to_numeric(elem_data['raw'].str[40:50],errors='coerce')
    elem_data['X']=         pd.to_numeric(elem_data['raw'].str[50:60],errors='coerce')
    elem_data['Y']=         pd.to_numeric(elem_data['raw'].str[60:70],errors='coerce')
    elem_data['Z']=         pd.to_numeric(elem_data['raw'].str[70:80],errors='coerce')
    elem_data['elem_activity']=    elem_data['raw'].str[80:]
    
    return elem_data



def read_PARAM(param_str):       
    import pandas as pd
    from io import StringIO
    
    if type(param_str)==list:
        param_str='\n'.join(param_str)
        
    
    param_str=param_str.strip()
    param_str=param_str.split('\n')
    param_str=param_str[1:]
    
    param_1_header=['Max_NumNRIterations', 'OutputOption', 'Max_NumTimeSteps', 'iCPU_MaxTime', 'PRINT_frequency',
                     'MOP','BaseDiffusionCoef', 'DiffusionExpon', 'DiffusionStrength','SAVE_frequency', 
                     'TimeSeries_frequency' ]
    param_1_widths=[2]*2+3*[4]+[24]+3*[10]+2*[5]

    # param_1_dtype=dict(zip(param_1_header,6*[np.int32]+3*[np.float64]+2*[np.int32]))

    param_1=pd.read_fwf(StringIO(param_str[0]),header=None,names=param_1_header,widths=param_1_widths,dtype=str).loc[0]
    param_1['MOP_items']=dict(zip(range(1,len(param_1.MOP)+1),list(str(param_1.MOP))))


    # param_1.MOP=

    
    
    param_2_header=['TimeOrigin', 'SimulationTimeEnd', 'InitialTimeStep',
                    'MaxTimeStep', 'TrackElemName:', 'empty_space', 'Gravity',
                    'Dt_reducer','Scale']
    
    param_2_widths=4*[10]+2*[5]+3*[10]
    param_2=pd.read_fwf(StringIO(param_str[1]),header=None,names=param_2_header,widths=param_2_widths,dtype=str).loc[0]

    param_3_header=['rel_convergence_crit', 'abs_convergence_crit', 'U_p', 'W_upstream',
                   'W_NRIteration', 'derivative_increment', 'W_implicitness', 'empty1',
                   'DefaultStateIndex', 'empty2', 'P_overshoot', 'T_overshoot',
                   'S_overshoot']
    param_3_widths=7*[10]+[2]+[3]+[5]+3*[10]
    param_3=pd.read_fwf(StringIO(param_str[2]),header=None,names=param_3_header,widths=param_3_widths,dtype=str).loc[0]
    
    
    param_4_header=['var1', 'var2', 'var3', 'var4', 'var5', 'var6']
    param_4_widths=6*[20]
    param_4=pd.read_fwf(StringIO(param_str[3]),header=None,names=param_4_header,widths=param_4_widths,dtype=str).loc[0]
    
    param = pd.concat([param_1,param_2,param_3,param_4],keys=['P_1','P_2','P_3','P_4'])
    
    return param



def read_INDOM(INDOM_data):
    import pandas as pd
    
    #Clean up header and tail
    if type(INDOM_data)==list:
        INDOM_data='\n'.join(INDOM_data)
        

    INDOM_data=INDOM_data.strip()
    INDOM_data=INDOM_data.split('\n')
    INDOM_data=INDOM_data[1:]
    
    """
    INDOM_1
    Format (A5, 2x, A3)
    ElName5C, NSEQ, NADD, porosity, StateIndex, (perm(i), i=1,3)
    """
    INDOM_1=pd.DataFrame(data=INDOM_data[::2],columns=['raw1'])
    INDOM_1['Rk_name']=INDOM_1['raw1'].str[:5]
    INDOM_1['StateIndex']=INDOM_1['raw1'].str[7:10]
    INDOM_1=INDOM_1.set_index('Rk_name',drop=False)
    
    cols_1=list(INDOM_1.columns)
    

    """
    INDOM_2
    Format (6E20.13)
    X(i), i = 1,NumCom+1
    """
    INDOM_2=pd.DataFrame(data=INDOM_data[1::2],columns=['raw2'])
    INDOM_2['var1']=INDOM_2['raw2'].str[:20]
    INDOM_2['var2']=INDOM_2['raw2'].str[20:40]
    INDOM_2['var3']=INDOM_2['raw2'].str[40:60]
    INDOM_2['var4']=INDOM_2['raw2'].str[60:80]
    INDOM_2['ElName']=INDOM_1.index
    INDOM_2=INDOM_2.set_index('ElName')

    cols_2=list(INDOM_2.columns)

    cols=cols_1[1:]+cols_2[1:]+cols_1[:1]+cols_2[:1]
    
    INDOM=pd.concat([INDOM_1,INDOM_2],axis=1,sort=False)
    
    INDOM=INDOM[cols]
    return INDOM


def read_MESH_file(file):
    import io
    
    #Read file
    if isinstance(file,io.StringIO):
        text=file.getvalue()
    
    else:
        with open(file,'r',encoding='cp1250') as f:
            text = f.read()
        
    conne_idx=text.find('CONNE')
    
    eleme_data=text[:conne_idx]
    conne_data=text[conne_idx:]
    
    eleme=read_ELEME(eleme_data)
    conne=read_CONNE(conne_data)
    
    return eleme,conne


#%% Write Functions

def check_MEMORY(data):
    import numpy as np
    import pandas as pd
    #Read HYDRATE section
    if 'EQUILIBRIUM' in data.HYDRATE.processed['H_11']['Reaction_Type']:
        data.MEMORY.processed['M_2']['EOS_Name'] = "'HYDRATE-EQUILIBRIUM'"
        if data.HYDRATE.processed['H_9']['inhibitor_flag']:
            NumCom,NumEqu,NumPhases = 3,4,4
        else:
            NumCom,NumEqu,NumPhases = 2,3,4
    else:
        data.MEMORY.processed['M_2']['EOS_Name'] = "'HYDRATE-KINETIC'"
        if data.HYDRATE.processed['H_9']['inhibitor_flag']:
            NumCom,NumEqu,NumPhases = 3,4,4
        else:
            NumCom,NumEqu,NumPhases = 4,5,4

    
    data.MEMORY.processed['M_3']['NumCom'] = NumCom
    data.MEMORY.processed['M_3']['NumEq'] = NumEqu
    data.MEMORY.processed['M_3']['NumPhases'] = NumPhases

    # Read ROCKS section
    data.MEMORY.processed['M_6']['Max_NumMedia'] = data.ROCKS.processed.shape[0]
    
    #Read grid specs from ELEME and CONNE sections
    Max_NumElem = 5*(1+data.ELEME.processed.shape[0]//5)
    Max_NumConx = 5*(1+data.CONNE.processed.shape[0]//5)
    
    data.MEMORY.processed['M_4']['Max_NumElem'] = Max_NumElem
    data.MEMORY.processed['M_4']['Max_NumConx'] = Max_NumConx

    if type(data.GENER.processed) is pd.DataFrame and len(data.GENER.raw)-3 > 0:
        MaxNum_SS = data.GENER.processed.shape[0]
        data.MEMORY.processed['M_5']['Max_NumSS'] = MaxNum_SS

    data.MEMORY.modified = True

    return data


def write_MEMORY(memory):

    fmt_lst=['{:s}',
            '{:5d}{:5d}{:5d}  .{:s}.',
            '{:s} {:6d}{:6d}{:6d} .{:s}. .{:s}.',
            '{:4d}',
            '{:4d}',
            '.{:s}. .{:s}. .{:s}. {:s}',
            '.{:s}. {:s}',
            '.{:s}. {:s} {:s}{:5d}']


    fmts=dict(zip(memory.keys(),fmt_lst))

    str_lst=['MEMORY']


    #iterate over dictionary
    for key,vals in memory.items():
        val_lst=[]

        subvals=list(memory[key].values())
    #     print('\nkey:',key)
        #Iterate over each sub dictionary and extract values
        for val in subvals:
    #         print(val)

            #Process iterables
            if type(val)==list:
                val_lst+=val


            else:
                #Check for flags
                if type(val)==bool:

                    flag = val
                    val=str(val).upper()


                #Update list of iterables
                val_lst.append(val)
            
    #     print(val_lst)

        #Format iterable
    #     print(fmts[key])
        line=fmts[key].format(*val_lst)

        #Append formatted string to string list
        str_lst.append(line) 

    extra_str = ['! NumCom, NumEqu, NumPhases, binary_diffusion',
                '! coordinate_system, Max_NumElem, Max_NumConx, ElemNameLength, active_conx_only, boundaries_in_matrix',
                '! MaxNum_SS',
                '! MaxNum_Media',
                '! element_by_element_properties, porosity_perm_dependence, scaled_capillary_pressure, Option_tortuosity_CompuMethod',
                "! coupled_geochemistry, property_update [= 'Continuous', 'Iteration', 'Timestep']",
                '! coupled_geomechanics, geomechanical_code_name, property_update, num_geomech_param']    

    final_str_lst = str_lst[0:2]

    for item, text in zip(str_lst[2:],extra_str):
        final_str_lst.append(item + ' '*(48-len(item)) + text)

    return final_str_lst  

def write_HYDRATE(hydrate):
    
    fmt_lst=['{:5d}',
            '{:5}{:9.1E}{:10.2E}',
            '{:5d}',
            '{:8.1E}',
            '{:5d}',
            '{:8.1E}',
            '{:5d}',
            '{:8.1E}',
            '.{:s}.{:9.1E}{:11.3E}{:12.4E}{:9.1E}{:12.4E}{:12.4E}{:12.4E}{:12.4E}',
            '{:5d}',
            '{:s}',
            '{:10.2E}{:13.5E}{:10.2E}'
            ]


    fmts=dict(zip(hydrate.keys(),fmt_lst))

    if hydrate['H_7']['N_Rho'] == 0:
        print("hydrate density calculated from Ballard [2002]")
        hydrate.pop('H_8', None)
        fmts.pop('H_8', None)


    str_lst=['HYDRATE']

        
    #iterate over dictionary
    for key,vals in hydrate.items():
        val_lst=[]

        subvals=list(hydrate[key].values())

        if 'EQUILIBRIUM' in str(subvals[0]):
            val_lst.append(subvals[0])
            line=fmts[key].format(*val_lst)
            str_lst.append(line)
            break

        #Iterate over each sub dictionary and extract values
        for val in subvals:

            #Process iterables
            if type(val)==list:
                val_lst+=val


            else:
                #Check for flags
                if type(val)==bool:

                    flag = val
                    val=str(val).upper()
                    if not flag:
                        val_lst.append(val)
                        fmts[key]=fmts[key][:6]
                        #Break loop if there inhibitor flag is False
                        print('inhibitor_flag is',val,'Other values will be ignored.')
                        break


                #Update list of iterables
                val_lst.append(val)

        #Format iterable
        line=fmts[key].format(*val_lst)

        #Append formatted string to string list
        str_lst.append(line)    
    
    return str_lst
    

def write_ROCKS(table):
    
    if table.columns[0]=='Type':
        str_table=table.iloc[:,1:].copy()
    else:
        str_table=table.copy()
    
    str_table.loc[str_table['Perm2'].isnull(),'Perm2']=str_table['Perm1'][str_table['Perm2'].isnull()]
    str_table.loc[str_table['Perm3'].isnull(),'Perm3']=str_table['Perm1'][str_table['Perm3'].isnull()]
    
    fmts=['{:5}','{:5}']+8*['{:10.3E}'] \
        + 11*['{:10.3E}'] \
        +[('{:5}'+5*' ')]+7*['{:10.3E}'] \
        +[('{:5}'+5*' ')]+7*['{:10.3E}'] \
        +[('{:5}'+5*' ')]+7*['{:20.12E}'] \
        +5*['{:10.3E}']
        
    for col_i, col_n in enumerate(str_table.columns):
        # print(col_n, str_table.iloc[0,col_i], fmts[col_i])
        str_table[col_n]=str_table[col_n].map(fmts[col_i].format)
            
    
    str_rocks='ROCKS----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n'
    
    for i,row in str_table.iterrows():
        a=0
        str_rocks+=''.join(row.iloc[a:a+10].values.astype(str))+'\n'
        a+=10
        
        if table.iloc[i]['NAD']>=1:
            str_rocks+=''.join(row.iloc[a:a+11].values.astype(str))+'\n'
            a+=11
        if table.iloc[i]['NAD']>=2:
            str_rocks+=''.join(row.iloc[a:a+8].values.astype(str))+'\n'
            a+=8
            str_rocks+=''.join(row.iloc[a:a+8].values.astype(str))+'\n'
            a+=8
        if table.iloc[i]['NAD']==5:
            str_rocks+=''.join(row.iloc[a:a+8].values.astype(str))+'\n'
            a+=8
        if table.iloc[i]['NAD']==6:
            str_rocks+=''.join(row.iloc[a:a+5].values.astype(str))+'\n'            
    
    str_rocks=str_rocks.replace('NAN','   ')
    
    str_rocks=str_rocks.split('\n')

    return(str_rocks)
    

def write_RPCAP(rpcap):
    rpcap_str = 'RPCAP----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n'

    fmts = [('{:5}'+5*' ')]+7*['{:10.3E}']

    str_table = rpcap.copy()

    for col,fmt in zip(str_table.iteritems(),fmts):
        str_table[col[0]] = col[1].map(fmt.format)

    str_rpcap=['RPCAP----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8']
        
    for row in str_table.iterrows():
        line = ''.join(row[1].values)
        line = line.replace('NAN','   ')
        str_rpcap += [line]
        
    str_rpcap += ['\n']    

    
    return str_rpcap                
    
    


def write_PARAM(param):

    #Make a copy of input DF/series
    str_param=param.copy()

    #Replace nan values with blanks
    str_param=str_param.fillna('')


    def param_line(line,fmt,fmt_blanks):
        idx=0
        param_l=''
        for item,value in line.items():

            if isinstance(value, str):

                str_value=fmt_blanks[idx].format(value)
                param_l+=str_value
            else:

                str_value=fmt[idx].format(value)
                param_l+= str_value

            
            idx+=1

        return param_l

    #PARAM.1
    fmt_1=2*['{:2d}']+ 3*['{:4d}']+ 1*['{:24}']+ 3*['{:10.3E}'] + 2*['{:5d}']
    fmt_1_blanks=2*['{:>2}']+ 3*['{:>4}']+ 1*['{:>24}']+ 3*['{:>10}'] + 2*['{:>5}']
    str_param['P_1'].MOP=''.join(str_param['P_1'].MOP_items.values())

    param_1=param_line(str_param.P_1.iloc[:-1],fmt_1,fmt_1_blanks)

    #PARAM.2
    fmt_2=4*['{:10.3E}']+2*['{:>5s}'] + 3*['{:10.3E}']
    fmt_2_blanks=4*['{:>10}']+2*['{:>5s}'] + 3*['{:>10}']
    param_2=param_line(str_param.P_2,fmt_2,fmt_2_blanks)



    #PARAM.3
    fmt_3=7*['{:10.3E}']+['{:>2s}'] + ['{:>3s}'] + ['{:>5s}'] + 3*['{:10.3E}']
    fmt_3_blanks=7*['{:>10}']+['{:>2s}'] + ['{:>3s}'] + ['{:>5s}'] + 3*['{:>10}']
    param_3=param_line(str_param.P_3,fmt_3,fmt_3_blanks)


    #PARAM.4
    fmt_4=6*['{:20.5E}']
    fmt_4_blanks=6*['{:>20}']
    param_4=param_line(str_param.P_4,fmt_4,fmt_4_blanks)

    param_lst=['PARAM----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8----*----9----*---10',
            param_1,param_2,param_3,param_4,'\n'] 




    return param_lst

def write_INCON(incon):
    incon_header='INCON----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8----*----9----*---10'

    try:
        str_incon=incon.copy()
        
        
        
        col_names=['ElName', 'NSEQ', 'NADD', 'porosity', 'StateIndex', 'comment',
                'PermX', 'PermY', 'PermZ', 'phi', 'var1', 'var2', 'var3', 'var4',
                'var5', 'var6']
        col_format=['{:5s}']+2*['{:5.0f}']+['{:15.6E}  ']+['{:3s}']+['{:36s}']+4*['{:15.6E}']+6*['{:20.5E}']
        col_blanks=['{:5s}']+2*['{:5s}']+['{:15s}  ']+['{:3s}']+['{:36s}']+4*['{:15s}']+6*['{:20s}']

        
        for col, fmt, blk in zip(col_names,col_format,col_blanks):
            # print('formatting', col)
            try:
                str_incon[col]=incon[col].map(fmt.format)
            except:
                str_incon[col]=incon[col].map(blk.format)

        
        cols_1=col_names[:10]
        cols_2=col_names[10:]
        
        
        str_incon['raw1']=str_incon[cols_1].apply(lambda row: ''.join(row.values), axis=1)
        str_incon['raw1'] = str_incon['raw1'].str.replace('nan', '   ')
        str_incon['raw1'] = str_incon['raw1'].str.replace('NAN', '   ')

        str_incon['raw2']=str_incon[cols_2].apply(lambda row: ''.join(row.values), axis=1)
        str_incon['raw2'] = str_incon['raw2'].str.replace('nan', '   ')
        str_incon['raw2'] = str_incon['raw2'].str.replace('NAN', '   ')


        raw_incon=[None]*2*len(str_incon)
        raw_incon[::2]=str_incon['raw1'].values
        raw_incon[1::2]=str_incon['raw2'].values

        raw_incon=[incon_header]+raw_incon+['\n']

    
    except:
        raw_incon=[incon_header]+['\n']

        pass

    return raw_incon


def write_CONNE(conne):

    col_format=2*['{:5s}']+4*['{:5d}']+5*['{:10.3E}']
    col_blanks=2*['{:5s}']+4*['{:5s}']+5*['{:10s}']
    col_names = ['ConxName1', 'ConxName2', 'NSEQ', 'NAD1', 'NAD2', 'ConxKi', 'ConxD1', 'ConxD2', 'ConxArea', 'ConxBeta', 'emissivity']

    str_conne = conne.copy()
    str_conne = str_conne.fillna('')

    for col, fmt, blk in zip(col_names,col_format,col_blanks):
        # print('formatting', col)
        try:
            str_conne[col]=str_conne[col].map(fmt.format)
        except:
            str_conne[col]=str_conne[col].map(blk.format)
            
            
    str_conne['raw']=str_conne[col_names].apply(lambda row: ''.join(row.values), axis=1)

    raw_conne = str_conne['raw'].to_list()

    conne_header = 'CONNE'

    raw_conne=[conne_header]+raw_conne+['\n']

    return raw_conne

def write_SUBDOMAINS(subdomains):
    import pandas as pd
    SD_list=['SUBDOMAINS']
    heading='&Subdomain_General_Info  {:s} = {:d} /'
             
    subdm_head=[ 3*' '+'&Individual_Subdomain_Specifics  {:s} = {:s},',
            ' '*36+'{:s} = {:d}\n'+' '*36+'/']
    
    def print_region(reg_series):
        first_fmt=6*' '+'&Region_Specifics  {:s} = {:s},'
        lines_fmt=25*' '+'{:s} = {:s},'
        last_fmt=lines_fmt+'\n'+25*' '+'/'
        
        reg_str=[]
        
        
        
        line_ct=0
        for idx,item in reg_series.iteritems():
            if isinstance(item,float):
                str_item='{:.2e}'.format(item)
            elif isinstance(item,int):
                str_item='{:d}'.format(item)
            else:
                str_item=item
            
            
            if line_ct==0:
                reg_str.append(first_fmt.format(idx,str_item))
            elif idx=='elements':
                continue
            elif line_ct==len(reg_series)-1:
                reg_str.append(last_fmt.format(idx,str_item))

            else:
                reg_str.append(lines_fmt.format(idx,str_item))
            line_ct+=1
            
        if 'elements' in reg_series:
            reg_str.append(', '.join(map(str,reg_series['elements'])))
                           
        return reg_str
            
    
    for key1,item1 in subdomains.items():
        if isinstance(item1,dict):
            subdomain=item1
            for key2,item2 in subdomain.items():
                if isinstance(item2,pd.DataFrame):
                    for idx,row in item2.iterrows():
                        if row.definition_mode=="'Geometry'":
                            dft_idx=[0,6]
                            if row.region_shape=="'Rectangle'":
                                f_row=row.iloc[dft_idx+[8,9,10,11,12,13]]
                                
                            elif row.region_shape=="'Cylinder'":
                                f_row=row.iloc[dft_idx+[12,13,14,15]]
                            elif row.region_shape=="'Sphere'":
                                f_row=row.iloc[dft_idx+[14,15,18]]
                            
                            else:
                                print('Wrong definition of geometry')
                            
                            
                                
                        elif row.definition_mode=="'Sequence'":
                            #retrieve mode of defining
                            idx_str=row.iloc[3:5][(~row.iloc[2:5].isna())].index
                            idx_no=list(row.index).index(idx_str)
                            dft_idx=[0,1,idx_no,5]
                            
                            f_row=row.iloc[dft_idx]
                            
                        elif row.definition_mode=="'NumberList'" or row.definition_mode=="'NameList'":
                            dft_idx=[0,1,2,7]
                            f_row=row.iloc[dft_idx]

                        else:
                            print("wrong definition name")
                
                    str_row=print_region(f_row)
                    SD_list+=str_row
                
                else:
                    line_fmt=list(subdomain.keys()).index(key2)
                    SD_list.append(subdm_head[line_fmt].format(key2,item2))
                    

        else:
            SD_list.append(heading.format(key1,item1))
        SD_list.append('\n')

        
    return SD_list

def write_GENER(gener):

    def GENER_record(ElName5C,SS_name,NSEQ='',NADD='',NADS='',LTAB='',SS_Type='',ITAB='',
          GX='', EX='', HX='',WellResponse='', PresLimits='',RateStepChange='', 
          RateLimit='',table=''):
        """
        Format (A5, A5, 4I5, 5X,A4, A1, 3E10.4 ,A4, 6x, 3(E10.4))
        """
        gen_line_1=''

        #A5
        gen_line_1+='{:5}'.format(ElName5C)

        #A5
        gen_line_1+='{:5}'.format(SS_name)

        #4I5
        try:
            gen_line_1+='{:<5d}'.format(NSEQ)
        except:
            gen_line_1+='{:5}'.format(NSEQ)

        try:
            gen_line_1+='{:<5d}'.format(NADD)
        except:
            gen_line_1+='{:5}'.format(NADD)
            
        try:
            gen_line_1+='{:<5d}'.format(NADS)
        except:
            gen_line_1+='{:5}'.format(NADS)
        
        try:
            gen_line_1+='{:<5d}'.format(LTAB)
        except:
            gen_line_1+='{:5}'.format(LTAB)
            
        #5X
        gen_line_1+='{:5}'.format('')
        
        #A4
        gen_line_1+='{:4}'.format(SS_Type)
        
        #A1
        gen_line_1+='{:1}'.format(ITAB)
        
        #3E10.4 
        try:
            gen_line_1+='{:10.2e}'.format(GX)
        except:
            gen_line_1+='{:10}'.format(GX)
        
        try:
            gen_line_1+='{:10.3e}'.format(EX)
        except:
            gen_line_1+='{:10}'.format(EX)

        try:
            gen_line_1+='{:10.3e}'.format(HX)
        except:
            gen_line_1+='{:10}'.format(HX)
        
        #A4    
        gen_line_1+='{:4}'.format(WellResponse)
        
        #6x
        gen_line_1+='{:6}'.format('')

        #3(E10.4)
        try:
            gen_line_1+='{:10.3e}'.format(PresLimits)
        except:
            gen_line_1+='{:10}'.format(PresLimits)
        
        try:
            gen_line_1+='{:10.3e}'.format(RateStepChange)
        except:
            gen_line_1+='{:10}'.format(RateStepChange)

        try:
            gen_line_1+='{:10.3e}'.format(RateLimit)
        except:
            gen_line_1+='{:10}'.format(RateLimit)
        
        
        if LTAB>1:
            gen_line_2='\n'

            for line in table:
                line_start = 0
                for j in range(1+(len(line)-1)//4):
                    line_end = (j+1)*4
                    subline = line[line_start:line_end]
                    gen_line_2+=('{:14.5E}'*len(subline)).format(*subline)+'\n'
                    line_start = line_end
            
            gen_line_2=gen_line_2[:-1]
            gen_line = gen_line_1+gen_line_2
        else:
            gen_line = gen_line_1
            
        gen_line = gen_line.replace('nan', '   ')
        
        return gen_line

    gener_str = ['GENER']

    for row in gener.iterrows():
        gener_str.append(GENER_record(*row[1]))

    gener_str.append('\n')
    
    return gener_str


def MESH_file_to_str(file):
    

    with open(file) as mesh_f:
        mesh_data=mesh_f.read()


    conne_idx=mesh_data.find('CONNE')

    eleme=mesh_data[:conne_idx-1]
    conne=mesh_data[conne_idx:-1]


    eleme=eleme.split('\n')
    conne=conne.split('\n')

    return eleme,conne


#%%Read DATA INPUT T+H file
def read_TH_data(file):
    import pandas as pd
    import re
    import io
    KEYWORDS=["TITLE","MEMORY","ROCKS","RPCAP","HYDRATE","START","PARAM","MEDIA","RANDO","WETTABILITY","DIFFUSION",
              "ELEME","CONNE","GENER","INDOM","INCON","EXT-INCON",">>>BOUNDARIES","SOLVR",
              "TIMES","SUBDOMAINS","INTERFACES","SS_GROUPS","ENDCY","ENDFI"]
    
    

    
    
    # data=dict.fromkeys(KEYWORDS,dict())
    
    data={keyword:{'active':False} for keyword in KEYWORDS}
    
    
    proc_func={'MEMORY':read_MEMORY, 'HYDRATE':read_HYDRATE, "ROCKS":read_ROCKS, "RPCAP":read_RPCAP, 'INCON':read_INCON, 'CONNE':read_CONNE, 'ELEME':read_ELEME, 'INDOM':read_INDOM,
               'PARAM':read_PARAM}
    
    sections_idx=pd.DataFrame(index=KEYWORDS,columns=['starts','ends'])
    
    sections_idx.loc['TITLE']=[0,1]
    
    #Read file
    if isinstance(file,io.StringIO):
        text=file.getvalue()
    
    else:
        with open(file,'r',encoding='cp1250') as f:
            text = f.read()
    
    #Split file in lines
    text = text.split('\n')
    
    # str_iter=iter(text)
    

    
    for i,line in enumerate(text):
        
       
        #Make match pattern to identify lines
        pattern=re.compile(r"(\b{}\b)".format("|".join(KEYWORDS[1:])))
        
        match_line=pattern.match(line)
        
       
        if match_line:
            keyword=match_line.group(1)
            sections_idx.loc[keyword]['starts']=i

        if line.startswith('ENDCY') or line.startswith('ENDFI'):
            end_line=i+1
            break        
        
    sections_idx=sections_idx.dropna(how='all')
    sections_idx=sections_idx.sort_values(by='starts')
    sections_idx['ends'].iloc[1:-1]=sections_idx['starts'].iloc[2:].values
    sections_idx['ends'].iloc[-1]=end_line
    
    for section,vals in sections_idx.iterrows():
        data[section]['raw']=text[vals[0]:vals[1]]
        data[section]['active']=True
        if section == 'MEMORY':
            data[section]['modified']=True
        else:
            data[section]['modified']=False
        
        
        
        if section in proc_func.keys():
            data[section]['processed']=proc_func[section](data[section]['raw'])
            if section=='ELEME' and 'ina' in data[section]['raw'] :
                print('ELEME section in old format. ELEME section stored will be updated')
                updated_raw=data[section]['processed']['raw'].tolist()
                updated_raw=['ELEME']+updated_raw+['']
                data[section]['raw']=updated_raw

    data=pd.DataFrame(data)
 
    return data


def update_TH_data(data,mesh_file=False,mesh_path=''):

    if mesh_file:
        eleme,conne=MESH_file_to_str(mesh_path)
        data['ELEME']['raw'] = eleme
        data['ELEME']['processed'] = read_ELEME(eleme)
        data['ELEME']['active'] = True

        data['CONNE']['raw'] = conne
        data['CONNE']['processed'] = read_ELEME(conne)
        data['CONNE']['active'] = True


    data = check_MEMORY(data)

    updt_func={"MEMORY":write_MEMORY, "HYDRATE":write_HYDRATE, "ROCKS":write_ROCKS, 
                "RPCAP":write_RPCAP, "PARAM":write_PARAM, "CONNE":write_CONNE, 
                "INCON":write_INCON, "SUBDOMAINS":write_SUBDOMAINS, "GENER":write_GENER}
    
    for section in data:
    
        if data[section]['active']:
            if data[section]['modified']:
                if section in section in updt_func.keys():
                    print("Updating with function", section)
                    data[section]['raw']=updt_func[section](data[section]['processed'])
                else:
                    print("Writing function not defined for section ",section)
    


    return data


def write_TH_data(data,file,ow=False):
    import os
    import pathlib
    
    list_data=[]
    
    for key in data:
        if data[key]['active']:
            list_data+=data[key]['raw']
    
    str_data='\n'.join(list_data)
    
    
    if os.path.isfile(file):

        if ow:
            print(file,'already exists and will be overwritten')
            os.remove(file)
        else:
            print(file,'already exists and cannot be overwritten')
            return
    else:
        print(file,'does not exists and will be created')
        dirname=os.path.dirname(file)
        if len(dirname)==0:
            print(file,"created in ",os.getcwd())
        else:
            if os.path.isdir(dirname):
                print(file,'created in ',dirname)
            else:
                path=pathlib.Path(file)
                path.parent.mkdir(parents=True,exist_ok=True)
                print('creating',dirname,'directory')
    
    
    
    with open(file,'x') as f:
        f.write(str_data)
        
        
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
        
        
