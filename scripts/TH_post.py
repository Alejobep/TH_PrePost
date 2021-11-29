from TH_proc import *
import itertools
import numpy as np
import pandas as pd
import os
import io

#%%Auxiliar functions
def read_in_chunks(f,chunk_size):
    while True:
        chunk = f.read(chunk_size)
        if not chunk:
            break
        yield chunk

#%%Processing functions
def read_Time_Series(path):

    old_cwd = os.getcwd()

    if os.path.isdir(path):
        folder = path
    else:
        folder = os.path.dirname(path)

    os.chdir(folder)


    interfs=dict()
    subdoms=dict()
    ss_groups=dict()


    for f in os.listdir():
        if r'Time_Series' in f or r'Hydrate_Status' in f:
            with open(f) as ts_file:
                first_line = ts_file.readline()


            ts_df = pd.read_fwf(f,skiprows = 1)

            name = f[:-12]

            if 'Subdomain' in first_line:
                subdoms[name] = ts_df

            elif 'Interface' in first_line:
                interfs[name] = ts_df

            elif 'SSGroup' in first_line:
                ss_groups[name] = ts_df

            else:
                hydrate_status = pd.read_fwf(f)


    print('{:d} subdomain(s)\n{:d} interface(s)\n1 Hydrate_Status file\n{:d} source/sink group(s)'.format(len(subdoms.keys()), len(interfs.keys()), len(ss_groups.keys())))
    os.chdir(old_cwd)

    return subdoms, interfs, hydrate_status, ss_groups
            


    
    
    

def ELEME_geometry(elem_data):
    #Remove raw column
    elem_data = elem_data.iloc[:,1:]

    elem_data=elem_data.sort_values(by=['X','Y','Z'],ascending=[True,True,False]) #Sort by X and Z coordinates
    

    Xc_val=elem_data.X.drop_duplicates().values
    Yc_val=elem_data.Y.drop_duplicates().values
    Zc_val=elem_data.Z.drop_duplicates().values

    Tdim=elem_data.shape[0]
    Xdim=Xc_val.shape[0] #Calculate amount of grid cells in X direction
    Ydim=Yc_val.shape[0] #Calculate amount of grid cells in X direction
    Zdim=Zc_val.shape[0] #Calculate amount of grid cells in Z direction

    if Xdim*Ydim*Zdim>Tdim:
        print("Grid with unactive cells")

    print(Xdim*Ydim*Zdim,' cells. ',Tdim,' cells are active.')


    all_coord=list(itertools.product(Xc_val,Yc_val,Zc_val)) #All possible coordinates
    act_coord=np.array(elem_data[['X','Y','Z']]) #Coordinates defined in MESH
    act_coord=[tuple(coord) for coord in act_coord]

    inact_coord=list(set(all_coord)-set(act_coord))

    inact_coord=pd.DataFrame(data=inact_coord,columns=['X','Y','Z'])
    inact_coord['active']=False
    inact_coord['MA12']='INA'

    elem_data=elem_data.append(inact_coord, ignore_index=True,sort=False)

    elem_data=elem_data.sort_values(by=['X','Y','Z'],ascending=[True,True,False])
    elem_data=elem_data.reset_index(drop=True)

    elem_data.loc[:,'elem_vol']=elem_data.loc[:,'elem_vol'].interpolate(method='linear',limit_direction='both')


    Xcoord=np.array(elem_data.X).reshape(Xdim,Ydim,Zdim) #X coordinates in mxn matrix form
    Ycoord=np.array(elem_data.Y).reshape(Xdim,Ydim,Zdim) #X coordinates in mxn matrix form
    Zcoord=np.array(elem_data.Z).reshape(Xdim,Ydim,Zdim) #Z coordinates in mxn matrix form


    vol=np.array(elem_data.elem_vol).reshape(Xdim,Ydim,Zdim) #grid cell volume in mxn matrix form

    #Create arrays for using in pcolormesh
    pmesh_dim=(Xdim+1,Ydim+1,Zdim+1)
    X=np.zeros(pmesh_dim) #Empty X matrix to build mesh for plotting
    Z=np.zeros(pmesh_dim) #Empty Z matrix to build mesh for plotting
    Y=np.zeros(pmesh_dim) #Empty Y matrix to build mesh for plotting

    Xsize=np.zeros(Xcoord.shape)
    Ysize=np.zeros(Xcoord.shape)
    Zsize=np.zeros(Xcoord.shape)


        
    delta=0
    corner=0
    for i,coord in enumerate(Xcoord[:,0,0]):
        delta=np.round((coord-corner)*2,5)
        # print(coord,corner,delta)
        corner=np.round(corner+delta,5)
        X[i+1,:,:]=corner
        Xsize[i,:,:]=delta

    delta=0
    corner=0
    for i,coord in enumerate(Ycoord[0,:,0]):
        delta=(coord-corner)*2
        corner+=delta
        Y[:,i+1,:]=corner
        Ysize[:,i,:]=delta


    elem_data['dX']= Xsize.flatten()
    elem_data['dY']= Ysize.flatten()

    shallowest=elem_data.loc[elem_data['Z'].idxmax()]
    dZ_sh=shallowest.elem_vol/(shallowest.dX*shallowest.dY)

    delta=0
    corner=shallowest.Z+dZ_sh/2
    Z+=corner
    for i,coord in enumerate(Zcoord[0,0,:]):
        delta=(coord-corner)*2
        corner+=delta
        Z[:,:,i+1]=corner

    Zsize=vol/(Ysize*Xsize)
    elem_data['dZ']= Zsize.flatten()


    x_c=np.unique(elem_data.X.values)
    y_c=np.unique(elem_data.Y.values)
    z_c=np.unique(elem_data.Z.values)

    x_c=np.sort(x_c)
    y_c=np.sort(y_c)
    z_c=np.sort(z_c)[::-1]

    I=pd.Series(index=elem_data.index,dtype=int)
    J=pd.Series(index=elem_data.index,dtype=int)
    K=pd.Series(index=elem_data.index,dtype=int)

    for i,X_val in enumerate(x_c):
        I[(elem_data['X']==X_val)]=i

    for j,Y_val in enumerate(y_c):
        J[(elem_data['Y']==Y_val)]=j
        
    for k,Z_val in enumerate(z_c):
        K[(elem_data['Z']==Z_val)]=k
                
    elem_data['I'],elem_data['J'],elem_data['K']=I,J,K

    return elem_data
    # return elem_data[~(elem_data['active']==False)].copy()

def process_init(file_path):

    import os



    data = read_TH_data(file_path)

    eleme_data = data.ELEME.processed.copy()

    mesh = ELEME_geometry(eleme_data)


    #Set reaction type flag
    if 'KINETIC' in data.HYDRATE.processed['H_11']['Reaction_Type']:
        equil_flag = False
    else:
        equil_flag = True

    #Set inhibitor flag
    inh_flag = data.HYDRATE.processed['H_9']['inhibitor_flag']

    #Create Default State Index table
    State_Id=["Gas","Aqu","AqG","IcG","GsH","AqH","AqI","IcH","AGH","AIG","AIH","IGH","QuP"]
    
    if equil_flag:
        State_var=['phases','var1','var2','var3','var4','n_phases','phases']
        StateIndex=pd.DataFrame(index=State_Id,columns=State_var,dtype=str)
        
        
        # StateIndex.var1=["P_gas","P","P_gas","P_gas","P_gas","P","P","P","S_gas","P_gas","P","S_gas","S_gas"]
        StateIndex.var1=["P","P","P","P","P","P","P","P","S_gas","P","P","S_gas","S_gas"]
        StateIndex.var2=["Y_m_G","X_m_A","S_aqu","S_ice","S_gas","S_aqu","S_aqu","S_ice","S_aqu","S_aqu","S_aqu","S_ice","S_aqu"]
        StateIndex.var3='X_inh'
        StateIndex.var4=["T","T","T","T","T","T","X_m_A","T","T","S_gas","S_ice","T","S_ice"]
        StateIndex.n_phases=2*[1]+6*[2]+4*[3]+[4]
        StateIndex.phases=["G","A","AG","IG","HG","AH","AI","IH","AH","AIG","AIH","IHG","IHAG"]
    else:
        State_var=['phases','var1','var2','var3','var4','var5','n_phases','phases']
        StateIndex=pd.DataFrame(index=State_Id,columns=State_var,dtype=str)
        
        
        # StateIndex.var1=["P_gas","P","P_gas","P_gas","P_gas","P","P","P","S_gas","P_gas","P","S_gas","S_gas"]
        StateIndex.var1=["P","P","P","P","P","P","P","P","P","P","P","P","P"]
        StateIndex.var2=['Y_m_G','X_m_A','S_aqu','S_ice','S_gas','S_aqu','S_aqu','S_ice','S_aqu','S_aqu','S_aqu','S_gas','S_aqu']
        StateIndex.var3=['S_hyd','S_hyd','S_hyd','S_hyd','S_ice','X_m_A','X_m_A','S_gas','S_gas','S_hyd','S_ice','S_ice','S_gas']
        StateIndex.var4='X_inh'
        StateIndex.var5=['T','T','T','T','T','T','T','T','T','S_gas','T','T','S_ice']
        StateIndex.n_phases=2*[1]+6*[2]+4*[3]+[4]
        StateIndex.phases=["G","A","AG","IG","HG","AH","AI","IH","AH","AIG","AIH","IHG","IHAG"]

    
    
    #Table of State Index
    table_Vars=StateIndex.iloc[:,1:-2]

    #Extract default State Index from DATA file
    default_SI = data.PARAM.processed.P_3.DefaultStateIndex
    default_IC = data.PARAM.processed.P_4.astype(float)

    #Set Default StateIndex into mesh
    mesh['StateIdx'] = default_SI
    mesh['IC_var1']=default_IC.var1
    mesh['IC_var2']=default_IC.var2
    mesh['IC_var3']=default_IC.var3

    if equil_flag:
        if inh_flag:
            mesh['IC_var4']=default_IC.var4
    else:
        mesh['IC_var4']=default_IC.var4
        if inh_flag:
            mesh['IC_var5']=default_IC.var5
            
    #Set INCON StateINdex into mesh
    incon = data.INCON.processed

    if len(incon)==0:
        print('No explictit INCON file')
        print('Read default initial contidions from Input file')
        IC=False
    else:
        IC=True
        print('INCON element-specific initial conditions')

    
    if IC:
        mesh['old_IDX']=mesh.index
        mesh=mesh.set_index('ElName', drop=False)
        
        print("Processing element specific initial conditions")
        
        mesh.loc[incon.index,'StateIdx']=incon['StateIndex']
        mesh.loc[incon.index,'IC_var1']=incon['var1'].astype(float)
        mesh.loc[incon.index,'IC_var2']=incon['var2'].astype(float)
        mesh.loc[incon.index,'IC_var3']=incon['var3'].astype(float)
        
        if equil_flag:
            if inh_flag:
                 mesh.loc[incon.index,'IC_var4']=incon['var4'].astype(float)
                 
        else:
            mesh.loc[incon.index,'IC_var4']=incon['var4'].astype(float)
            if inh_flag:
                mesh.loc[incon.index,'IC_var5']=incon['var5'].astype(float)
        
        mesh=mesh.set_index(mesh['old_IDX'].values)        

    mesh = Set_IC_Vars(mesh,inh_flag,equil_flag)

        
    return mesh


def readPDE_file(path):
    from sys import exit
    #Extract variable header lines
    # print("Get variable header for "+path)

    if os.path.isdir(path):
        case = path
    else:
        case = os.path.dirname(path)

        if path.endswith(r'.in'):
            ip_file = path
        else:
            files = os.listdir(case)
            for f in files:
                if f.endswith(r'.in'):
                    ip_file = os.path.join(case,f)
                    

    if r'Plot_Data_Elem' not in path:
        path = os.path.join(case,r'Plot_Data_Elem')


        
    with open(path) as f:
        var_header=f.readline()
        var_header=var_header.strip().split()[2:]


    #Get dimenstions of Data Structure = time_steps x cells x variables
    str_line=44 #Top line of each Time Step ZONE T= .....
    rec_width=15 #Width of the column in Plot_Data_Elem file
    n_var=len(var_header)
   
    ip_data = read_TH_data(ip_file)
    eleme = ip_data.ELEME.processed
    
    eleme = eleme.sort_values(by=['X','Y','Z'],ascending=[True,True,False])
    Tdim = len(eleme)

    chunk_size=str_line+Tdim*n_var*rec_width+3
    # print("Get number of time steps")
    n_time_steps=0
    with open(path) as f:
        f.readline()
        for chunk in read_in_chunks(f,chunk_size):
            n_time_steps+=1

    # print(n_time_steps)



    df = pd.read_fwf(path, names = var_header, skiprows = 1, widths = [14]+[15]*(len(var_header)-1))

    #Create numpy array to index lines with time step info and lines to be dropped
    drop_lines = np.arange(0,len(df),Tdim+4)
    tstep_lines = df.loc[drop_lines].copy()

    drop_lines = np.concatenate((drop_lines, np.arange(Tdim+1,len(df),Tdim+4)))
    drop_lines = np.concatenate((drop_lines, np.arange(Tdim+2,len(df),Tdim+4)))
    drop_lines = np.concatenate((drop_lines, np.arange(Tdim+3,len(df),Tdim+4)))

    drop_lines.sort()

    df = df.drop(index=drop_lines)
    df = df.apply(pd.to_numeric)

    #Process time step lines
    tstep_lines = tstep_lines.astype(str)

    tstep_lines = tstep_lines.apply(''.join, axis=1)

    tstep_lines = tstep_lines.str.strip('nan')

    tstep_lines = tstep_lines.str[9:24]

    tsteps = pd.to_numeric(tstep_lines).values

    #Create multi index
    time = np.repeat(tsteps, Tdim)
    ElName = np.tile(eleme.ElName.values, len(tsteps))
    MA12 = np.tile(eleme.MA12, len(tsteps))

    df['time'] = time
    df = df.sort_values(by=['time','x','y','z'], ascending=[True,True,True,False])
    df['ElName'] = ElName
    df['MA12'] = MA12

    df = df.set_index(['time', 'ElName'])

    return df










def get_output(*file, write=False):
    """
    Reads initialization parameters
    Process Plot_Data_Elem
    Returns a dataframe with both initialization and input
    """


    if len(file)==1:
        pklpath = os.path.dirname(file[0])

    else:
        pklpath = os.path.dirname(os.path.dirname(file[0]))


    pklfile = os.path.join(pklpath,'GRID_data.pkl')

    pkl_flag = False
    if os.path.isfile(pklfile):
        pkl_flag = True
        
    
    if pkl_flag:
        print('Read pickle data')
        full_tdata = pd.read_pickle(pklfile)

    else:
        print('No pickle data available. Read and process')

        if len(file)==1:
            print("Processing single file")

            data = read_TH_data(file[0])
            mesh = process_init(file[0])
            tdata = readPDE_file(file[0])

            pklpath = os.path.dirname(file[0])

            

        else:
            print("Processing multiple files")        
            n_files = len(file)
            sorted_files = pd.DataFrame(index = range(n_files), columns=('Name','t_min','ctime'))
            
            sorted_files['Name']=file
            sorted_files['ctime']=sorted_files.Name.map(os.path.getctime)

            sorted_files = sorted_files.sort_values(by='Name')
        
            
            for idx,row in sorted_files.iterrows():
                
                try:
                    tdata = readPDE_file(row.Name)

                    min_t = tdata.index.get_level_values(0).min()
                    
                    sorted_files.loc[idx,'t_min'] = min_t
                    try:
                        merged = pd.concat([merged,tdata], sort = False)
                    
                    except:
                        merged = pd.concat([tdata], sort = False)

                    # print(idx, row.Name)
                    last_file = row.Name

            
                except:
                    print('Last file processed:',last_file)
                    continue

            
            tdata = merged
                    
            sorted_files = sorted_files.sort_values(by='t_min')
            
            sorted_files = sorted_files.reset_index()   

            # for f_idx, f_row in sorted_files.iterrows():
            #     print('iteration No. {:d}. {:s}'.format(f_idx,f_row['Name']))


            file_0 = sorted_files.loc[0,'Name']
            mesh = process_init(file_0)
            data = read_TH_data(file_0)

            pklpath = os.path.dirname(os.path.dirname(file[0]))


        tsteps = tdata.index.get_level_values(0).drop_duplicates()
        t1 = tsteps[0]

        print(len(tsteps),'steps')

        tdata_0 = tdata.loc[t1].copy()
        tdata_0.iloc[:,3:] = np.nan

        #Set porosity and permeability from ROCKS
        rocks = data.ROCKS.processed
        rocks = rocks.set_index('Name')

        #Reindex INCON mesh to ELNAme
        mesh = mesh.set_index('ElName')
        
        try:
            eleme_rock_idx = pd.to_numeric(mesh.MA12)
            print('Rocks referenced in ELEME as indexes')
            tdata_0.porosity = rocks.iloc[1-eleme_rock_idx.values]['Poros'].values
            tdata_0.perm_abs = rocks.iloc[1-eleme_rock_idx.values]['Poros'].values
            
        except:
            print('Rocks referenced as string names')
            tdata_0.porosity = rocks.loc[mesh.MA12,'Poros'].values
            tdata_0.perm_abs = rocks.loc[mesh.MA12,'Perm1'].values 

        tdata_0.loc[:,['MA12','P','T','S_hyd','S_gas','S_aqu','S_icd','X_inh']]=mesh[['MA12','P','T','S_hyd','S_gas','S_aqu','S_icd','X_inh']]

        #Append 1st time step and create new multi-Index
        tsteps = np.append(0,tsteps)

        MI_iterables = [tsteps,tdata_0.index]
        df_MIndex = pd.MultiIndex.from_product(MI_iterables, names=['time', 'ElName'])

        #Concatenate

        full_tdata = pd.concat([tdata_0,tdata])
        full_tdata = full_tdata.set_index(df_MIndex)

        full_tdata['I'] = np.tile(mesh.I.values, len(tsteps))
        full_tdata['J'] = np.tile(mesh.J.values, len(tsteps))
        full_tdata['K'] = np.tile(mesh.K.values, len(tsteps))
        full_tdata['dx'] = np.tile(mesh.dX.values, len(tsteps))
        full_tdata['dy'] = np.tile(mesh.dY.values, len(tsteps))
        full_tdata['dz'] = np.tile(mesh.dZ.values, len(tsteps))



        #Reorganize columns
        cols = full_tdata.columns.to_list()
        cols = cols[:3]+cols[17:]+cols[16:17]+cols[3:16]
        full_tdata = full_tdata[cols]

        #Change dtypes to float
        float_cols = cols[:3]+cols[10:]
        dtype_dict = dict(zip(float_cols,['float']*len(float_cols)))
        full_tdata = full_tdata.astype(dtype_dict)

        #Convert P from Pa to bar
        full_tdata['P'] /= 1e5
        full_tdata['P_cap'] /= 1e5


        #Calculate Phase Pressures
        full_tdata['P_aqu'] = full_tdata['P'] + full_tdata['P_cap']

        #Convert permeability from m2 to mD
        full_tdata['perm_abs'] *= 1013249965828144.8


        #Reset index to time only
        full_tdata = full_tdata.reset_index(level=[1])
        
        if write:
            pkl_file = os.path.join(pklpath,'GRID_data.pkl')
            print('storing data into pickle file: {:s}'.format(pkl_file))

            full_tdata.to_pickle(pkl_file)


    return full_tdata




def Set_IC_Vars(mesh,Inh,equil):
    from Aux_Functions import v_Equil_P    

    StIdxs=list(set(mesh.StateIdx.values))

    mesh['P'], mesh['T'], mesh['S_hyd'], mesh['S_aqu'], mesh['S_gas'], mesh['S_icd'], mesh['X_inh'] = [np.nan]*7

    
    for StIdx in StIdxs:
        indexes=(mesh['StateIdx']==StIdx)
        if equil:
            if StIdx == 'Gas':
                print('Gas: 1 phase Gas')
                mesh.loc[indexes,'P']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_aqu']=0
                mesh.loc[indexes,'S_gas']=1
                mesh.loc[indexes,'S_hyd']=0
                mesh.loc[indexes,'S_icd']=0
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var3']
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var4']
                else:
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var3']
            
            elif StIdx == 'Aqu':
                print('Aqu: 1 phase Aqueous')
                mesh.loc[indexes,'P']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_aqu']=1
                mesh.loc[indexes,'S_gas']=0
                mesh.loc[indexes,'S_hyd']=0
                mesh.loc[indexes,'S_icd']=0
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var3']
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var4']
                else:
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var3']
                    
                    
            elif StIdx == 'AqG':
                print('AqG: 2 phases Aqueous+Gas')
                mesh.loc[indexes,'P']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_aqu']=mesh.loc[indexes,'IC_var2']
                mesh.loc[indexes,'S_gas']=1-mesh.loc[indexes,'IC_var2']
                mesh.loc[indexes,'S_hyd']=0
                mesh.loc[indexes,'S_icd']=0
    
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var3']
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var4']
                else:
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var3']
                    
            elif StIdx == 'IcG':
                print('IcG: 2 phases Ice+Gas')
                mesh.loc[indexes,'P']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_aqu']=0
                mesh.loc[indexes,'S_gas']=1-mesh.loc[indexes,'IC_var2']
                mesh.loc[indexes,'S_hyd']=0
                mesh.loc[indexes,'S_icd']=mesh.loc[indexes,'IC_var2']
    
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var3']
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var4']
                else:
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var3']
                    
            elif StIdx == 'GsH':
                print('GsH: 2 phases Hydrate+Gas')
                mesh.loc[indexes,'P']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_aqu']=0
                mesh.loc[indexes,'S_gas']=mesh.loc[indexes,'IC_var2']
                mesh.loc[indexes,'S_hyd']=1-mesh.loc[indexes,'IC_var2']
                mesh.loc[indexes,'S_icd']=0
    
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var3']
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var4']
                else:
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var3']
        
    
            elif StIdx == 'AqH':
                print('AqH: 2 phases Aqueous+Hydrate')
                mesh.loc[indexes,'P']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_aqu']=mesh.loc[indexes,'IC_var2']
                mesh.loc[indexes,'S_gas']=0
                mesh.loc[indexes,'S_hyd']=1-mesh.loc[indexes,'IC_var2']
                mesh.loc[indexes,'S_icd']=0
    
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var3']
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var4']
                else:
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var3']
    
    
            elif StIdx == 'AqI':
                print('AqI: 2 phases Aqueous+Ice')
                mesh.loc[indexes,'P']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_aqu']=mesh.loc[indexes,'IC_var2']
                mesh.loc[indexes,'S_gas']=0
                mesh.loc[indexes,'S_hyd']=0
                mesh.loc[indexes,'S_icd']=1-mesh.loc[indexes,'IC_var2']
    
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var3']
                    
                    
            elif StIdx == 'IcH':
                print('IcH: 2 phases Ice+Hydrate')
                mesh.loc[indexes,'P']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_aqu']=0
                mesh.loc[indexes,'S_gas']=0
                mesh.loc[indexes,'S_hyd']=1-mesh.loc[indexes,'IC_var2']
                mesh.loc[indexes,'S_icd']=mesh.loc[indexes,'IC_var2']
    
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var3']
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var4']
                else:
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var3']    
                    
                    
            elif StIdx == 'AGH':
                print('AGH: 3 phases Aqueous+Hydrate+Gas')
                mesh.loc[indexes,'S_aqu']=mesh.loc[indexes,'IC_var2']
                mesh.loc[indexes,'S_gas']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_hyd']=1-(mesh.loc[indexes,'IC_var2']+mesh.loc[indexes,'IC_var1'])
                mesh.loc[indexes,'S_icd']=0
    
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var3']
                                        
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var4']
                    mesh.loc[indexes,'P']=v_Equil_P(mesh.loc[indexes,'IC_var4'],C=True)
                    
                else:
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var3']    
                    mesh.loc[indexes,'P']=v_Equil_P(mesh.loc[indexes,'IC_var3'],C=True)
    
    
            elif StIdx == 'AIG':
                print('AIG: 3 phases Aqueous+Ice+Gas')
                mesh.loc[indexes,'P']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_aqu']=mesh.loc[indexes,'IC_var2']
                mesh.loc[indexes,'S_hyd']=0
                
    
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var3']
                    mesh.loc[indexes,'S_gas']=mesh.loc[indexes,'IC_var4']
                    mesh.loc[indexes,'S_icd']=1-(mesh.loc[indexes,'IC_var2']+mesh.loc[indexes,'IC_var4'])
                    
                else:
                    mesh.loc[indexes,'S_gas']=mesh.loc[indexes,'IC_var3']
                    mesh.loc[indexes,'S_icd']=1-(mesh.loc[indexes,'IC_var2']+mesh.loc[indexes,'IC_var3'])
    
            elif StIdx == 'AIH':
                print('AIH: 3 phases Aqueous+Ice+Hydrate')
                mesh.loc[indexes,'P']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_aqu']=mesh.loc[indexes,'IC_var2']
                mesh.loc[indexes,'S_gas']=0
    
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var3']
                    mesh.loc[indexes,'S_icd']=mesh.loc[indexes,'IC_var4']
                    mesh.loc[indexes,'S_hyd']=1-(mesh.loc[indexes,'IC_var2']+mesh.loc[indexes,'IC_var4'])
                    
                else:
                    mesh.loc[indexes,'S_icd']=mesh.loc[indexes,'IC_var3']
                    mesh.loc[indexes,'S_hyd']=1-(mesh.loc[indexes,'IC_var2']+mesh.loc[indexes,'IC_var3'])
    
            elif StIdx == 'IGH':
                print('IGH: 3 phases Ice+Gas+Hydrate')
                mesh.loc[indexes,'S_gas']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_icd']=mesh.loc[indexes,'IC_var2']
                mesh.loc[indexes,'S_aqu']=0
                mesh.loc[indexes,'S_hyd']=1-(mesh.loc[indexes,'IC_var2']+mesh.loc[indexes,'IC_var1'])
    
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var3']
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var4']
                    mesh.loc[indexes,'P']=v_Equil_P(mesh.loc[indexes,'IC_var4'],C=True)
                    
                    
                else:
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var3']
                    mesh.loc[indexes,'P']=v_Equil_P(mesh.loc[indexes,'IC_var3'],C=True)
                    
            elif StIdx == 'QuP':
                print('QuP: 4 phases Ice+Gas+Hydrate+Aqueous')
                mesh.loc[indexes,'S_gas']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_aqu']=mesh.loc[indexes,'IC_var2']
    
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var3']
                    mesh.loc[indexes,'S_icd']=mesh.loc[indexes,'IC_var4']
                    mesh.loc[indexes,'S_hyd']=1-(mesh.loc[indexes,'IC_var4']+mesh.loc[indexes,'IC_var2']+mesh.loc[indexes,'IC_var1'])
                    
                else:
                    mesh.loc[indexes,'S_icd']=mesh.loc[indexes,'IC_var3']
                    mesh.loc[indexes,'S_hyd']=1-(mesh.loc[indexes,'IC_var3']+mesh.loc[indexes,'IC_var2']+mesh.loc[indexes,'IC_var1'])
                
        else:
            if StIdx == 'Gas':
                print('Gas: 1 phase Gas')
                mesh.loc[indexes,'P']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_aqu']=0
                mesh.loc[indexes,'S_gas']=1-mesh.loc[indexes,'IC_var3']
                mesh.loc[indexes,'S_hyd']=mesh.loc[indexes,'IC_var3']
                mesh.loc[indexes,'S_icd']=0
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var4']
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var5']
                else:
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var4']
            
            elif StIdx == 'Aqu':
                print('Aqu: 1 phase Aqueous')
                mesh.loc[indexes,'P']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_aqu']=1-mesh.loc[indexes,'IC_var3']
                mesh.loc[indexes,'S_gas']=0
                mesh.loc[indexes,'S_hyd']=mesh.loc[indexes,'IC_var3']
                mesh.loc[indexes,'S_icd']=0
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var4']
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var5']
                else:
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var4']
                    
                    
            elif StIdx == 'AqG':
                print('AqG: 2 phases Aqueous+Gas')
                mesh.loc[indexes,'P']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_aqu']=mesh.loc[indexes,'IC_var2']
                mesh.loc[indexes,'S_gas']=1-(mesh.loc[indexes,'IC_var2']+mesh.loc[indexes,'IC_var3'])
                mesh.loc[indexes,'S_hyd']=mesh.loc[indexes,'IC_var3']
                mesh.loc[indexes,'S_icd']=0
    
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var4']
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var5']
                else:
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var4']
                    
            elif StIdx == 'IcG':
                print('IcG: 2 phases Ice+Gas')
                mesh.loc[indexes,'P']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_aqu']=0
                mesh.loc[indexes,'S_gas']=1-mesh.loc[indexes,'IC_var2']-mesh.loc[indexes,'IC_var3']
                mesh.loc[indexes,'S_hyd']=mesh.loc[indexes,'IC_var3']
                mesh.loc[indexes,'S_icd']=mesh.loc[indexes,'IC_var2']
    
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var4']
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var5']
                else:
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var4']
                    
            elif StIdx == 'GsH':
                print('GsH: 2 phases Hydrate+Gas')
                mesh.loc[indexes,'P']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_aqu']=0
                mesh.loc[indexes,'S_gas']=mesh.loc[indexes,'IC_var2']
                mesh.loc[indexes,'S_hyd']=1-mesh.loc[indexes,'IC_var2']-mesh.loc[indexes,'IC_var3']
                mesh.loc[indexes,'S_icd']=mesh.loc[indexes,'IC_var3']
    
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var4']
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var5']
                else:
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var4']
        
    
            elif StIdx == 'AqH':
                print('AqH: 2 phases Aqueous+Hydrate')
                mesh.loc[indexes,'P']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_aqu']=mesh.loc[indexes,'IC_var2']
                mesh.loc[indexes,'S_gas']=0
                mesh.loc[indexes,'S_hyd']=1-mesh.loc[indexes,'IC_var2']
                mesh.loc[indexes,'S_icd']=0
    
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var4']
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var5']
                else:
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var4']
    
    
            elif StIdx == 'AqI':
                print('AqI: 2 phases Aqueous+Ice')
                mesh.loc[indexes,'P']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_aqu']=mesh.loc[indexes,'IC_var2']
                mesh.loc[indexes,'S_gas']=0
                mesh.loc[indexes,'S_hyd']=0
                mesh.loc[indexes,'S_icd']=1-mesh.loc[indexes,'IC_var2']
    
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var4']
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var5']
                else:
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var4']
                    
                    
            elif StIdx == 'IcH':
                print('IcH: 2 phases Ice+Hydrate')
                mesh.loc[indexes,'P']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_aqu']=0
                mesh.loc[indexes,'S_gas']=mesh.loc[indexes,'IC_var3']
                mesh.loc[indexes,'S_hyd']=1-mesh.loc[indexes,'IC_var2']-mesh.loc[indexes,'IC_var3']
                mesh.loc[indexes,'S_icd']=mesh.loc[indexes,'IC_var2']
    
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var4']
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var5']
                else:
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var4']    
                    
                    
            elif StIdx == 'AGH':
                print('AGH: 3 phases Aqueous+Hydrate+Gas')
                mesh.loc[indexes,'P']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_aqu']=mesh.loc[indexes,'IC_var2']
                mesh.loc[indexes,'S_gas']=mesh.loc[indexes,'IC_var3']
                mesh.loc[indexes,'S_hyd']=1-(mesh.loc[indexes,'IC_var2']+mesh.loc[indexes,'IC_var3'])
                mesh.loc[indexes,'S_icd']=0
    
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var4']
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var5']
                else:
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var4']    
    
    
            elif StIdx == 'AIG':
                print('AIG: 3 phases Aqueous+Ice+Gas')
                mesh.loc[indexes,'P']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_aqu']=mesh.loc[indexes,'IC_var2']
                mesh.loc[indexes,'S_hyd']=mesh.loc[indexes,'IC_var3']
    
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var4']
                    mesh.loc[indexes,'S_gas']=mesh.loc[indexes,'IC_var5']
                    mesh.loc[indexes,'S_icd']=1-(mesh.loc[indexes,'IC_var2']+mesh.loc[indexes,'IC_var3']+mesh.loc[indexes,'IC_var5'])
                    
                else:
                    mesh.loc[indexes,'S_gas']=mesh.loc[indexes,'IC_var4']
                    mesh.loc[indexes,'S_icd']=1-(mesh.loc[indexes,'IC_var2']+mesh.loc[indexes,'IC_var3']+mesh.loc[indexes,'IC_var4'])
    
            elif StIdx == 'AIH':
                print('AIH: 3 phases Aqueous+Ice+Hydrate')
                mesh.loc[indexes,'P']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_aqu']=mesh.loc[indexes,'IC_var2']
                mesh.loc[indexes,'S_icd']=mesh.loc[indexes,'IC_var3']
                mesh.loc[indexes,'S_hyd']=1-(mesh.loc[indexes,'IC_var2']+mesh.loc[indexes,'IC_var3'])
                mesh.loc[indexes,'S_gas']=0
    
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var4']
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var5']
                    
                else:
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var4']
    
            elif StIdx == 'IGH':
                print('IGH: 3 phases Ice+Gas+Hydrate')
                mesh.loc[indexes,'P']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_gas']=mesh.loc[indexes,'IC_var2']
                mesh.loc[indexes,'S_icd']=mesh.loc[indexes,'IC_var3']
                mesh.loc[indexes,'S_hyd']=1-(mesh.loc[indexes,'IC_var2']+mesh.loc[indexes,'IC_var3'])
                mesh.loc[indexes,'S_aqu']=0
    
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var4']
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var5']
                    
                else:
                    mesh.loc[indexes,'T']=mesh.loc[indexes,'IC_var4']
                    
            elif StIdx == 'QuP':
                print('QuP: 4 phases Ice+Gas+Hydrate+Aqueous')
                mesh.loc[indexes,'P']=mesh.loc[indexes,'IC_var1']
                mesh.loc[indexes,'S_aqu']=mesh.loc[indexes,'IC_var2']
                mesh.loc[indexes,'S_gas']=mesh.loc[indexes,'IC_var3']
    
                
                if Inh:
                    mesh.loc[indexes,'X_inh']=mesh.loc[indexes,'IC_var4']
                    mesh.loc[indexes,'S_icd']=mesh.loc[indexes,'IC_var5']
                    mesh.loc[indexes,'S_hyd']=1-(mesh.loc[indexes,'IC_var5']+mesh.loc[indexes,'IC_var2']+mesh.loc[indexes,'IC_var1'])
                    
                else:
                    mesh.loc[indexes,'S_icd']=mesh.loc[indexes,'IC_var4']
                    mesh.loc[indexes,'S_hyd']=1-(mesh.loc[indexes,'IC_var4']+mesh.loc[indexes,'IC_var2']+mesh.loc[indexes,'IC_var1'])
            
    return mesh

