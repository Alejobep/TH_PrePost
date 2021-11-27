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
def read_Time_Series(path,t='days',P='Pa'):


    if os.path.isdir(path):
        folder = path
    else:
        folder = os.path.dirname(path)

    os.chdir(folder)

    for f in os.listdir():
        if r'Time_Series' in f:
            with open(f) as ts_file:
                first_line = ts_file.readline()

            ts_df = pd.read_fwf(f,skiprows = 1)

            if t == 'h':
                ts_df['Time [days]']*=24

                ts_df = ts_df.rename(columns = {'Time [days]':'Time [h]'})

            if P == 'bar':
                ts_df['P[Pa]']/=1e5 
                
                ts_df = ts_df.rename(columns = {'P[Pa]':'P[bar]'})
                

            name = f[:-12]

            if 'Subdomain' in first_line:
                try:
                    subdoms[name] = ts_df

                except:
                    subdoms=dict()
                    subdoms[name] = ts_df
                

            if 'Interface' in first_line:
                try:
                    interfs[name] = ts_df

                except:
                    interfs=dict()
                    interfs[name] = ts_df

    return subdoms, interfs
            



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




def read_out_file(path):
    import re

    if os.path.isdir(path):
        case = path
    else:
        case = os.path.dirname(path)

    ip_data = read_TH_data(path)
    diff = ip_data.MEMORY.processed['M_3']['binary_diffusion']



    for file in os.listdir(case):
        if file.endswith('.out'):
            op_file = file

    path = os.path.join(case, op_file)


    #Retrieve indexes where tables are located in Output file
    idxs=[]

    txt_1 = r'TOTAL TIME'
    txt_2 = r'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
    txt_3 = r'MASS FLOW RATES'

    sep_str = [txt_1, txt_2, txt_3]

    tb1_bool = False
    test = []

    ts = []

    with open(path, 'r', encoding='latin') as f:

        for idx,line in enumerate(f):

            if any(text in line for text in sep_str):
                counter=0
                idxs.append(idx)

    idxs = np.array(idxs)

    #Sources/sinks flag
    ss_flag = len(''.join(ip_data.GENER.raw).strip())>len('GENER')

    if diff:
        if ss_flag:
            idxs = idxs.reshape((len(idxs)//6),6)
        else:
            idxs = idxs.reshape((len(idxs)//5),5)
    else:
        if ss_flag:
            idxs = idxs.reshape((len(idxs)//4),4)
        else:
            idxs = idxs.reshape((len(idxs)//4),4)




    #Read time steps and headers
    with open(path, 'r', encoding='latin') as f:

        for idx,line in enumerate(f):
            if idx in idxs[:,0]+1:
                ts.append(float(line[:20]))

            if idx == idxs[0,0]+6:
                line = line.replace('S_Hydrate', 'S_hyd')
                line = line.replace('S_aqueous', 'S_aqu')
                line = line.replace('S_Ice', 'S_icd')
                line = line.replace('Pressure', 'P')
                line = line.replace('Temperature', 'T')
                line = line.replace('X_Inhibitor', 'X_inh')
                tb1_header = line.split()

            if idx == idxs[0,1]+8:
                line = line.replace('Phi', 'porosity')
                line = line.replace('PC_w', 'P_cap')
                tb2_header = line.split()

            if idx == idxs[0,2]+8:
                tb3_header = line

            if diff and idx == idxs[0,3]+2:
                tb4_header = line

    ts = np.array(ts)


    tb3_header = tb3_header.strip()


    tb3_header=re.sub('[\s]{2,}','*', tb3_header)
    tb3_header=re.sub('Flo CH4','Flo*CH4', tb3_header)
    tb3_header=re.sub('[\s]','_', tb3_header)
    tb3_header=tb3_header.split('*')

    if diff:
        tb4_header = tb4_header.strip()
        tb4_header = tb4_header.split()
        tb4_header = tb4_header[:2]+tb4_header[3::2]
        tb4_header[2]+='_CH4'
        tb4_header[3]+='_H20'
        tb4_header[4]+='_NaCl'

    #Get clean indexes that mark first and last lines of tables   
    clean_idxs = np.zeros((idxs.shape[0],idxs.shape[1]+idxs.shape[1]-2))
    clean_idxs[:,0] = idxs[:,0]+10
    clean_idxs[:,1] = idxs[:,1]-1

    clean_idxs[:,2] = idxs[:,1]+12
    clean_idxs[:,3] = idxs[:,2]-1

    clean_idxs[:,4] = idxs[:,2]+12
    clean_idxs[:,5] = idxs[:,3]-1

    if diff:
        clean_idxs[:,5] = idxs[:,3]-4
        clean_idxs[:,6] = idxs[:,3]+4
        clean_idxs[:,7] = idxs[:,4]-2

    clean_idxs = clean_idxs.astype(int)   


    #Extract tables
    def tb_logic(index):
        for col in range(clean_idxs.shape[1]//2):
            tb_idx = clean_idxs[:,col*2:col*2+2]

            if any(index>=idx[0] and index<idx[1] for idx in tb_idx):
                return False

        return True

    tb = pd.read_csv(path, skiprows= lambda x: tb_logic(x), names = ['raw'])

    tb_sizes = np.diff(clean_idxs, axis=1)[0,::2]
    chunk_size = tb_sizes.sum()

    idx_list = tb.index.values
    eval_list = idx_list-chunk_size*(idx_list//chunk_size)
    tb1 = tb.iloc[eval_list<tb_sizes.cumsum()[0]].copy()
    tb2 = tb.iloc[(eval_list>=tb_sizes.cumsum()[0]) & (eval_list<tb_sizes.cumsum()[1])].copy()
    tb3 = tb.iloc[(eval_list>=tb_sizes.cumsum()[1]) & (eval_list<tb_sizes.cumsum()[2])].copy()
    if diff:
        tb4 = tb.iloc[(eval_list>=tb_sizes.cumsum()[2]) & (eval_list<tb_sizes.cumsum()[3])].copy()


    #Parse tables    
    tb1 = tb1['raw'].str.split(expand=True)
    tb1 = tb1.rename(columns=dict(zip(tb1.columns.values, tb1_header)))

    for j, col in tb1.iteritems():
        query = (col.str.contains("[0-9][-][0-9]|[0-9][+][0-9]"))
        col.loc[query] = col.loc[query].str[:-4]+'E'+col.loc[query].str[-4:]


    tb1 = tb1.astype(dict(zip(tb1_header, ['str', 'int32']+['float64']*(len(tb1_header)-2))))

    tb2 = tb2['raw'].str.split(expand=True)
    tb2 = tb2.rename(columns=dict(zip(tb2.columns.values, tb2_header)))
    for j, col in tb2.iteritems():
        query = (col.str.contains("[0-9][-][0-9]|[0-9][+][0-9]"))
        col.loc[query] = col.loc[query].str[:-4]+'E'+col.loc[query].str[-4:]

    tb2 = tb2.astype(dict(zip(tb2_header, ['str', 'int32']+['float64']*(len(tb2_header)-2))))

    tb3 = tb3['raw'].str.split(expand=True)
    tb3 = tb3.rename(columns=dict(zip(tb3.columns.values, tb3_header)))
    for j, col in tb3.iteritems():
        query = (col.str.contains("[0-9][-][0-9]|[0-9][+][0-9]"))
        col.loc[query] = col.loc[query].str[:-4]+'E'+col.loc[query].str[-4:]

    tb3 = tb3.astype(dict(zip(tb3_header, ['str', 'str']+['float64']*(len(tb3_header)-2))))

    if diff:
        tb4 = tb4['raw'].str.split(expand=True)
        tb4 = tb4.rename(columns=dict(zip(tb4.columns.values, tb4_header)))
        for j, col in tb4.iteritems():
            query = (col.str.contains('\-|\+')) & (~col.str.contains('E'))
            col.loc[query] = col.loc[query].str[:-4]+'E'+col.loc[query].str[-4:]
        tb4 = tb4.astype(dict(zip(tb4_header, ['str', 'str']+['float64']*(len(tb4_header)-2))))


    #Drop repeated columns
    tb1 = tb1.drop(columns='INDEX')
    tb1 = tb1.reset_index(drop=True)

    tb2 = tb2.drop(columns=['ELEM', 'INDEX'])
    tb2 = tb2.reset_index(drop=True)

    tb3 = tb3.reset_index(drop=True)

    if diff:
        tb4 = tb4.drop(columns=['ELEM1', 'ELEM2'])
        tb4 = tb4.reset_index(drop=True)


    #Produce final tables
    data_grid = pd.concat([tb1,tb2], axis=1)

    if diff:
        data_conne = pd.concat([tb3,tb4], axis=1)
    else:
        data_conne = tb3

    #Clean ELEM column    
    data_grid['ELEM'] = data_grid['ELEM'].str.strip().str.strip('*')

    #Retrieve grid dimensions
    eleme = ip_data.ELEME.processed
    eleme = eleme.set_index('ElName')

    Tdim = len(eleme)
    Cdim = len(ip_data.CONNE.processed)

    time_grid = np.repeat(ts, Tdim)
    time_conne = np.repeat(ts, Cdim)

    MA12 = eleme.loc[data_grid['ELEM'], 'MA12'].values


    data_grid['time'] = time_grid
    data_conne['time'] = time_conne
    data_grid['MA12'] = MA12

    data_grid['x'] = eleme.loc[data_grid.ELEM, 'X'].values
    data_grid['y'] = eleme.loc[data_grid.ELEM, 'Y'].values
    data_grid['z'] = eleme.loc[data_grid.ELEM, 'Z'].values

    data_grid = data_grid.sort_values(by=['time','x','y','z'], ascending=[True,True,True,False])

    data_grid = data_grid.rename(columns={'ELEM': 'ElName'})

    data_grid = data_grid.set_index(['time', 'ElName'])

    data_conne = data_conne.rename_axis('C_IDX')



    return data_grid, data_conne


def tile_out_file(files_df):

    for idx,row in files_df.iterrows():

        try:
            iter_path = os.path.join(row.Name, fname)
            tdata_gr, tdata_conne = read_out_file(iter_path)

            min_t = tdata_gr.index.get_level_values(0).min()

            files_df.loc[idx,'t_min'] = min_t
            try:
                merged_gr = pd.concat([merged_gr,tdata_gr], sort = False)
                merged_conne = pd.concat([merged_conne,tdata_conne], sort = False)

            except:
                merged_gr = pd.concat([tdata_gr], sort = False)
                merged_conne = pd.concat([tdata_conne], sort = False)

            print(idx, row.Name)
            last_file = row.Name


        except:
            print('Last file processed:',last_file)
            continue


    return merged_gr, merged_conne



def get_output(*file):
    """
    Reads initialization parameters
    Process Plot_Data_Elem
    Returns a dataframe with both initialization and input
    """

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

    full_tdata.to_pickle(os.path.join(pklpath,'GRID_data.pkl'))


    return full_tdata



def get_output_full(*file):
    """
    Reads initialization parameters
    Process .out file
    Returns a dataframe with both initialization and input
    """

    if len(file)==1:
        print("Processing single file")

        data = read_TH_data(file[0])
        mesh = process_init(file[0])
        tdata_gr, tdata_conne = read_out_file(file[0])

        pklpath = os.path.dirname(file[0])

    else:
        print("Processing multiple files")        
        n_files = len(file)
        sorted_files = pd.DataFrame(index = range(n_files), columns=('Name','t_min','ctime'))
        
        sorted_files['Name']=file
        sorted_files['ctime']=sorted_files.Name.map(os.path.getctime)

        sorted_files = sorted_files.sort_values(by='ctime')
    
        
        for idx,row in sorted_files.iterrows():
            
            try:
                tdata_gr, tdata_conne = read_out_file(row.Name)

                min_t = tdata_gr.index.get_level_values(0).min()
                
                sorted_files.loc[idx,'t_min'] = min_t
                try:
                    merged_gr = pd.concat([merged_gr,tdata_gr], sort = False)
                    merged_conne = pd.concat([merged_conne,tdata_conne], sort = False)
                
                except:
                    merged_gr = pd.concat([tdata_gr], sort = False)
                    merged_conne = pd.concat([tdata_conne], sort = False)

                # print(idx, row.Name)
                last_file = row.Name

           
            except:
                print('Last file processed:',last_file)
                continue

        
        tdata_gr = merged_gr
        tdata_conne = merged_conne
                
        sorted_files = sorted_files.sort_values(by='t_min')
        
        sorted_files = sorted_files.reset_index()   

        # for f_idx, f_row in sorted_files.iterrows():
        #     print('iteration No. {:d}. {:s}'.format(f_idx,f_row['Name']))


        file_0 = sorted_files.loc[0,'Name']
        mesh = process_init(file_0)
        data = read_TH_data(file_0)

        pklpath = os.path.dirname(os.path.dirname(file[0]))



    tsteps = tdata_gr.index.get_level_values(0).drop_duplicates()
    t1 = tsteps[0]

    tdata_gr_0 = tdata_gr.loc[t1].copy()
    tdata_gr_0.iloc[:,3:] = np.nan
    tdata_gr_0[['x', 'y', 'z']] = tdata_gr.loc[t1, ['x', 'y', 'z']]


    

    #Set porosity and permeability from ROCKS
    rocks = data.ROCKS.processed
    rocks = rocks.set_index('Name')

    #Reindex INCON mesh to ELNAme
    mesh = mesh.set_index('ElName')

#     tdata_gr_0.porosity = rocks.loc[mesh.MA12,'Poros'].values
#     tdata_gr_0['perm_abs'] = rocks.loc[mesh.MA12,'Perm1'].values 

    tdata_gr_0.loc[:,['P','T','S_hyd','S_gas','S_aqu','S_icd','X_inh','MA12']]=mesh[['P','T','S_hyd','S_gas','S_aqu','S_icd','X_inh', 'MA12']]

    #Append 1st time step and create new multi-Index
    tsteps = np.append(0,tsteps)

    print(len(tsteps),'steps')

    MI_iterables = [tsteps,tdata_gr_0.index]
    df_MIndex = pd.MultiIndex.from_product(MI_iterables, names=['time', 'ElName'])

    #Concatenate
    
    
    full_tdata_gr = pd.concat([tdata_gr_0,tdata_gr])
    full_tdata_gr = full_tdata_gr.set_index(df_MIndex)

    full_tdata_gr['I'] = np.tile(mesh.I.values, len(tsteps))
    full_tdata_gr['J'] = np.tile(mesh.J.values, len(tsteps))
    full_tdata_gr['K'] = np.tile(mesh.K.values, len(tsteps))
    full_tdata_gr['dx'] = np.tile(mesh.dX.values, len(tsteps))
    full_tdata_gr['dy'] = np.tile(mesh.dY.values, len(tsteps))
    full_tdata_gr['dz'] = np.tile(mesh.dZ.values, len(tsteps))


    # full_tdata_gr['MA12'] = np.tile(mesh.MA12.values, len(tsteps))

    #Reorganize columns
    cols = full_tdata_gr.columns.to_list()
    
    cols = cols[-10:]+cols[:-10]
    cols = cols[1:10]+cols[:1]+cols[10:]
    full_tdata_gr = full_tdata_gr[cols]

#     #Change dtypes to float
#     float_cols = cols[:3]+cols[10:]
#     dtype_dict = dict(zip(float_cols,['float']*len(float_cols)))
#     full_tdata_gr = full_tdata_gr.astype(dtype_dict)

    #Convert P from Pa to bar
    full_tdata_gr['P'] /= 1e5
    full_tdata_gr['P_CH4'] /= 1e5
    full_tdata_gr['P_EqHydr'] /= 1e5
    full_tdata_gr['P_SatWat'] /= 1e5
    full_tdata_gr['P_cap'] /= 1e5

    full_tdata_gr = full_tdata_gr.rename(columns={"C-CH4inGas": "C_CH4inGas", "C-CH4inAqu": "C_CH4inAqu"})

    #Reset index to time only
    full_tdata_gr = full_tdata_gr.reset_index(level=[1])

    full_tdata_gr.to_pickle(os.path.join(pklpath,'GRID_data_full.pkl'))
    tdata_conne.to_pickle(os.path.join(pklpath,'CONNE_data.pkl'))

    return full_tdata_gr





   






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

