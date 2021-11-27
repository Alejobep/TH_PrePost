

def RunMeshMaker(file_path):
    import subprocess
    import os
    print("Run MESHMAKER")    
    mm_path= r'"/mnt/c/Users/gpb/OneDrive - Equinor/ET/GasHydrates/TOUGH/MeshMaker/MeshMaker/Executable_MeshMaker/MM_LNX.bin"'
    
    start_cwd = os.getcwd()

    file_in=os.path.basename(file_path)
    file_out=file_in[:-2]+r'out'

    cmd=mm_path+r' <'+file_in+r'> '+file_out

    wd=os.path.dirname(file_path)
    
    p=subprocess.Popen(cmd,shell=True,cwd=wd,stdout=subprocess.PIPE)


    while p.wait() is None:
        print("Still working...")

        
def gridSpecs(w,d,h,TBoun=False,BBoun=False,bound_h=1e-4):
    """
    w=Array containing the horizontal width in X direction of each cell
    d=Array containing the horizontal depth in Y direction of each cell
    h=Array containing the vertical height in Z direction of each cell
    TBoun= Top boundary
    BBoun= Bottom boundary
    bound_h=height of boundary in m
    """
    import numpy as np
    
    w=np.array(w)
    d=np.array(d)
    h=np.array(h)
    
      
    
    
    # w_cat=np.unique(w)
    # d_cat=np.unique(d)
 
    w_counts=[]
    w_cat=[]
    w_idx=[]
    it='na'
    
    for idx,w_val in enumerate(w):
        if w_val != it:
            w_cat.append(w_val)
            w_idx.append(idx)
            it=w_val
    
    w_counts=w_idx.copy()
    w_counts.append(len(w))
    w_counts=np.array(w_counts)
    w_counts=np.diff(w_counts)  




    d_counts=[]
    d_cat=[]
    d_idx=[]
    it='na'
    
    for idx,d_val in enumerate(d):
        if d_val != it:
            d_cat.append(d_val)
            d_idx.append(idx)
            it=d_val
    
    d_counts=d_idx.copy()
    d_counts.append(len(d))
    d_counts=np.array(d_counts)
    d_counts=np.diff(d_counts)    
    
    
    h_counts=[]
    h_cat=[]
    h_idx=[]
    it='na'
    
    for idx,h_val in enumerate(h):
        if h_val != it:
            h_cat.append(h_val)
            h_idx.append(idx)
            it=h_val
    
    h_counts=h_idx.copy()
    h_counts.append(len(h))
    h_counts=np.array(h_counts)
    h_counts=np.diff(h_counts)
    
    
    
    XYZ_2=''

    print('Grid dimensions')
    
    if len(w_cat)>10:
        XYZ_2+='NX   '+'{:5d}'.format(len(w))
        for i,INC in enumerate(w):
            if i%8==0:
                XYZ_2+='\n'
            XYZ_2+='{:10.2e}'.format(INC)
        XYZ_2+='\n'   
        print('{:d}'.format(len(w))+'\tcolumn(s) of variable width')
    else:
        for i,cat in enumerate(w_cat):
            XYZ_2+='NX   '+'{:5d}'.format(w_counts[i])+' '+'{:.3e}'.format(cat)+'\n'
        print('{:d}'.format(len(w))+'\tcolumn(s) of '+'{:.1f}'.format(w_cat[0])+' m in X direction')
        
    
   
    if len(d_cat)>5:
        XYZ_2+='\nNY   '+'{:5d}'.format(len(d))
        for i,INC in enumerate(d):
            if i%8==0:
                XYZ_2+='\n'
            XYZ_2+='{:10.2e}'.format(INC)
        
        print('{:d}'.format(len(d))+'\tcolumn(s) of variable depth')
    else:
        for i,cat in enumerate(d_cat):
            XYZ_2+='NY   '+'{:5d}'.format(d_counts[i])+' '+'{:.3e}'.format(cat)+'\n'
        print('{:d}'.format(len(d))+'\tcolumn(s) of '+'{:.1f}'.format(d_cat[0])+' m in Y direction')

    
    if len(h_cat)>5:
        XYZ_2+='\nNY   '+'{:5d}'.format(len(d))
        for i,INC in enumerate(d):
            if i%8==0:
                XYZ_2+='\n'
            XYZ_2+='{:10.2e}'.format(INC)
        
        print('{:d}'.format(len(d))+'\tcolumn(s) of variable depth')
    else:
        for i,cat in enumerate(h_cat):
            XYZ_2+='NZ   '+'{:5d}'.format(h_counts[i])+' '+'{:.3e}'.format(cat)+'\n'
        print('{:d}'.format(len(h))+'\tcolumn(s) of '+'{:.1f}'.format(h_cat[0])+' m in Z direction')
    
       
        
    return XYZ_2