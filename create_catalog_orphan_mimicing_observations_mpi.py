import numpy as np
from tqdm import tqdm
import pandas as pd
from scipy.spatial import cKDTree as tree
from astropy.table import Table,vstack
import time
import os
import glob
from mpi4py import MPI
comm=MPI.COMM_WORLD
rank=comm.Get_rank()
size=comm.Get_size()

h=0.7
Lbox=250.
Mlimit=h*1e14   #hinvMsun    (in catalog given in Msun)
names='ID DescID UPID Flags Uparent_Dist X Y Z VX VY VZ M V MP VMP R Rank1 Rank2 RA RARank SM ICL SFR obs_SM obs_SFR SSFR SM/HM obs_UV'.split()  #column names
names+=["Parent_"+ name for name in names]

fp=pd.read_csv("/mnt/home/student/camit/Project_universemachine/DataStore/sfr_catalog_1.002310.txt",header=None,comment="#",sep=" ",names=np.array(names)[:int(len(names)/2)])
idx_centrals= (fp["UPID"]==-1) & (fp["M"]>Mlimit)
idx_orphans=(fp['ID'].values>1e15) & (fp['UPID']!=-1) #halos with haloid>1e15 are orphan
fp_centrals=fp[idx_centrals]
fp_orphans=fp[idx_orphans]
print("Making tree")

sx,sy,sz=fp_orphans["X"].values,fp_orphans["Y"].values,fp_orphans["Z"].values
cx,cy,cz=fp_centrals["X"].values,fp_centrals["Y"].values,fp_centrals["Z"].values
orphan_coordinates=np.column_stack((sx,sy,sz))
cluster_coordinates=np.column_stack((cx,cy,cz))
orphan_tree=tree(np.column_stack((sx,sy)),boxsize=Lbox) 
print("Tree making done")
data_type=np.concatenate((fp.dtypes.values,fp.dtypes.values))
df_store=Table(names=names,dtype=data_type)
chunk_size=int(cx.size/(size-1))
if not os.path.isdir("./temp"):
    craete_temp_dir=os.system("mkdir -p ./temp")
if rank!=0:
    if rank==size-1:
        start,end=chunk_size*(rank-1),cx.size
    else:
        start,end=chunk_size*(rank-1),chunk_size*rank
    print(start,end)
    for i in np.arange(start,end):
        idx=orphan_tree.query_ball_point((cx[i],cy[i]),1)
        if len(idx)<1:
            continue
        zp=np.abs(sz[idx]-cz[i])
        pbc_zp=zp>Lbox/2.
        zp[pbc_zp]=Lbox-zp[pbc_zp]
        idx_zp=zp<50   
        orphans_info=fp_orphans.iloc[idx][idx_zp].values
        parent_info=fp_centrals.iloc[i].values
        if rank==size-1:
           #print(f"Done {(i-start)*100/(end-start)}%")
           print(len(orphans_info))
        for ii in orphans_info:
            df_store.add_row(np.concatenate((ii,parent_info)))
        df_store.write("./temp/%s.fits"%rank,overwrite=True)
        req=comm.isend(rank,dest=0)
        req.wait()
"""
if rank==0:
    stacked_table=Table(names=names,dtype=data_type)
    for i in range(1,size):
       rank_p=comm.irecv(source=i)
       rank_p.wait()
    time.sleep(2) 
    flist=np.sort(glob.glob("./temp/*.fits"))
    for files in flist:
        tmp_table=Table.read(files)
        print(files)
        stacked_table=vstack([stacked_table,tmp_table]) 
    Proj_R_x=np.abs(stacked_table["X"]-stacked_table["Parent_X"])
    Proj_R_y=np.abs(stacked_table["Y"]-stacked_table["Parent_Y"])
    pbc_Proj_x=Proj_R_x>Lbox/2.
    pbc_Proj_y=Proj_R_y>Lbox/2.
    Proj_R_x[pbc_Proj_x]=Lbox-Proj_R_x[pbc_Proj_x]
    Proj_R_y[pbc_Proj_y]=Lbox-Proj_R_y[pbc_Proj_y]
    stacked_table.add_column(np.sqrt( Proj_R_x**2 + Proj_R_y**2),name="Proj_R",index=0)
    stacked_table.write("/mnt/home/student/camit/Project_universemachine/DataStore/parallel_orphan_catalog_mimicing_observations_sfr_catalog_1.002310.fits",overwrite=True)
#        
"""
            

