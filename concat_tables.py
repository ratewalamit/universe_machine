import numpy as np
from tqdm import tqdm
import pandas as pd
from scipy.spatial import cKDTree as tree
from astropy.table import Table,vstack
import time
import os
import glob
from mpi4py import MPI
    
names='ID DescID UPID Flags Uparent_Dist X Y Z VX VY VZ M V MP VMP R Rank1 Rank2 RA RARank SM ICL SFR obs_SM obs_SFR SSFR SM/HM obs_UV'.split()  #column names
names+=["Parent_"+ name for name in names]
fp=pd.read_csv("/mnt/home/student/camit/Project_universemachine/DataStore/sfr_catalog_0.055623.txt",header=None,comment="#",sep=" ",names=np.array(names)[:int(len(names)/2)])
data_type=np.concatenate((fp.dtypes.values,fp.dtypes.values))
Lbox=250

stacked_table=Table(names=names,dtype=data_type)
for gal_type in "orphan","satellite":
    flist=np.sort(glob.glob(f"./temp/{gal_type}*"))
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
stacked_table.write(f"/mnt/home/student/camit/Project_universemachine/DataStore/combined_satellites_and_orphan_catalog_mimicing_observations_sfr_catalog_1.002310.fits")
    

