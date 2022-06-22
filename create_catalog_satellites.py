import numpy as np
from tqdm import tqdm
import pandas as pd
from scipy.spatial import KDTree as tree
from astropy.table import Table
h=0.7
Mlimit=h*1e14   #hinvMsun    (in catalog given in Msun)
fp=pd.read_csv("/mnt/home/student/camit/Project_universemachine/DataStore/sfr_catalog_1.002310.txt",header=None,comment="#",sep=" ")
idx_satellites=(fp[0].values<1e15) & (fp[2]!=-1) #halos with haloid>1e15 are orphan
fp_idx=fp[idx_satellites]

parent_id=fp[0].values
satellites_parent_id=fp_idx[2].values
print("Making tree")
my_tree=tree(parent_id.reshape(-1,1))
print("Tree making done")
counter=0
names='ID DescID UPID Flags Uparent_Dist X Y Z VX VY VZ M V MP VMP R Rank1 Rank2 RA RARank SM ICL SFR obs_SM obs_SFR SSFR SM/HM obs_UV Parent_ID Parent_DescID Parent_UPID Parent_Flags Parent_Uparent_Dist Parent_X Parent_Y Parent_Z Parent_VX Parent_VY Parent_VZ Parent_M Parent_V Parent_MP Parent_VMP Parent_R Parent_Rank1 Parent_Rank2 Parent_RA Parent_RARank Parent_SM Parent_ICL Parent_SFR Parent_obs_SM Parent_obs_SFR Parent_SSFR Parent_SM/HM Parent_obs_UV'
column_names=names.split()
data_type=np.concatenate((fp.dtypes.values,fp.dtypes.values))
df_store=Table(names=column_names,dtype=data_type)
for i in tqdm(range(satellites_parent_id.size)):
    idx=my_tree.query_ball_point(satellites_parent_id[i].reshape(-1,1),0.0001)
    if len(idx[0])==0:
        continue
    satellites_info=fp_idx.iloc[i].values
    parent_info=fp.iloc[idx[0]].values[0]
    if parent_info[11]<Mlimit:
        continue
    if (satellites_info[2] != parent_info[0]):
        print("Satellite parent not found in original catalog") 
        break
    #print(satellites_info,parent_info)
    counter+=1
    df_store.add_row(np.concatenate((satellites_info,parent_info)))

idx_parent_not_central=df_store["Parent_UPID"]!=-1
df_save=df_store[~idx_parent_not_central]
print(f"# of satellites satellites= {counter}") 
df_save.write("/mnt/home/student/camit/Project_universemachine/DataStore/satellites_catalog_sfr_catalog_1.002310.fits",overwrite=True)
