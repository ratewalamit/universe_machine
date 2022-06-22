from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import sys

fo=fits.open("/mnt/home/student/camit/Project_universemachine/DataStore/orphan_catalog_sfr_catalog_1.002310.fits")[1].data 
fs=fits.open("/mnt/home/student/camit/Project_universemachine/DataStore/satellites_catalog_sfr_catalog_1.002310.fits")[1].data 
Lbox=250.

ssfr_fo_idx=fo["SSFR"]<1e-11
ssfr_fs_idx=fs["SSFR"]<1e-11
fo=fo[ssfr_fo_idx]
fs=fs[ssfr_fs_idx]

ox,oy,oz=fo["x"],fo["y"],fo["z"]
opx,opy,opz=fo["Parent_x"],fo["Parent_y"],fo["Parent_z"]
sx,sy,sz=fs["x"],fs["y"],fs["z"] 
spx,spy,spz=fs["Parent_x"],fs["Parent_y"],fs["Parent_z"]

odx=np.absolute(opx-ox)
ody=np.absolute(opy-oy)
odz=np.absolute(opz-oz)
sdx=np.absolute(spx-sx)
sdy=np.absolute(spy-sy) 
sdz=np.absolute(spz-sz) 

odx[odx>Lbox/2]=Lbox-odx[odx>Lbox/2]
ody[ody>Lbox/2]=Lbox-ody[ody>Lbox/2]
odz[odz>Lbox/2]=Lbox-odz[odz>Lbox/2]  
sdx[sdx>Lbox/2]=Lbox-sdx[sdx>Lbox/2] 
sdy[sdy>Lbox/2]=Lbox-sdy[sdy>Lbox/2] 
sdz[sdz>Lbox/2]=Lbox-sdz[sdz>Lbox/2] 
odz_idx=np.abs(odz)<1.
sdz_idx=np.abs(sdz)<1.

oR=np.sqrt(odx**2+ody**2+odz**2)[odz_idx]
sR=np.sqrt(sdx**2+sdy**2+sdz**2)[sdz_idx]
oRp=np.sqrt(odx**2+ody**2)[odz_idx]
sRp=np.sqrt(sdx**2+sdy**2)[sdz_idx]

#plt.hist(oR[oR<5],histtype="step",bins=30,label="orphan")
#plt.hist(sR[sR<5],histtype="step",bins=30,label="satellite")
#plt.legend()
for i in((0.1,0.3),(0.3,0.5),(0.5,0.7),(0.7,0.9)):
    bmin,bmax=i[0],i[1]
    idxo=(oR>bmin) & (oR<bmax)
    idxs=(sR>bmin) & (sR<bmax)
    idxo_p=(oRp>bmin) & (oRp<bmax)
    idxs_p=(sRp>bmin) & (sRp<bmax)
   
    print(oR[idxo].size/sR[idxs].size,oRp[idxo_p].size/sRp[idxs_p].size)
    print(f"{oR[idxo].size}/{sR[idxs].size},{oRp[idxo_p].size}/{sRp[idxs_p].size}")
