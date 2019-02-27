import cdms2 as cdms, MV2 as MV, cdutil
import numpy as np

region_mask_list = { 0:'us_mask_sCF_0.1',   1:'us_mask_sCF_0.2',   2:'us_mask_sCF_0.3',   3:'us_mask_sCF_0.4',   4:'us_mask_sCF_0.5',   5:'us_mask_sCF_0.6', 
                     6:'us_mask_sCF_0.7',   7:'us_mask_sCF_0.8',   8:'us_mask_sCF_0.9',   9:'us_mask_sCF_1.0'}

f_mask = cdms.open('SWGDN.nc')
v=f_mask('SWGDN')
lat=v.getAxis(1)
lon=v.getAxis(2)
f_mask.close()

fm = cdms.open('selected_us_mask.nc')
fs = cdms.open('')

scf = MV.array(fs('scf',squeeze=1))
scf[scf<0] = 0.
fs.close()

print ('I can open it here 1')

# use NetCDF3 Classic format
cdms.setNetcdfShuffleFlag(0)      # netcdf3 classic...
cdms.setNetcdfDeflateFlag(0)      # netcdf3 classic...
cdms.setNetcdfDeflateLevelFlag(0) # netcdf3 classic...

#"""
g=cdms.open('us_averaged_scf_series2015.nc','w')
len_axis = len(scf.getAxis(0))
for idx in range(10):
    print (idx)
    new_data = MV.array(np.zeros(len_axis))
    new_data.id = 'averaged_' + region_mask_list[idx]
    mask_idx = MV.array(fm(region_mask_list[idx]))
    for i in range(len_axis):
        scf_idx = scf[i] * mask_idx
        scf_idx.setAxis(0,lat)
        scf_idx.setAxis(1,lon)
        new_data[i] = cdutil.averager(scf_idx, axis='yx')
    g.write(new_data)
    del(scf_idx)
    del(mask_idx)
    del(new_data)
g.close()
#""" 

