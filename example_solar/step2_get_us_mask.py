from matplotlib import pyplot as plt
import regionmask, numpy as np, MV2 as MV
import cdms2 as cdms

f_axis = cdms.open('SWGDN.nc')
v=f_axis('SWGDN')
lat=v.getAxis(1)
lon=v.getAxis(2)
f_axis.close()

f_land_mask = cdms.open('land_sea_mask_merra.nc4')
land_mask_tmp = f_land_mask('FRLAND',squeeze=1)
land_mask_tmp[land_mask_tmp>=0.5] = 1.
land_mask_tmp[land_mask_tmp<0.5]  = 0.
land_mask = MV.array(MV.masked_equal(land_mask_tmp,0.))
land_mask.setAxis(0,lat);      land_mask.setAxis(1,lon)
f_land_mask.close()

countries = regionmask.defined_regions.natural_earth.countries_110
mask2 = countries.mask(lon, lat)
mask2[:,(lon[:]<=-125)] = 999.
mask_ma = MV.masked_not_equal(mask2,4)
print (mask_ma.shape)

mask_ma = mask_ma*0+1.

US_mask = MV.array(mask_ma * land_mask)
US_mask.id='US_mask'
US_mask.setAxis(0,lat)
US_mask.setAxis(1,lon)

g=cdms.open('US_mask.nc','w')
g.write(US_mask)
g.close()
