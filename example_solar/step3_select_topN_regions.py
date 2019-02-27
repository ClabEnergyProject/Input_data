#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 14:01:15 2018

@author: leiduan
"""

import cdms2 as cdms,MV2 as MV,numpy as np
import cdutil, regionmask
#from func import cal_gridarea

region_mask_list = { 0:'US_mask'}

def set_axes(var):
    var_new = MV.array(var)
    var_new.setAxis(0,lat)
    var_new.setAxis(1,lon)
    return var_new

def select_CF(sCF, wCF, idx, *i):
    if idx == 1.:
        s_avg_tmp = cdutil.averager(sCF,axis='yx')
        w_avg_tmp = cdutil.averager(wCF,axis='yx')
        
    elif idx == 2.:
        # The numbers of 0.15 for solar and wind come from Sterl et al., 2018
        # The numders of 0.01 for solar 0.10 for wind come from Jacobson et al., 2017
        s_thresholds = 0.10
        w_thresholds = 0.10
        sCF[sCF<s_thresholds] = 0.
        wCF[wCF<w_thresholds] = 0.
        masked_sCF = set_axes(MV.masked_equal(sCF,0.))
        masked_wCF = set_axes(MV.masked_equal(wCF,0.))
        s_avg_tmp = cdutil.averager(masked_sCF,axis='yx')
        w_avg_tmp = cdutil.averager(masked_wCF,axis='yx')
        
    elif idx ==3.:
        # thresholds are selected to produce similar US annual mean values with EIA data
        p_threshold_s = (int(i[0])+1) * 0.1
        p_threshold_w = (int(i[0])+1) * 0.1
        print p_threshold_s, p_threshold_w
        
        s_reshape = np.sort(np.array(MV.filled(sCF,0)).ravel())[::-1]
        w_reshape = np.sort(np.array(MV.filled(wCF,0)).ravel())[::-1]
        s_no_zero = s_reshape[s_reshape!=0.]
        w_no_zero = w_reshape[w_reshape!=0.]
        num_s = int(len(s_no_zero)*p_threshold_s)
        num_w = int(len(w_no_zero)*p_threshold_w)
        s_thresholds = s_no_zero[num_s-1]
        w_thresholds = w_no_zero[num_w-1]
        sCF[sCF<s_thresholds] = 0.
        wCF[wCF<w_thresholds] = 0.
        masked_sCF = set_axes(MV.masked_equal(sCF,0.))
        masked_wCF = set_axes(MV.masked_equal(wCF,0.))
        s_avg_tmp = cdutil.averager(masked_sCF,axis='yx')
        w_avg_tmp = cdutil.averager(masked_wCF,axis='yx')
        
        masked_sCF2 = MV.array(masked_sCF/masked_sCF)        
        masked_sCF2.id = 'us_mask_sCF_' + str(p_threshold_s)
        masked_wCF2 = MV.array(masked_wCF/masked_wCF)
        masked_wCF2.id = 'us_mask_wCF_' + str(p_threshold_w)
        
        g.write(masked_sCF2)
        g.write(masked_wCF2)
        
    return s_avg_tmp, w_avg_tmp

p_threshold_list = np.arange(100)*0.01+0.01

def test_thresholds(sCF, wCF):
    s_avg_tmp = np.zeros([len(p_threshold_list)])
    w_avg_tmp = np.zeros([len(p_threshold_list)])
    for i in range(len(p_threshold_list)):
        tmp1 = np.copy(sCF)
        tmp2 = np.copy(wCF)  
        p_threshold = p_threshold_list[i]
        s_reshape = np.sort(tmp1.ravel())[::-1]
        w_reshape = np.sort(tmp2.ravel())[::-1]
        s_no_zero = s_reshape[s_reshape!=0.]
        w_no_zero = w_reshape[w_reshape!=0.]
        num_s = int(len(s_no_zero)*p_threshold)
        num_w = int(len(w_no_zero)*p_threshold)
        s_thresholds = s_no_zero[num_s-1]
        w_thresholds = w_no_zero[num_w-1]
        
        tmp1[tmp1<s_thresholds] = 0.
        tmp2[tmp2<w_thresholds] = 0.
        masked_sCF = set_axes(MV.masked_equal(tmp1,0.))
        masked_wCF = set_axes(MV.masked_equal(tmp2,0.))
        s_avg_tmp[i] = cdutil.averager(masked_sCF,axis='yx')
        w_avg_tmp[i] = cdutil.averager(masked_wCF,axis='yx')
    return s_avg_tmp, w_avg_tmp

######## ######## ######## ######## ######## ######## ######## ######## 
######## start calculation
######## ######## ######## ######## ######## ######## ######## ######## 

f_mask = cdms.open('SWGDN.nc')
v=f_mask('SWGDN'); lat=v.getAxis(1); lon=v.getAxis(2) #; gridarea = cal_gridarea(lat,lon)
f_mask.close()

f_land_mask = cdms.open('land_sea_mask_merra.nc4')
land_mask_tmp = f_land_mask('FRLAND',squeeze=1)
land_mask_tmp[land_mask_tmp>=0.5] = 1.
land_mask_tmp[land_mask_tmp<0.5]  = 0.
land_mask = MV.array(MV.masked_equal(land_mask_tmp,0.))
land_mask.setAxis(0,lat);      land_mask.setAxis(1,lon)
f_land_mask.close()

fsolar = cdms.open('MERRA2_400.tavg1_2d_rad_Nx.2015_scf_annual.nc') # need to change
s = MV.array(fsolar('scf_annual')); s.setAxis(0,lat); s.setAxis(1,lon)
fsolar.close()
fwind = cdms.open('MERRA2_400.tavg1_2d_slv_Nx.2015_wcf100m031225_annual.nc') # need to change
w = MV.array(fwind('wcf_annual')); w.setAxis(0,lat); w.setAxis(1,lon)
fwind.close()

f_us_mask = cdms.open('US_mask.nc')
us_mask = f_us_mask('US_mask') * 0. + 1. 

"""
s_us = np.array( MV.filled(s*us_mask*land_mask,0) )
w_us = np.array( MV.filled(w*us_mask*land_mask,0) )
s_US, w_US = test_thresholds(s_us, w_us)
from matplotlib import pyplot as plt
plt.plot(p_threshold_list,s_US,c='firebrick')
plt.plot(p_threshold_list,w_US,c='royalblue')
plt.show()
#plt.savefig('us_CF_change.ps')
plt.clf()
"""

s_us_avg=np.zeros(10)
w_us_avg=np.zeros(10)
g = cdms.open('selected_us_mask.nc','w')
for i in range(10):
    s_us = MV.array(s * us_mask * land_mask); s_us.setAxis(0,lat); s_us.setAxis(1,lon)
    w_us = MV.array(w * us_mask * land_mask); w_us.setAxis(0,lat); w_us.setAxis(1,lon)
    s_us_avg[i], w_us_avg[i] = select_CF(s_us, w_us, 3., i)
g.close()
