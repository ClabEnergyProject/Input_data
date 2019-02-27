import cdms2 as cdms,MV2 as MV, cdutil, numpy as np

f_mask = cdms.open('SWGDN.nc')
v=f_mask('SWGDN'); lat=v.getAxis(1); lon=v.getAxis(2)
f_mask.close()

region_mask_slist = { 0:'averaged_us_mask_sCF_0.1',   1:'averaged_us_mask_sCF_0.2',   2:'averaged_us_mask_sCF_0.3',   
                      3:'averaged_us_mask_sCF_0.4',   4:'averaged_us_mask_sCF_0.5',   5:'averaged_us_mask_sCF_0.6', 
                      6:'averaged_us_mask_sCF_0.7',   7:'averaged_us_mask_sCF_0.8',   8:'averaged_us_mask_sCF_0.9',
                      9:'averaged_us_mask_sCF_1.0'}

region_mask_wlist = { 0:'averaged_us_mask_wCF_0.1',   1:'averaged_us_mask_wCF_0.2',   2:'averaged_us_mask_wCF_0.3',   
                      3:'averaged_us_mask_wCF_0.4',   4:'averaged_us_mask_wCF_0.5',   5:'averaged_us_mask_wCF_0.6', 
                      6:'averaged_us_mask_wCF_0.7',   7:'averaged_us_mask_wCF_0.8',   8:'averaged_us_mask_wCF_0.9',
                      9:'averaged_us_mask_wCF_1.0'}

f_scf = cdms.open('us_averaged_scf_series2015.nc')
f_wcf = cdms.open('us_averaged_wcf_series2015.nc')

scf = np.zeros([10,8760])
wcf = np.zeros([10,8760])

for idx in range(10):
    scf[idx] = f_scf(region_mask_slist[idx], squeeze=1)
    wcf[idx] = f_wcf(region_mask_wlist[idx], squeeze=1)

scf_avg = cdutil.averager(scf,axis=1,weights='equal')
wcf_avg = cdutil.averager(wcf,axis=1,weights='equal')

from matplotlib import pyplot as plt
ind = np.arange(10)
width = 0.2
fig, ax = plt.subplots()
rects1 = ax.bar(ind - width/2, scf_avg, width*0.8, yerr=0., color='firebrick')
rects2 = ax.bar(ind + width/2, wcf_avg, width*0.8, yerr=0., color='royalblue')
plt.ylim(0,0.5)
plt.savefig('srex_region_annual.ps')
#plt.show()
plt.clf()

axis = np.arange(365)+1
for idx in range(10):
    scf_reshape = cdutil.averager(scf[idx].reshape(-1,24),axis=1,weights='equal')
    wcf_reshape = cdutil.averager(wcf[idx].reshape(-1,24),axis=1,weights='equal')
    plt.plot(axis, scf_reshape, c='firebrick')
    plt.plot(axis, wcf_reshape, c='royalblue')
    plt.ylim(0,1)
    plt.savefig('us_CF'+str( (idx+1)*0.1 )+'.ps')
    #plt.show()
    plt.clf()




