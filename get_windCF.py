import cdms2 as cdms,MV2 as MV, numpy as np,cdutil
import os,sys

###################### control variables
atm_year=8760
case_name='MERRA300.prod.assim.tavg1_2d_slv_Nx.2014'
######################

data_path="/lustre/scratch/leiduan/data_needed/"
fd_cam=os.listdir(data_path)

u_ci = 3 # cut in speed in m/s
u_r = 12 # rated speed in m/s
u_co = 25 # cut out speed in m/s

#fmask = cdms.open()

month_days={1:31,2:28,3:31,4:30,5:31,6:30,7:31,8:31,9:30,10:31,11:30,12:31}
month_days_array = np.array([31,28,31,30,31,30,31,31,30,31,30,31])

wcf = MV.array(np.zeros([atm_year,361,540])); wcf.id='wcf'; wcf.units='1'; wcf.missing_value = 1e20

count_num=1
for file in fd_cam:
        if file[:-7] == case_name and file[-3:]=='.nc':
                print (file, count_num)
                f=cdms.open(data_path+file)
                month = int(file[-7:-5])
                days = int(file[-5:-3])
                if month == 1:
                    position = (0 + days-1) * 24
                else:
                    position = (np.sum(month_days_array[:month-1]) + (days-1)) *24
                print (month, days, position)

                u50m_tmp = f('U50M')  #(24,361,540)
                v50m_tmp = f('V50M')  #(24,361,540)
                u10m_tmp = f('U10M')  #(24,361,540)
                v10m_tmp = f('V10M')  #(24,361,540)
                ws10m = np.hypot(u10m_tmp,v10m_tmp)
                ws50m = np.hypot(u50m_tmp,v50m_tmp)
                wsc = (np.log(ws50m)-np.log(ws10m))/(np.log(50.)-np.log(10.))
                ws100m_tmp = ws10m * (100./10.)**wsc   #(24,361,540)
                ws100m = MV.filled(ws100m_tmp,0.)               
 
                for hr_idx in range(24):
                    wcf[position+hr_idx,  ws100m[hr_idx] <  u_ci ] = 0.
                    wcf[position+hr_idx, (ws100m[hr_idx] >= u_ci) & (ws100m[hr_idx] <  u_r)  ] = ws100m[hr_idx, (ws100m[hr_idx] >= u_ci) & (ws100m[hr_idx] < u_r)]**3 / (u_r**3)
                    wcf[position+hr_idx, (ws100m[hr_idx] >= u_r)  & (ws100m[hr_idx] <= u_co) ] = 1.
                f.close
                count_num = count_num+1

# use NetCDF3 Classic format
cdms.setNetcdfShuffleFlag(0)      # netcdf3 classic...
cdms.setNetcdfDeflateFlag(0)      # netcdf3 classic...
cdms.setNetcdfDeflateLevelFlag(0) # netcdf3 classic...

fout=cdms.open(case_name+'_wcf100m12_duanlei.nc','w')
fout.write(wcf)
fout.close()

wcf_annual = cdutil.averager(wcf,axis=0,weights='equal')
wcf_annual.id='wcf_annual'
gout=cdms.open(case_name+'_wcf100m12_annual_duanlei.nc','w')
gout.write(wcf_annual)
gout.close()

