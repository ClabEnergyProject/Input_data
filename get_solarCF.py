import cdms2 as cdms,MV2 as MV, numpy as np,cdutil, os

###################### control variables
hour_in_years=8760
case_name='MERRA300.prod.assim.tavg1_2d_rad_Nx.2014'
max_clearness_index=1.0
tilt_pv_input = 0.  # in units of degree
azim_pv_input = 0.
######################

data_path="/lustre/scratch/leiduan/data_needed/"
#data_path="/Volumes/Untitled/Project_CarbonPrice/Step3_cal_WindSolarCF_onMEMEX/"
fd_cam=os.listdir(data_path)

month_days={1:31,2:28,3:31,4:30,5:31,6:30,7:31,8:31,9:30,10:31,11:30,12:31}
month_days_array = np.array([31,28,31,30,31,30,31,31,30,31,30,31])

scf = MV.array(np.zeros([hour_in_years,361,540])); scf.id='scf'; scf.units='1'; scf.missing_value = 1e20

R_TAMB = 20.  # Reference ambient temperature (degC)
R_TMOD = 25.  # Reference module temperature (degC)
R_IRRADIANCE = 1000.  # Reference irradiance (W/m2)
HF_free = 0.035  # Free-standing module, assuming no wind
HF_inte = 0.05   # Building-integrated module

k_1 = -0.017162
k_2 = -0.040289
k_3 = -0.004681
k_4 = 0.000148
k_5 = 0.000169
k_6 = 0.000005

degree_to_radius = np.pi / 180.
radius_to_degree = 180. / np.pi

# pay attention to the defition of sign
def cal_solar_angles(lat, lon, days_ord, hr_idx):
    lat_radius = np.array(lat) * degree_to_radius
    lat_use = np.ones([361,540]) * lat_radius[:,None]
    
    sun_decli = 23.45 * np.sin(np.pi*2*(284+days_ord+1)/365) * degree_to_radius  # positive in NH, negative in SH
    ha_ew = -1 * np.tan(sun_decli) * np.tan(lat_radius)
    ha_ew_deg = np.zeros(361)
    ha_ew_deg[ (ha_ew>=-1)&(ha_ew<=1) ] = np.arccos( ha_ew[(ha_ew>=-1)&(ha_ew<=1)] ) * radius_to_degree
    ha_ew_deg[ (ha_ew<-1) ] = 180.
    ha_ew_deg[ (ha_ew >1) ] = 0.
    ha_ew_use_deg = np.ones([361,540]) * ha_ew_deg[:,None] # 0 -> always dark; 180 -> always day; 0-180 -> sunrise or sunset angle;
    
    local_time = np.zeros(540)
    lon_array = np.array(lon)
    local_time[lon_array<0.] = (24. + lon_array[lon_array<0.]/15. + hr_idx) + 0.5
    local_time[lon_array>=0.] = (lon_array[lon_array>=0.]/15 + hr_idx) + 0.5
    local_time[local_time>24] = local_time[local_time>24]-24.
    ha_deg = (local_time-12.)*(-15.)
    ha_rad = (local_time-12.)*(-15.) * degree_to_radius 
    ha_use_rad = np.ones([361,540]) * ha_rad
    ha_use_deg = np.ones([361,540]) * ha_deg # positive -> before noon; negative -> after noon
    
    cosine_zenith = np.sin(sun_decli) * np.sin(lat_use) + np.cos(sun_decli) * np.cos(lat_use) * np.cos(ha_use_rad)
    #cosine_zenith[cosine_zenith<0.] = 0.
    zenith_rad = np.arccos( cosine_zenith )
    #zenith_deg = zenith_rad * radius_to_degree
    sine_zenith = np.sin(zenith_rad)
    
    term1 = np.zeros([361,540]) # if 
    term2 = np.zeros([361,540]) # if facing toward the equator
    term3 = np.zeros([361,540]) 
    
    term0 = np.arcsin( (np.sin(ha_use_rad)*np.cos(sun_decli)) / (sine_zenith) ) #* mask1 * mask2
    term1[np.abs(ha_use_deg)<=ha_ew_use_deg] = 1.
    term1[np.abs(ha_use_deg) >ha_ew_use_deg] = -1.                
    term2[lat_use*(lat_use-sun_decli)>=0] = 1
    term2[lat_use*(lat_use-sun_decli) <0] = -1                
    term3[ha_use_rad>=0] = 1
    term3[ha_use_rad <0] = -1
    solar_azi_rad = term0*term1*term2 + (1-term1*term2)/2.*term3*np.pi
    #solar_azi_deg = solar_azi_rad * radius_to_degree
    
    return zenith_rad, solar_azi_rad, ha_rad


def cal_incidence_angles(zenith, solar_azi, tilt_pv, azim_pv, pv_type):
    max_azim_angle = 360.*degree_to_radius # for 1-axis vertical PV
    max_turn_angle = 45.*degree_to_radius # for 1-axis horizontal (with and without tilt)
    incidence_rad = np.zeros([361,540])
    if pv_type == 'f':
        print ('Using fixed tilt and azim solar panel')
        incidence_rad = np.arccos(np.cos(zenith)*np.cos(tilt_pv) + 
                                  np.sin(zenith)*np.sin(tilt_pv)*np.cos(solar_azi-azim_pv))
        
    elif pv_type == 'v':
        print ('Using 1-axis tracking panel, with vertical axis')
        incidence_rad = np.where( (solar_azi<=max_azim_angle)&(zenith<=tilt_pv)&(solar_azi>=0.), tilt_pv - zenith, incidence_rad) 
        incidence_rad = np.where( (solar_azi>=(-1)*max_azim_angle)&(zenith<=tilt_pv)&(solar_azi<=0.), tilt_pv - zenith, incidence_rad)
        incidence_rad = np.where( (solar_azi<=max_azim_angle)&(zenith>tilt_pv)&(solar_azi>=0.), zenith - tilt_pv, incidence_rad) 
        incidence_rad = np.where( (solar_azi>=(-1)*max_azim_angle)&(zenith>tilt_pv)&(solar_azi<=0.), zenith - tilt_pv, incidence_rad)
        incidence_rad = np.where( (solar_azi>max_azim_angle), np.arccos(np.cos(zenith)*np.cos(tilt_pv)+np.sin(zenith)*np.sin(tilt_pv)*np.cos(solar_azi-max_azim_angle)), incidence_rad )
        incidence_rad = np.where( (solar_azi<(-1)*max_azim_angle), np.arccos(np.cos(zenith)*np.cos(tilt_pv)+np.sin(zenith)*np.sin(tilt_pv)*np.cos(solar_azi-max_azim_angle)), incidence_rad )
        
    elif pv_type == 'h':
        print ('Using 1-axis tracking panel, with horizontal axis, no tilt, horizontal axis')
        # adjusted tilt is the tilt after the panel rotated on its axis 
        # the tilt is the angle between surface normal and panel surface normal
        adjusted_azim = np.zeros([361,540])
        adjusted_azim[ (solar_azi-azim_pv)>=0. ] = azim_pv[ (solar_azi-azim_pv)>=0. ] + np.pi/2.
        adjusted_azim[ (solar_azi-azim_pv) <0. ] = azim_pv[ (solar_azi-azim_pv) <0. ] - np.pi/2.
        adjusted_tilt = np.arctan(np.tan(zenith)*np.cos(adjusted_azim-solar_azi))
        adjusted_tilt[adjusted_tilt<0.] = adjusted_tilt[adjusted_tilt<0.]+np.pi
        adjusted_tilt[adjusted_tilt>max_turn_angle] = max_turn_angle
        adjusted_tilt[adjusted_tilt<max_turn_angle*(-1)] = max_turn_angle * (-1)
        incidence_rad = np.arccos(np.cos(zenith)*np.cos(adjusted_tilt) + np.sin(zenith)*np.sin(adjusted_tilt)*np.cos(solar_azi-adjusted_azim))
        
    elif pv_type == 'ht':
        print ('Using 1-axis tracking panel, with horizontal axis, with tilt')
        incidence_rad_tmp = np.arccos(np.cos(zenith)*np.cos(tilt_pv) + np.sin(zenith)*np.sin(tilt_pv)*np.cos(solar_azi-azim_pv))
        term1 = np.zeros([361,540])
        term2 = np.zeros([361,540])
        term0 = azim_pv + np.arctan( np.sin(zenith)*np.sin(solar_azi-azim_pv)/np.cos(incidence_rad_tmp)/np.sin(tilt_pv) )
        term1[(term0-azim_pv)*(solar_azi-azim_pv)>=0.] = 0.
        term1[(term0-azim_pv)*(solar_azi-azim_pv) <0.] = 1.
        term2[(solar_azi-azim_pv)>=0.] = 1.
        term2[(solar_azi-azim_pv) <0.] = -1.
        adjusted_azim = term0 + term1*term2*np.pi
        
        term00 = np.arctan(np.tan(tilt_pv)/np.cos(term0-azim_pv))
        term11 = np.zeros([361,540])
        term11[term00>=0] = 0.
        term11[term00 <0] = 1.
        adjusted_tilt = term00+term11*np.pi 
        # cos_adjusted_tilt = cos_adjust_tilt * cos_turn_angle
        incidence_rad = np.arccos(np.cos(zenith)*np.cos(adjusted_tilt) + np.sin(zenith)*np.sin(adjusted_tilt)*np.cos(solar_azi-adjusted_azim))
        
    elif pv_type == '2a':
        print ('Using 2-axis tracking panela')
        incidence_rad = 0.
    #incidence_deg = incidence_rad * radius_to_degree
    return incidence_rad

def replace_nan(var):
    whereisNan = np.isnan(var)
    var[whereisNan] = 0.
    return var

tilt_pv = np.zeros([361,540])+tilt_pv_input*degree_to_radius
azim_pv = np.zeros([361,540])+azim_pv_input*degree_to_radius

df = np.zeros([361,540])
count_num=1
for file in fd_cam:
        if file[:-7] == case_name and file[-3:]=='.nc':
            f=cdms.open(data_path+file)
            month = int(file[-7:-5])
            days = int(file[-5:-3])
            if month == 1:
                days_ord = 0 + days
                position = (0 + (days-1)) *24
            else:
                days_ord = np.sum(month_days_array[:month-1]) + days
                position = (np.sum(month_days_array[:month-1]) + (days-1)) *24
            print (month, days, days_ord, position)
            
            # get data
            lat = f.getAxis('YDim_EOSGRID')
            lon = f.getAxis('XDim_EOSGRID')
            swgdn_tmp = MV.filled(f('SWGDN',squeeze=1),0.) #(24,361,540)
            swtdn_tmp = MV.filled(f('SWTDN',squeeze=1),0.) #(24,361,540)
            t = f('TS',squeeze=1)-273.15
            
            # separate direct and diffuse, based on Erbs et al., 1982 and https://pvlib-python.readthedocs.io/en/latest/
            kt = np.zeros([24,361,540])
            kt[swtdn_tmp != 0.] = (swgdn_tmp[swtdn_tmp != 0.]/swtdn_tmp[swtdn_tmp != 0.])
            kt[kt<0.] = 0.
            df = np.zeros([24, 361,540]) # error1
            df = np.where( kt<= 0.22, 1 - 0.09 * kt, df)
            df = np.where((kt > 0.22) & (kt <= 0.8), 0.9511 - 0.1604*kt + 4.388*kt**2 - 16.638*kt**3 + 12.336*kt**4, df)
            df = np.where( kt > 0.8, 0.165, df)  # 
            df = np.where( kt== 0.0, 1.0, df)    # where no TOA SW, all diffuse flux
            
            dhi = df * swgdn_tmp                            # diffuse radiation
            dni = (swgdn_tmp - dhi)                         # direct radiation            
            
            for hr_idx in range(24):
                zenith, solar_azi, ha = cal_solar_angles(lat, lon, days_ord, hr_idx)                
                mask1 = MV.filled(MV.masked_equal(swtdn_tmp[hr_idx],0)*0+1,0)
                mask2 = MV.filled(MV.masked_greater_equal(zenith,np.pi/2)*0+1,0)
                
                incidence_rad = cal_incidence_angles(zenith, solar_azi, tilt_pv, azim_pv, 'h') #Based on Braun and Mitchell, 1983
                cosine_zenith = np.cos(zenith)* mask1*mask2
                cosine_incide = np.cos(incidence_rad) * mask1*mask2
                adjust_factor_dni = replace_nan(cosine_incide / cosine_zenith) 
                adjust_factor_dni[(adjust_factor_dni<1.)]=1.
                dni_adjust = dni[hr_idx] * adjust_factor_dni  
                dhi_adjust = dhi[hr_idx]*(1+np.cos(ha))/2. + 0.3*(dni[hr_idx]+dhi[hr_idx])*(1-np.cos(ha))/2.
                rad_adjust = replace_nan(dni_adjust + dhi_adjust)
                
                # Now we use the beer-lambert law to avoid extremely small cosine zenith angle
                maximum_index = np.argwhere( swgdn_tmp[hr_idx]== np.max(swgdn_tmp[hr_idx]) )
                base = swgdn_tmp[hr_idx][maximum_index[0,0]][maximum_index[0,1]]/swtdn_tmp[hr_idx][maximum_index[0,0]][maximum_index[0,1]]
                cons = swtdn_tmp[hr_idx][maximum_index[0,0]][maximum_index[0,1]]
                potential_max_solar = np.zeros([361,540])
                potential_max_solar = np.where(cosine_zenith!=0.,
                                               cons*base**(1.0/cosine_zenith),
                                               potential_max_solar)
                
                rad_adjust = np.minimum(rad_adjust, potential_max_solar)
                
                T_ = np.array((1*t[hr_idx] + HF_free * rad_adjust) - R_TMOD )
                G_ = np.array(rad_adjust / R_IRRADIANCE)

                index = int(position+hr_idx)
                scf[index,G_==0.] = 0.
                scf[index,G_ >0.] = G_[G_>0.] * (1 + k_1 * np.log(G_[G_>0.]) + k_2 * (np.log(G_[G_>0.])) ** 2 +
                   T_[G_>0.] * (k_3 + k_4 * (np.log(G_[G_>0.])) + k_5 * (np.log(G_[G_>0.])) ** 2) + k_6 * (T_[G_>0.] ** 2))

# use NetCDF3 Classic format
cdms.setNetcdfShuffleFlag(0)      # netcdf3 classic...
cdms.setNetcdfDeflateFlag(0)      # netcdf3 classic...
cdms.setNetcdfDeflateLevelFlag(0) # netcdf3 classic...

fout=cdms.open(case_name+'_scf_duanleiNEW2.nc','w')
fout.write(scf)
fout.close()

scf_annual = cdutil.averager(scf,axis=0,weights='equal')
scf_annual.id='scf_annual'
gout=cdms.open(case_name+'_scf_annual_duanleiNEW2.nc','w')
gout.write(scf_annual)
gout.close()
