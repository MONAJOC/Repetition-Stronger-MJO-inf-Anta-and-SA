
diro1 = "/home/yinqc/SMJO/Data/CFSR/"
DATAN = "ghp200"
diro2 = "/home/yinqc/SMJO/Data/CFSv2/"

varclimf = addfile(diro1+"gph200"+"/"+DATAN+".gdas.1979-2010.grb2.nc", "r")
vardayf  = addfile(diro2+"gph200"+".cdas1.201503.grb2.nc", "r")

time   =  varclimf->time
TIME   =  cd_calendar(time, 0)
year   =  toint(TIME(:,0))
mon    =  toint(TIME(:,1))
 
var_mon  = varclimf->HGT_L100_Avg
var_hour = vardayf->HGT_L100

lat = var_hour&lat
lon = var_hour&lon
ntime = dimsizes(time)/2.
;retain the first 3 harmonics.
var_fft = ezfftf_n(var_mon, 0)
var_fft(:,3::,:,:) = 0.
var_re = ezfftb_n(var_fft, 0., 1)
copy_VarCoords(var_mon, var_re)
nt       = ind(mon .eq. 3)	
var_3m   = dim_avg_n_Wrap(var_re(nt,:,:),0)
var_mean = dim_avg_n_Wrap(var_mon(nt,:,:),0)
;calc anomalies.
var_anom = var_hour
time := var_hour&time
ntime = dimsizes(time)
do i = 0,ntime-1
    var_anom(i,:,:) = var_hour(i,:,:) - var_mean
end do
copy_VarCoords(var_hour, var_anom)

system("/bin/rm -f "+diro2+DATAN+".anom.201503.nc")
ncdf = addfile(diro2+DATAN+".anom.201503.nc","c")
ncdf->gph = var_anom

