
diro1 = "/home/yinqc/SMJO/Data/CFSR/"
DATAN = "uwind"
diro2 = "/home/yinqc/SMJO/Data/CFSv2/"

varclimf = addfile(diro1+"uwind"+"/"+DATAN+".gdas.1979-2010.grb2.nc", "r")

time   =  varclimf->time
TIME   =  cd_calendar(time, 0)
mon    =  toint(TIME(:,1))
 
var_mon  = varclimf->U_GRD_L100_Avg
lat   = var_mon&lat
lon   = var_mon&lon
ntime = dimsizes(time)/2.
;retain the first 3 harmonics.
var_fft = ezfftf_n(var_mon, 0)
var_fft(:,3::,:,:) = 0.
var_re = ezfftb_n(var_fft, 0., 1)
copy_VarCoords(var_mon, var_re)
nt       = ind(mon .eq. 3)	
var_3m   = dim_avg_n_Wrap(var_re(nt,:,:),0)
var_mean = dim_avg_n_Wrap(var_mon,0)
varm = var_3m+var_mean
copy_VarCoords(var_mon(0,:,:),varm)
system("/bin/rm -f "+diro2+DATAN+".clim.201503.nc")
ncdf = addfile(diro2+DATAN+".clim.201503.nc","c")
ncdf->uwind = varm

