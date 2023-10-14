;This script calculate 2-D wave vector. reference: TN 2001.
;uwnd_mon: mean values of uwind. (nmon,nlat,nlon)
;vwnd_mon: mean values of vwind. (nmon,nlat,nlon)
;uwnd_mean: mean values of uwind at March.
;vwnd_mean: mean values of vwind at March.
;psi_anom: daily values of psi anomalies.
;magn_wind: maganitude of uwind and vwind March mean.
;dphidx, dphidy ... : partitial of psi200-hPa.
;Wx,Wy: Wave vector in lonigtude and latitude.
;constant
pi = atan(1.0)*4.
r  = 6371000.


DATADIR="/home/yinqc/SMJO/Data/"
phif = addfile(DATADIR+"CFSv2/psi200.anom.201503.nc","r")
uwnf = addfile(DATADIR+"CFSR/uwind/uwind.gdas.1979-2010.grb2.nc","r")
vwnf = addfile(DATADIR+"CFSR/vwind/vwind.gdas.1979-2010.grb2.nc","r")

;climatology variables including uwind and vwind.
uwnd_mon = uwnf->U_GRD_L100_Avg
vwnd_mon = vwnf->V_GRD_L100_Avg

;calculate climatology
ufft = ezfftf_n(uwnd_mon,0)
ufft(:,3::,:,:) = 0.
uwnd_re = ezfftb_n(ufft,0.,1)
uwnd_m  = dim_avg_n_Wrap(uwnd_mon, 0)

vfft = ezfftf_n(vwnd_mon,0)
vfft(:,3::,:,:) = 0.
vwnd_re = ezfftb_n(vfft,0.,1)
vwnd_m  = dim_avg_n_Wrap(vwnd_mon, 0)


time = uwnd_mon&time
TIME = cd_calendar(time, 0)
year = toint(TIME(:,0))
mon  = toint(TIME(:,1))
nt   = ind(mon.eq.3)
do i = 0, dimsizes(time)-1
    uwnd_re(i,:,:) = uwnd_re(i,:,:) + uwnd_m
    vwnd_re(i,:,:) = vwnd_re(i,:,:) + vwnd_m
end do
uwnd_mean = dim_avg_n_Wrap(uwnd_re(nt,:,:),0)
vwnd_mean = dim_avg_n_Wrap(vwnd_re(nt,:,:),0)
delete([/time,TIME,year,mon,nt/])

;psi 200hPa
psi_hour = phif->psi
psi_anom = calculate_daily_values(psi_hour,"avg",0,False)
lat  = psi_anom&lat
lon  = psi_anom&lon
nlat = num(lat)
nlon = num(lon)
time = psi_anom&time
TIME = cd_calendar(time,0)
year = toint( TIME(:,0) )
mon  = toint( TIME(:,1) )
day  = toint( TIME(:,2) )

;calculate horizontal wind mganitude.
magn_wind = wind_speed(uwnd_mean,vwnd_mean)
copy_VarCoords(uwnd_mean,magn_wind)

;calculate dphidx, dphidy, dphi2dx2, dphi2dxdy, dphi2dy2
dlon   = abs(lon(1) - lon(2))
dlat   = abs(lat(1) - lat(2))
dX     = r*cos(lat*pi/180)*pi/180.*dlon
dY     = r*dlat*pi/180.

dphidx = new((/4,nlat+1,nlon+1/),float)
copy_VarCoords(psi_anom(0,:,:), dphidx(0,:,:))
dphidx    = 0.
dphidy    = dphidx
dphi2dx2  = dphidx
dphi2dxdy = dphidx
dphi2dy2  = dphidx
do i = 0,3
    ihead   = i*3 + 16
    iend    = i*3 + 18 
    nt      = ind( day.ge.ihead.and.day.le.iend )
    phi_tmp = dim_avg_n_Wrap(psi_anom(nt,:,:),0)
    do j = 0, nlat-1
        dphidx(i,j:j,:)    = center_finite_diff_n( phi_tmp(j:j,:),dX(j),True,0,1)
        dphi2dx2(i,j:j,:)  = center_finite_diff_n( dphidx(i,j:j,:),dX(j),True,0,1)
    end do
    dphidy(i,:,:)    = center_finite_diff_n( phi_tmp,dY,False,0,0)
    dphi2dxdy(i,:,:) = center_finite_diff_n( dphidx(i,:,:),dY,False,0,0)
    dphi2dy2(i,:,:)  = center_finite_diff_n( dphidy(i,:,:),dY,False,0,0)
    delete([/ihead,iend,nt,phi_tmp/])
end do

;calculate wave vector.
Wx = dphidx
Wx = 0.
Wy = Wx

tm1 = Wx
tm2 = Wx
tm3 = Wx
tm4 = Wx
do i = 0,3
    ihead   = i*3 + 16
    iend    = i*3 + 18 
    nt      = ind( day.ge.ihead.and.day.le.iend )
    phi_tmp = dim_avg_n_Wrap(psi_anom(nt,:,:),0)
    tm1(i,:,:) = uwnd_mean*( dphidx(i,:,:)+dphidx(i,:,:) - phi_tmp*dphi2dx2(i,:,:) )/(2.*magn_wind)
    tm2(i,:,:) = vwnd_mean*( dphidx(i,:,:)*dphidy(i,:,:) - phi_tmp*dphi2dxdy(i,:,:) )/(2.*magn_wind)
    tm3(i,:,:) = uwnd_mean*( dphidx(i,:,:)*dphidy(i,:,:) - phi_tmp*dphi2dxdy(i,:,:) )/(2.*magn_wind)
    tm4(i,:,:) = vwnd_mean*( dphidx(i,:,:)*dphidx(i,:,:) - phi_tmp*dphi2dy2(i,:,: ))/(2.*magn_wind)
    delete([/ihead,iend,nt,phi_tmp/])
end do

Wx = tm1 + tm2
Wy = tm3 + tm4
Wx@long_name = "Wave vector in longitude"
Wy@long_name = "Wave vector in latitude"
Wx@unit = "m^2s^-2"
Wy@unit = "m^2s^-2"

;mask where u less than 5 m/s.
cuvar = conform_dims(dimsizes(Wx), uwnd_mean, (/1,2/))
Wx = mask(Wx, cuvar.lt.5,False)
Wy = mask(Wy, cuvar.lt.5,False)
;mask 10S-10N.
Wx(:,{-10:10},:) = Wx@_FillValue
Wy(:,{-10:10},:) = Wy@_FillValue

;write to netcdf file.
diro = DATADIR
filo = "W_vector"+".nc"
system("/bin/rm -f "+diro+filo)
ncdf = addfile(diro+filo,"c")
ncdf->Wx = Wx
ncdf->Wy = Wy

