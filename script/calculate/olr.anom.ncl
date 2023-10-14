;caculate olr anomalies by subtracting daily climatology from 1979 to 2010.


olr_climf = addfile("/home/yinqc/SMJO/Data/OLR/OLR.daily.mean.1979-1995.nc","r")
olr_dayf  = addfile("/home/yinqc/SMJO/Data/OLR/OLR.uninterpolated.daily.201503.nc","r")

time  = olr_climf->time
TIME    = cd_calendar(time,0)
year    = toint( TIME(:,0) )
mon     = toint( TIME(:,1) )
day     = toint( TIME(:,2) )
ddd     = day_of_year(year, mon, day)
yyyyddd = year*1000 + ddd

olr_clim = olr_climf->olr
olr_day  = olr_dayf->olr

lat = olr_day&lat
lon = olr_day&lon

;caculate daily mean.
olr_mean  = clmDayTLL(olr_clim, yyyyddd)
;caculate daily anomaly
olr_anom  = olr_day
do i = 61,91,1
   olr_anom(i-61,:,:) = olr_day(i-61,:,:) - olr_mean(i,:,:)
end do
copy_VarCoords(olr_day,olr_anom)
diro = "/home/yinqc/SMJO/Data/OLR/"
filo = "olr.daily.anom.201503"+".nc"
system("/bin/rm -f  "+diro+filo)
ncdf = addfile(diro+filo,"c")
ncdf->olr = olr_anom
