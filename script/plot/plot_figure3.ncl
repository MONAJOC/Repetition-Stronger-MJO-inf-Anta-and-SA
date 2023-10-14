;This plot the figure3.
;(a) gph200 average among March 22-24.
;(b) Ks.
;Information
;Date: 2023-10-05
load "/home/yinqc/script/raytracing.ncl"
DATADIR = "/home/yinqc/SMJO/Data"

gphf = addfile(DATADIR+"/CFSv2/ghp200.anom.201503.nc","r")
uwnf = addfile(DATADIR+"/CFSR/uwind/uwind.gdas.1979-2010.grb2.nc","r")

;calculate mean uwind.
time   =  uwnf->time
TIME   =  cd_calendar(time, 0)
mon    =  toint(TIME(:,1))
var_mon  = uwnf->U_GRD_L100_Avg
lat = var_mon&lat
lon = var_mon&lon
ntime = dimsizes(time)/2.
var_fft = ezfftf_n(var_mon, 0)
var_fft(:,3::,:,:) = 0.
var_re = ezfftb_n(var_fft, 0., 1)
copy_VarCoords(var_mon, var_re)
nt       = ind(mon .eq. 3)	
var_3m   = dim_avg_n_Wrap(var_re(nt,:,:),0)
var_mean = dim_avg_n_Wrap(var_mon,0)
uwmean   = var_3m + var_mean
uwmean   = dim_avg_n_Wrap(var_mon(nt,:,:),0)
copy_VarCoords(var_mon(0,:,:), uwmean)
delete([/time,TIME,mon,var_mon,lat,lon,ntime,var_fft,var_re,nt,var_3m,var_mean/])

;Main code.
gph  = gphf->gph
gph := dim_avg_n_Wrap(gph,0)

;caculate Ks.
;;calculate beta terms in Mercator
pi     = atan(1.0)*4.
a      = 6371000.
omega  = 2*pi/(24*60*60)
lat    = uwmean&lat
lon    = uwmean&lon
nlat   = dimsizes(lat)
phirad = lat*pi/180.
phi  = conform_dims(dimsizes(uwmean), phirad, 0)
YY   = phi
do i = 0, nlat-1
    YY(i,:)   = a*log((1.+sin(phi(i,:)))/cos(phi(i,:))) 
end do
; Mercator coordinates         
umz  = uwmean/cos(phi)                      ; uwind in Mercator.
umz(0,:) = 0.
umz(nlat-1,:) = 0.
copy_VarCoords(uwmean, umz)
b = (1.+cos(2*phi))/2.                     ; cos^2
d = new(dimsizes(b),float)                 ; 1/cos^2
do i=1,nlat-2
    d(i,:) = 1/b(i,:)
end do
c       = umz*b
dfdy    = c
tm2     = c
d2udy2  = c
;1st term: df/dy
dfdy = 2*omega*b/a
;2nd term: d2u/dy2
tm2  = center_finite_diff_n(c, YY, False, 0, 0)
tm2 := tm2*d
d2udy2 = center_finite_diff_n(tm2, YY, False, 0, 0)
;betam = df/dy - d2u/dy2
beta_m = dfdy - d2udy2
copy_VarCoords(uwmean, beta_m)
delete([/dfdy,tm2,d2udy2/])

Kst = new(dimsizes(uwmean),float,1.e20)
do i = 0, nlat-1
    do j = 0, dimsizes(lon)-1
        Kst(i,j) = sqrt(abs(beta_m(i,j))*cos(lat(i)*pi/180.)/abs(uwmean(i,j)))
    end do
end do
Kst  := mask(Kst,(umz*beta_m).le.0,False)
Ks  = Kst*6371000.
Ks := smth9_Wrap(Ks, 0.5, -0.25, 0)
Ks := smth9_Wrap(Ks, 0.5, -0.25, 0)
copy_VarCoords(gph,Ks)
;Ks := a*cos(phi)*Ks
;Plot.
res = True
res@gsnFrame           = False
res@gsnDraw            = False
res@gsnMaximize        = True
res@gsnAddCyclic       = True
res@gsnLeftString      = ""
res@gsnRightString     = ""
res@gsnCenterStringFontHeightF = 0.025
res@mpFillOn        = False
res@mpProjection    = "Satellite"
res@mpCenterLonF    = 250.0       
res@mpCenterLatF    = -35.        
res@mpLandFillColor = "white"
res@cnFillOn             = True
res@cnLinesOn            = False
res@cnLineLabelsOn       = False
res@lbLabelStride        = 8
;contour
;contour plot
resc = True
resc@gsnFrame     = False
resc@gsnDraw      = False
resc@gsnMaximize  = False
resc@gsnAddCyclic = True
resc@gsnContourNegLineDashPattern = 1
resc@gsnLeftString      = ""
resc@gsnRightString     = ""
resc@cnMonoLineColor    = True
resc@cnLinesOn          = True
resc@cnLineLabelsOn     = False
resc@cnLineColor        = "black"
resc@cnLineThicknessF   = 1.0
resc@cnInfoLabelPerimOn = False
resc@cnInfoLabelOn      = False
resc@cnLevelSelectionMode = "Explicitlevels"
resc@cnLevels      = (/-100.,-50.,-150.,200.,300.,400./)

plot    = new(2,graphic)
plotc   = plot(0)
plotr   = plot(0)

type = "eps"

wks = gsn_open_wks(type,"/home/yinqc/SMJO/plot/Figure.3")
Label = (/"(a) Geopotential Height Anomaly 200 hPa","(b) Stationary Wavenumber K~B~s~N~ 200 hPa"/)

res@gsnCenterString = Label(0)
res@lbLabelBarOn      = True
res@cnLabelBarEndStyle = "ExcludeOuterBoxes"
res@lbOrientation     = "Vertical"
res@lbTitleString     = "[gpm]"
res@lbTitlePosition   = "Right"
res@lbTitleFontHeightF= 0.015
res@lbBoxLinesOn      = False
res@mpSatelliteDistF  = 0.5        
res@mpGridAndLimbDrawOrder = "PostDraw"
res@mpGridLineColor        = "Gray60"
res@mpGridLineThicknessF   = 1.0
res@cnFillPalette        = "nrl_sirkes"
res@cnLevelSelectionMode = "Manuallevels"
res@mpGridLineDashPattern      = 2          
res@mpGridAndLimbOn            = True
res@cnMinLevelValF       = -600.
res@cnMaxLevelValF       =  600.
res@cnLevelSpacingF      =  25.
plot(0) = gsn_csm_contour_map(wks, gph, res)
plotc   = gsn_csm_contour(wks,gph,resc)
overlay(plot(0),plotc)
;Raytracing
;f        = asciiread("/home/yinqc/script/Ray.csv", -1, "string")
;delim = ","
;x = str_get_field(f,2,delim)
;y = str_get_field(f,1,delim)
xlist = raytracing(uwmean,-160.,-15.,1000.,15.,3,0.)
ypt = xlist(0,:)            
xpt = xlist(1,:)+360.   
pres                  = True              
pres@gsLineThicknessF = 3.0             
dum = gsn_add_polyline(wks,plot(0),xpt,ypt,pres)     

res@lbLabelStride        = 1
res@gsnCenterString = Label(1)
res@lbLabelBarOn      = True
res@lbBoxLinesOn      = False
res@mpGridAndLimbDrawOrder = "PostDraw"
res@mpGridLineColor        = "Gray60"
res@mpGridAndLimbOn        = False
res@mpGridLineThicknessF   = 1.0
res@cnFillPalette        = "BlueDarkRed18"
res@cnLevelSelectionMode = "Manuallevels"
res@cnMinLevelValF       = 0.0
res@cnMaxLevelValF       = 10.0
res@cnLevelSpacingF      = 1.0
res@lbTitleOn            = False
plot(1) = gsn_csm_contour_map(wks, Ks, res)


resP   = True
resP@gsnMaximize = True
resP@gsnFrame    = False
resP@gsnDraw     = True
resP@gsnPanelLabelBar = False
resP@gsnPanelTop    = 0.98
resP@gsnPanelBottom = 0.02
resP@gsnPanelRowSpec   = True
gsn_panel(wks,plot,(/2/),resP)
delete(wks)
