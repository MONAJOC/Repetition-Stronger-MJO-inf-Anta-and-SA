;Figure a: TRMM rainfall, OLR anomaly and MJO index Amp.
;Figure b: psi200 anomaly and slp anomaly.

lats = -15.
latn = -5.
lonw = 180.
lone = 210.
pi   = atan(1.0)*4.


DATADIR = "/home/yinqc/SMJO/Data/"

olrf = addfile(DATADIR+"OLR/olr.daily.anom.201503.nc","r")
raif = addfile(DATADIR+"TRMM/3B42RT_Daily.201503.7.nc4.nc","R")

;MJO index and amplitude in csv files
RMMf  = DATADIR+"rmm.74toRealtime.csv"
delim = ","
rmmf  = asciiread(RMMf,-1,"string")
iyear = str_get_field(rmmf,1,delim)
imon  = str_get_field(rmmf,2,delim)
rmm1  = str_get_field(rmmf,4,delim)
rmm2  = str_get_field(rmmf,5,delim)
phase = str_get_field(rmmf,6,delim)
amp   = str_get_field(rmmf,7,delim)

iyear := toint(iyear(1::))
imon  := toint(imon(1::))
rmm1  := tofloat(rmm1(1::))
rmm2  := tofloat(rmm2(1::))
phase := toint(phase(1::))
amp   := tofloat(amp(1::))

nt = ind(iyear.eq.2015.and.imon.eq.3)
rmm1  := rmm1(nt)
rmm2  := rmm2(nt)
phase := phase(nt)
amp   := amp(nt)
delete([/rmmf,nt/])

;reading data
olr  = olrf->olr(:,{latn:lats},{lonw:lone})
rad  = olr&lat*pi/180.
clat = cos(rad)
olr_mean = wgt_areaave_Wrap(olr, clat, 1.0, 1)

pre   = raif->precipitation(time|:,lat|:,lon|:)
pre  := pre(:,{-15.125:-4.875},{-179.875:-149.875})
radp  = pre&lat*pi/180.
clatp = cos(radp)
pre_mean = wgt_areaave_Wrap(pre, clatp, 1.0, 1)

time  = pre&time
;plot
;resource
resL = True
resL@gsnDraw  = False
resL@gsnFrame = False
resL@gsnMaximize = True
resL@vpXF      = 0.05
resL@vpYF      = 0.85
resL@vpWidthF  = 0.9
resL@vpHeightF = 0.5
resL@gsnLeftString    = "(a)"
resL@gsnLeftStringOrthogonalPosF = -0.005
resL@gsnLeftStringFontHeightF = 0.03
resL@tiYAxisString    = "MJO Intensity RMM Index"
resL@tiYAxisFontColor = "deepskyblue"
resL@tiYAxisFontHeightF = 0.015
resL@tmYLMode        = "Manual"
resL@tmYLMinorPerMajor = 0
resL@tmYLTickSpacingF  = 1
resL@tmYLLabelFontHeightF = 0.015
resL@tmXBMode        = "Manual"
resL@tmXBMinorPerMajor = 0
resL@tmXBTickSpacingF  = 5
resL@tmXBLabelFontHeightF = 0.015
resL@tmXTMinorPerMajor = 0
resL@tmYLMajorLineColor = "deepskyblue"
resL@tmYLLabelFontColor = "deepskyblue"
resL@tmYLMajorLengthF   = -0.01
resL@tmXBMajorLengthF   = -0.01
resL@tmXTMajorLengthF   = 0.0
resL@trYMaxF  = 5
resL@trYMinF  = 0
resL@trXMaxF  = 30
resL@trXMinF  = 0
resL@tmXBMode = "Explicit"
x = ispan(1,31,1)
resL@tmXBValues = x-1
resL@tmXBLabels = x

resR = True
resR@gsnDraw  = False
resR@gsnFrame = False
resR@xyLineThicknessF = 6.0
resR@xyLineColors     = "darkorange2"
resR@trYMaxF  =  20
resR@trYMinF  = -80
resR@tiYAxisFontColor = "darkorange2"
resR@tiYAxisString    = "OLR Anomlay [W/m~S~2~N~]"
resR@tmYRMinorPerMajor = 0
resR@tmYRLabelFontHeightF = 0.015
resR@tmXBLabelFontHeightF = 0.01
resR@tiYAxisFontHeightF = 0.015
resR@tmYRMajorLengthF  = -0.01
resR@tmYRMajorLineColor = "darkorange2"
resR@tmYRLabelFontColor = "darkorange2"

resp = resR
resp@xyLineColors     = "black"
resp@trYMaxF  = 30
resp@trYMinF  =	0
resp@tiYAxisFontColor = "black"
resp@tiYAxisString    = "Precipitation [mm/day]"
resp@amOrthogonalPosF = 0.63
resp@tmYRLabelFontHeightF = 0.015
resp@tmYRMinorPerMajor = 0
resp@tmYRMajorLengthF  = -0.01
resp@tmYRMajorLineColor = "black"
resp@tmYRLabelFontColor = "black"

;colors
colors = (/"indianred3","darkorange2","darkgoldenrod2","lightgoldenrod2","lightgoldenrodyellow","darkseagreen3","darkturquoise","deepskyblue3"/)
colorf = new(dimsizes(time),string)
do i = 0, 30
    colorf(i) = colors(phase(i)-1)
end do    

lbres = True
lbres@lbAutoManage = False
lbres@lbFillColors = colors
lbres@lbMonoFillPattern  = True
lbres@lbLabelFontHeightF = 0.02
lbres@lbOrientation      = "Horizontal"
lbres@lbPerimOn          = False
lbres@lbTitleOn = True
lbres@lbTitleFontHeightF = 0.02
lbres@lbTitleString = "MJO phase"
lbres@lbTitleDirection  = "Across"
lbres@lbBoxMinorExtentF = 0.2
lbres@lbTitlePosition = "Left"
lbres@lbLabelOffsetF  = -0.14
lbres@vpWidthF  = 0.5
lbres@vpHeightF = 0.25
lbres@lbBoxLinesOn = False
lbres@lbFillOpacityF = 0.8


type = "eps"
;type@wkPaperWidthF   = 11.7
;type@wkPaperHeightF  = 8.3
;type@wkOrientation   = "landscape"
wks = gsn_open_wks(type,"/home/yinqc/SMJO/plot/figure2")
dum  = new(31,graphic)
tmp  = new(dimsizes(amp),float,-999)
plot = gsn_csm_xy3(wks,x-1,tmp,olr_mean,pre_mean,resL,resR,resp)

pmres = True
pmres@gsMarkerIndex = 16
pmres@gsMarkerSizeF = 20.
pmres@gsMarkerOpacityF = 0.8
do i = 0,30
    pmres@gsMarkerColor = colorf(i)
    dum(i) = gsn_add_polymarker(wks,plot,i,amp(i),pmres)
end do
lblb = ispan(1,8,1)
lbid = gsn_create_labelbar(wks,8,lblb,lbres)
amres = True
amres@amJust = "TopCenter"
amres@amParallelPosF   = 0.0
amres@amOrthogonalPosF = -0.9
annoid = gsn_add_annotation(plot,lbid,amres)

draw(plot)







