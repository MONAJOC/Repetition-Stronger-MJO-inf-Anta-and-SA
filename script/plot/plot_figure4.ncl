;This script is uesed to plot olr(shading), wave vector(vectors) and psi(contour).

DATADIR = "/home/yinqc/SMJO/Data"

pwf   = addfile(DATADIR+"/CFSv2/p_water.20150322-25.nc","r")
uvf   = addfile(DATADIR+"/CFSv2/uv.20150322-25.nc","R")
pvf   = addfile(DATADIR+"/CFSv2/pv.20150322-25.nc","r")

;primary 
pv  = calculate_daily_values(pvf->PVORT_L107,"avg",0,False)
pw  = calculate_daily_values(pwf->P_WAT_L200,"avg",0,False)
uw  = calculate_daily_values(uvf->U_GRD_L100,"avg",0,False)
vw  = calculate_daily_values(uvf->V_GRD_L100,"avg",0,False)

time = pw&time
TIME = cd_calendar(time, 0)
day  = toint(TIME(:,2))

;calculate 3-day mean anomaly to plot.
pv := smth9_Wrap(pv,0.50,0.25,False)
pv := smth9_Wrap(pv,0.50,0.25,False)
pv := smth9_Wrap(pv,0.50,0.25,False)
pv := smth9_Wrap(pv,0.50,0.25,False)
pv := smth9_Wrap(pv,0.50,0.25,False)
pv := smth9_Wrap(pv,0.50,0.25,False)
pv := smth9_Wrap(pv,0.50,0.25,False)
pv := smth9_Wrap(pv,0.50,0.25,False)
pv := smth9_Wrap(pv,0.50,0.25,False)
pv := smth9_Wrap(pv,0.50,0.25,False)
pv := smth9_Wrap(pv,0.50,0.25,False)
pv := smth9_Wrap(pv,0.50,0.25,False)
pv := smth9_Wrap(pv,0.50,0.25,False)
pv := smth9_Wrap(pv,0.50,0.25,False)
pv := smth9_Wrap(pv,0.50,0.25,False)
pv := smth9_Wrap(pv,0.50,0.25,False)
pv := smth9_Wrap(pv,0.50,0.25,False)

pv = pv*10^(6)
delete([/pw@standard_name,pw@units,uw@standard_name,uw@units,vw@standard_name,vw@units,pv@long_name,pv@units/])
;set resources
;shading plot
resolr = True
resolr@gsnFrame           = False
resolr@gsnDraw            = False
resolr@gsnMaximize        = True
resolr@gsnAddCyclic       = True
resolr@gsnMajorLatSpacing = 15
resolr@gsnMajorLonSpacing = 30
resolr@gsnCenterStringFontHeightF = 0.025
resolr@mpFillOn        = False
resolr@mpFillDrawOrder = "Postdraw"
;resolr@mpProjection= "Robinson"
resolr@mpMinLatF   = -85.
resolr@mpMaxLatF   = 20.
resolr@mpMinLonF   = 180.
resolr@mpMaxLonF   = 345.
resolr@mpLandFillColor = "white"
resolr@mpGridAndLimbDrawOrder = "PostDraw"
resolr@mpGridLineColor        = "Gray60"
resolr@mpGridLineThicknessF   = 1.0
resolr@cnFillOn             = True
;resolr@cnFillMode           = 
resolr@cnFillPalette        = "precip2_17lev"
resolr@cnFillOpacityF       = 1.0
resolr@cnLinesOn            = False
resolr@cnLineLabelsOn       = False
resolr@cnLevelSelectionMode = "Manuallevels"
resolr@cnMinLevelValF       =   0.
resolr@cnMaxLevelValF       =  60.
;resolr@cnLevelSpacingF     =  
resolr@lbLabelBarOn   = False
resolr@cnLabelBarEndStyle = "ExcludeOuterBoxes"

;contour plot
respsi = True
respsi@gsnFrame     = False
respsi@gsnDraw      = False
respsi@gsnMaximize  = False
respsi@gsnAddCyclic = True
respsi@gsnContourNegLineDashPattern = 0
respsi@cnMonoLineColor    = True
respsi@cnLinesOn          = True
respsi@cnLineLabelsOn     = False
respsi@cnLineColor        = "black"
respsi@cnLineThicknessF   = 2.5
respsi@cnInfoLabelPerimOn = False
respsi@cnInfoLabelOn      = False
respsi@cnLevelSelectionMode = "Explicitlevels"
respsi@cnLevels      = (/-2.0/)

;vector plot.
resw = True
resw@gsnFrame     = False
resw@gsnDraw      = False
resw@gsnMaximize  = False
resw@gsnAddCyclic = False
resw@vcGlyphStyle            = "LineArrow"
resw@vcLineArrowThicknessF   = 0.8
resw@vcMinDistanceF          = 0.015
resw@vcRefMagnitudeF         = 20.
resw@vcRefLengthF            = 0.05
resw@vcRefAnnoOrthogonalPosF = -1.0
resw@vcRefAnnoParallelPosF   = 1.18
resw@vcLineArrowColor        = "black"
resw@vcRefAnnoString2On      = False
resw@vcRefAnnoString1        = "$VMG$" +" m/s"
resw@vcRefAnnoPerimSpaceF    = 0.5
resw@vcRefAnnoPerimOn        = False
resw@vcRefAnnoFontHeightF    = 0.015
resw@vcRefAnnoSide           = "Left"
resw@vcRefAnnoJust           = "TopRight"
resw@vcRefAnnoBackgroundColor = -1
plot    = new(4,graphic)
plotpsi = plot
plotw   = plot

type = "eps"
wks = gsn_open_wks(type,"/home/yinqc/SMJO/plot/Figure.4")
Label = (/"(a) March 22 2015 00:00 UTC","(b) March 23 2015 00:00 UTC", "(c) March 24 2015 00:00 UTC", "(d) March 25 2015 00:00 UTC"/)
do i = 0,3
    resolr@gsnCenterString = Label(i)
    plot(i) = gsn_csm_contour_map(wks, pw(i,:,:), resolr)
    plotpsi(i) = gsn_csm_contour(wks,  pv(i,:,:), respsi)
    plotw(i) = gsn_csm_vector(wks, uw(i,::12,::12), vw(i,::12,::12), resw)
    overlay(plot(i),plotpsi(i))
    overlay(plot(i),plotw(i))
end do

resP   = True
resP@gsnMaximize = True
resP@gsnFrame    = False
resP@gsnDraw     = True
resP@gsnPanelLabelBar = True
resP@gsnPanelTop    = 0.98
resP@gsnPanelBottom = 0.02
resP@pmLabelBarOrthogonalPosF = 0.0
resP@pmLabelBarWidthF  = 0.85
resP@lbTitlePosition   = "Right"
resP@lbTitleDirection  = "Across"
resP@lbTitleFontHeightF= 0.015
resP@lbBoxLinesOn      = False
resP@gsnPanelRowSpec   = True
gsn_panel(wks,plot,(/2,2/),resP)




