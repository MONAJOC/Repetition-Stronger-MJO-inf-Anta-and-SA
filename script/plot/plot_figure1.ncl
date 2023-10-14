;This script is uesed to plot olr(shading), wave vector(vectors) and psi(contour).

DATADIR = "/home/yinqc/SMJO/Data"

olrf  = addfile(DATADIR+"/OLR/olr.daily.anom.201503.nc","r")
wvecf = addfile(DATADIR+"/W_vector.nc","R")
psif  = addfile(DATADIR+"/CFSv2/psi200.anom.201503.nc","r")

;primary 
psi = calculate_daily_values(psif->psi,"avg",0,False)
olr = olrf->olr
wx  = wvecf->Wx
wy  = wvecf->Wy

time = olr&time
TIME = cd_calendar(time, 0)
day  = toint(TIME(:,2))

;calculate 3-day mean anomaly to plot.
olr_p = new((/4,dimsizes(olr&lat),dimsizes(olr&lon)/),float,olr@_FillValue)
copy_VarCoords(olr(0,:,:),olr_p(0,:,:))
copy_VarAtts(olr,olr_p)
psi_p = wx

do i = 0,3
    ihead = i*3 + 16
    iend  = i*3 + 18
    nt = ind(day.ge.ihead.and.day.le.iend)
    olr_p(i,:,:) = dim_avg_n_Wrap(olr(nt,:,:),0)
    psi_p(i,:,:) = dim_avg_n_Wrap(psi(nt,:,:),0)
end do
olr_p = lonFlip(olr_p)
psi_p = lonFlip(psi_p)
wx    = lonFlip(wx)
wy    = lonFlip(wy)
psi_p = psi_p*10^(-7)
psi_p = mask(psi_p, psi_p.lt.0.5.and.psi_p.gt.-0.5,False)
delete([/olr_p@long_name,olr_p@units,psi_p@long_name,psi_p@unit,wx@long_name,wx@unit,wy@long_name,wy@unit,psi_p@standard_name,psi_p@units/])
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
resolr@cnLabelBarEndStyle = "ExcludeOuterBoxes"
resolr@mpFillOn        = False
resolr@mpFillDrawOrder = "Postdraw"
;resolr@mpProjection= "Robinson"
resolr@mpMinLatF   = -85
resolr@mpMaxLatF   = 40
resolr@mpMinLonF   = 160
resolr@mpCenterLonF= 180.
resolr@mpMaxLonF   = 340
resolr@mpLandFillColor = "white"
resolr@mpGridAndLimbDrawOrder = "PostDraw"
resolr@mpGridLineColor        = "Gray60"
resolr@mpGridLineThicknessF   = 1.0
resolr@cnFillOn             = True
resolr@cnFillMode           = "CellFill"
resolr@cnFillPalette        = "BlueWhiteOrangeRed"
resolr@cnFillOpacityF       = 0.6
resolr@cnLinesOn            = False
resolr@cnLineLabelsOn       = False
resolr@cnLevelSelectionMode = "Manuallevels"
resolr@cnMinLevelValF       = -80.
resolr@cnMaxLevelValF       =  80.
resolr@cnLevelSpacingF      =  10.
resolr@lbLabelBarOn   = False


;contour plot
respsi = True
respsi@gsnFrame     = False
respsi@gsnDraw      = False
respsi@gsnMaximize  = False
respsi@gsnAddCyclic = True
respsi@gsnContourNegLineDashPattern = 1
respsi@cnMonoLineColor    = True
respsi@cnLinesOn          = True
respsi@cnLineLabelsOn     = False
respsi@cnLineColor        = "black"
respsi@cnLineThicknesses  = 3.0
respsi@cnInfoLabelPerimOn = False
respsi@cnInfoLabelOn      = True
respsi@cnLevelSelectionMode = "Manuallevels"
respsi@cnMinLevelValF       = -10.
respsi@cnMaxLevelValF       =  10.
respsi@cnLevelSpacingF      =  0.5

;vector plot.
resw = True
resw@gsnFrame     = False
resw@gsnDraw      = False
resw@gsnMaximize  = False
resw@gsnAddCyclic = False
resw@vcGlyphStyle            = "LineArrow"
resw@vcLineArrowThicknessF   = 1.
resw@vcMinDistanceF          = 0.025
resw@vcRefMagnitudeF         = 300.
resw@vcRefLengthF            = 0.05
resw@vcRefAnnoOrthogonalPosF = -1.0
resw@vcRefAnnoParallelPosF   = 1.18
resw@vcLineArrowColor        = "purple"
resw@vcRefAnnoString2On      = False
resw@vcRefAnnoString1        = "$VMG$" +" m~S~2~N~/s~S~2~N~"
resw@vcRefAnnoPerimSpaceF    = 0.5
resw@vcRefAnnoPerimOn        = False
resw@vcRefAnnoFontHeightF    = 0.015
resw@vcRefAnnoSide           = "Left"
resw@vcRefAnnoJust           = "TopRight"
resw@vcRefAnnoBackgroundColor = -1

respoly = True
respoly@gsLineColor      = "black"
respoly@gsLineThicknessF = 1.0


plot    = new(4,graphic)
plotpsi = plot
plotw   = plot

type = "eps"
type@wkOrientation = "portrait"
type@wkPaperWidthF   = 8.3
type@wkPaperHeightF  = 11.4

wks = gsn_open_wks(type,"/home/yinqc/SMJO/plot/Figure.1")
Label = (/"(a) 16-18 March 2015","(b) 19-21 March 2015", "(c) 22-24 March 2015", "(d) 25-27 March 2015"/)
do i = 0,3
    resolr@gsnCenterString = Label(i)
    plot(i) = gsn_csm_contour_map(wks, olr_p(i,:,:), resolr)
    plotpsi(i) = gsn_csm_contour(wks, psi_p(i,:,:), respsi)
    plotw(i) = gsn_csm_vector(wks, wx(i,::4,::4), wy(i,::4,::4), resw)
    overlay(plot(i),plotpsi(i))
    overlay(plot(i),plotw(i))
end do

polyy = (/-15.,-15.,-5.,-5.,-15.,-5./)
polyx = (/-180.,-150.,-150.,-180.,-180.,-180./)
dum = gsn_add_polyline(wks,plot(0), polyx, polyy, respoly) 

resP   = True
resP@gsnMaximize = True
resP@gsnFrame    = False
resP@gsnDraw     = True
resP@gsnPanelLabelBar = True
resP@gsnPanelTop    = 0.98
resP@gsnPanelBottom = 0.02
resP@pmLabelBarOrthogonalPosF = 0.0
resP@pmLabelBarWidthF  = 0.85
resP@lbTitleString     = "OLR anom (W/m2)"
resP@lbTitlePosition   = "Right"
resP@lbTitleDirection  = "Across"
resP@lbTitleFontHeightF= 0.015
resP@lbBoxLinesOn      = False
resP@lbLabelStride  = 2
resP@lbFillOpacityF = 0.6
resP@gsnPanelRowSpec   = True
gsn_panel(wks,plot,(/2,2/),resP)




