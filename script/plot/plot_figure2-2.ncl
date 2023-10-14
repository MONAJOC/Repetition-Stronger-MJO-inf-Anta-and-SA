;This plot the figure3.
;(a) gph200 average among March 22-24.
;(b) Ks.
;Information
;Date: 2023-10-05

DATADIR = "/home/yinqc/SMJO/Data"

psif = addfile(DATADIR+"/CFSv2/psi200.anom.201503.nc","r")
slpf = addfile(DATADIR+"/CFSv2/slp.anom.201503.nc","r")

;calculate mean uwind.
psi = calculate_daily_values(psif->psi, "avg", 0, 0)
slp = calculate_daily_values(slpf->slp, "avg", 0, 0)
psi := psi*10^(-7)
;Plot.
res = True
res@gsnFrame           = False
res@gsnDraw            = False
res@gsnMaximize        = True
res@gsnAddCyclic       = True
res@gsnLeftString      = ""
res@gsnRightString     = ""
res@gsnCenterStringFontHeightF = 0.025      
res@cnFillOn             = True
res@cnLinesOn            = False
res@cnLineLabelsOn       = False
res@cnLabelBarEndStyle = "ExcludeOuterBoxes"
res@lbOrientation      = "Vertical"
res@lbTitleString      = "[gpm]"
res@lbTitlePosition    = "Right"
res@lbBoxLinesOn       = False
res@cnFillPalette        = "nrl_sirkes"
res@cnLevelSelectionMode = "Manuallevels"
res@cnMinLevelValF       = -5.
res@cnMaxLevelValF       =  4.
res@cnLevelSpacingF      =  1.
res@tmXBMode = "Explicit"
x = ispan(1,31,1)
res@tmXBValues = x-1
res@tmXBLabels = x
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
Label = (/"(b) Stationary Wavenumber K~B~s~N~ 200 hPa"/)

plot(0) = gsn_csm_lat_time(wks, psi{}, res)



draw(plot(0))
delete(wks)