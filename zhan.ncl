load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
  a = addfile("/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-20/var_arctic_zhan.nc","r")


  atm= a->var_atm(:,:)
  sfc = a->var_sfc(:,:)
  total = a->var_total(:,:)

; set coordinates for total data
  atm!0         = "lat"
  atm&lat       = a->latitude(:)
  atm!1         = "lon"
  atm&lon       = a->longitude(:)
  atm&lat@units="degree_north"
  atm&lon@units="degree_east"

  sfc!0         = "lat"
  sfc&lat       = a->latitude(:)
  sfc!1         = "lon"
  sfc&lon       = a->longitude(:)
  sfc&lat@units="degree_north"
  sfc&lon@units="degree_east"

  total!0         = "lat"
  total&lat       = a->latitude(:)
  total!1         = "lon"
  total&lon       = a->longitude(:)
  total&lat@units="degree_north"
  total&lon@units="degree_east"


  wks  = gsn_open_wks ("png", "zhan" )              ; send graphics to PNG file
  plot = new (6, graphic)
  gsn_define_colormap(wks, "WhViBlGrYeOrReWh")   ; choose colormap

  res                      = True                ; plot mods desired
  ; res@tiMainString         = "Variance of dRa at TOA (NN)" ; plot title
  res@cnFillOn             = True               ; turns on the color
  res@cnFillPalette        = "BlAqGrYeOrRe"             ; set color map
  res@cnLinesOn            = False              ; turn off contour lines
  res@mpFillOn             = False            ; turn off gray filled land
  res@lbLabelBarOn        = False            ; turn off individual label bars
  res@gsnPolar             = "NH"               ; specify the hemisphere
  res@mpMinLatF            = 70
  res@cnLevelSelectionMode =  "ManualLevels"
  res@cnMinLevelValF       = 0.
  res@cnMaxLevelValF       = 200.
  res@cnLevelSpacingF      = 10.

  plot(0) = gsn_csm_contour_map_polar(wks, atm, res)
  plot(1) = gsn_csm_contour_map_polar(wks, sfc, res)
  plot(2) = gsn_csm_contour_map_polar(wks, total, res)

  pres                 = True
  pres@gsnPanelRowSpec = True                   ; tell panel what order to plot
  pres@gsnPanelLabelBar    = True       ; common label bar
  gsn_panel(wks,plot,(/2,2,2/),pres)

end
