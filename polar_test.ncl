load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
  a = addfile("/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-17/var_arctic_june_nn.nc","r")
  b = addfile("/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-17/var_arctic_june_kernel.nc","r")

; load var of nn data
  ta_nn = a->var_toa_dra(:,:)
  tc_nn = a->var_toa_drc(:,:)
  tw_nn = a->var_toa_drw(:,:)
  sa_nn = a->var_sfc_dra(:,:)
  sc_nn = a->var_sfc_drc(:,:)
  sw_nn = a->var_sfc_drw(:,:)
; load var of kernel data
  ta_k = b->var_toa_dra(:,:)
  tc_k = b->var_toa_drc(:,:)
  tw_k = b->var_toa_drw(:,:)
  sa_k = b->var_sfc_dra(:,:)
  sc_k = b->var_sfc_drc(:,:)
  sw_k = b->var_sfc_drw(:,:)

    printVarSummary(tw_k)
    printVarSummary(tw_nn)
; get total
  t_total_nn = ta_nn+tc_nn+tw_nn
  s_total_nn = sa_nn+sc_nn+sw_nn
  t_total_k = ta_k+tc_k+tw_k
  s_total_k = sa_k+sc_k+sw_k

  print(tw_nn(:,4))
; set coordinates for total data
  t_total_nn!0         = "lat"
  t_total_nn&lat       = a->latitude(:)
  t_total_nn!1         = "lon"
  t_total_nn&lon       = a->longitude(:)
  t_total_nn&lat@units="degree_north"
  t_total_nn&lon@units="degree_east"

  s_total_nn!0         = "lat"
  s_total_nn&lat       = a->latitude(:)
  s_total_nn!1         = "lon"
  s_total_nn&lon       = a->longitude(:)
  s_total_nn&lat@units="degree_north"
  s_total_nn&lon@units="degree_east"

  t_total_k !0         = "lat"
  t_total_k &lat       = b->latitude(:)
  t_total_k !1         = "lon"
  t_total_k &lon       = b->longitude(:)
  t_total_k &lat@units="degree_north"
  t_total_k &lon@units="degree_east"

  s_total_k !0         = "lat"
  s_total_k &lat       = b->latitude(:)
  s_total_k !1         = "lon"
  s_total_k &lon       = b->longitude(:)
  s_total_k &lat@units="degree_north"
  s_total_k &lon@units="degree_east"
; set coordinates for nn data at toa
  ta_nn!0         = "lat"
  ta_nn&lat       = a->latitude(:)
  ta_nn!1         = "lon"
  ta_nn&lon       = a->longitude(:)
  ta_nn&lat@units="degree_north"
  ta_nn&lon@units="degree_east"


  tc_nn!0         = "lat"
  tc_nn&lat       = a->latitude(:)
  tc_nn!1         = "lon"
  tc_nn&lon       = a->longitude(:)
  tc_nn&lat@units="degree_north"
  tc_nn&lon@units="degree_east"

  tw_nn!0         = "lat"
  tw_nn&lat       = a->latitude(:)
  tw_nn!1         = "lon"
  tw_nn&lon       = a->longitude(:)
  tw_nn&lat@units="degree_north"
  tw_nn&lon@units="degree_east"

; set coordinates for kernel data at toa
  ta_k!0         = "lat"
  ta_k&lat       = b->latitude(:)
  ta_k!1         = "lon"
  ta_k&lon       = b->longitude(:)
  ta_k&lat@units="degree_north"
  ta_k&lon@units="degree_east"

  tc_k!0         = "lat"
  tc_k&lat       = b->latitude(:)
  tc_k!1         = "lon"
  tc_k&lon       = b->longitude(:)
  tc_k&lat@units="degree_north"
  tc_k&lon@units="degree_east"

  tw_k!0         = "lat"
  tw_k&lat       = b->latitude(:)
  tw_k!1         = "lon"
  tw_k&lon       = b->longitude(:)
  tw_k&lat@units="degree_north"
  tw_k&lon@units="degree_east"

; set coordinates for nn data at sfc
  sa_nn!0         = "lat"
  sa_nn&lat       = a->latitude(:)
  sa_nn!1         = "lon"
  sa_nn&lon       = a->longitude(:)
  sa_nn&lat@units="degree_north"
  sa_nn&lon@units="degree_east"

  sc_nn!0         = "lat"
  sc_nn&lat       = a->latitude(:)
  sc_nn!1         = "lon"
  sc_nn&lon       = a->longitude(:)
  sc_nn&lat@units="degree_north"
  sc_nn&lon@units="degree_east"

  sw_nn!0         = "lat"
  sw_nn&lat       = a->latitude(:)
  sw_nn!1         = "lon"
  sw_nn&lon       = a->longitude(:)
  sw_nn&lat@units="degree_north"
  sw_nn&lon@units="degree_east"

; set coordinates for kernel data at toa
  sa_k!0         = "lat"
  sa_k&lat       = b->latitude(:)
  sa_k!1         = "lon"
  sa_k&lon       = b->longitude(:)
  sa_k&lat@units="degree_north"
  sa_k&lon@units="degree_east"

  sc_k!0         = "lat"
  sc_k&lat       = b->latitude(:)
  sc_k!1         = "lon"
  sc_k&lon       = b->longitude(:)
  sc_k&lat@units="degree_north"
  sc_k&lon@units="degree_east"

  sw_k!0         = "lat"
  sw_k&lat       = b->latitude(:)
  sw_k!1         = "lon"
  sw_k&lon       = b->longitude(:)
  sw_k&lat@units="degree_north"
  sw_k&lon@units="degree_east"

; ready to plot
wks = gsn_open_wks("png","panel")
cmap = read_colormap_file("gui_default")
plot = new(4,graphic)                          ; create a plot array

res                     = True
res@gsnDraw             = False           ; don't draw
res@gsnFrame            = False           ; don't advance frame
res@cnInfoLabelOn       = False           ; turn off cn info label
res@cnFillOn            = True            ; turn on color
res@cnFillPalette       = cmap(0:18,:)
res@lbLabelBarOn        = False
res@gsnPolar            = "NH"                   ;select norther hemisphere
res@mpMaxLatF           = 90          ; choose subregion
res@mpMinLatF           = 70
res@cnLinesOn           = False              ; turn off contour lines
res@mpFillOn            = False   ; turn off individual cb's1002cdx

res@gsnTickMarksOn      = False           ;turn off all tick marks
res@mpGridAndLimbOn     = False              ; turn off lat/lon lines
res@mpPerimOn           = False
res@txFontHeightF     = 0.07
res@cnLevelSelectionMode =  "ManualLevels"
res@cnMinLevelValF       = 0.
res@cnMaxLevelValF       = 50.
res@cnLevelSpacingF      = 10.

; res@gsnCenterString        = "Variance of dRa at TOA (NN)"
; plot(0)=gsn_csm_contour_map_polar(wks,ta_k(:,:),res)
; res@gsnCenterString        = "Variance of dRc at TOA (NN)"
; plot(1)=gsn_csm_contour_map_polar(wks,tc_k(:,:),res)


plot(0)=gsn_csm_contour_map_polar(wks,tw_nn(:,:),res)
plot(1)=gsn_csm_contour_map_polar(wks,tw_k(:,:),res)
plot(2)=gsn_csm_contour_map_polar(wks,sw_nn(:,:),res)
plot(3)=gsn_csm_contour_map_polar(wks,sw_k(:,:),res)
;************************************************
; create panel
;************************************************
resP                     = True                ; modify the panel plot
; resP@gsnPanelMainString = "A plot with a common label bar"

resP@gsnPanelLabelBar    = True                ; add common colorbar
resP@lbLabelFontHeightF  = 0.03              ; make labels smaller

gsn_panel(wks,plot,(/4,1/),resP)               ; now draw as one plot
end
