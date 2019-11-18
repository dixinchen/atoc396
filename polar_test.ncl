load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
  a = addfile("/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-17/var_dra_drc_global_nn.nc","r")
  b = addfile("/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-17/var_dra_drc_global_kernel.nc","r")

  ta_nn = a->var_toa_dra(:,:)
  tc_nn = a->var_toa_drc(:,:)
  sa_nn = a->var_sfc_dra(:,:)
  sc_nn = a->var_sfc_drc(:,:)
  ta_k = b->var_toa_dra(:,:)
  tc_k = b->var_toa_drc(:,:)
  sa_k = b->var_sfc_dra(:,:)
  sc_k = b->var_sfc_drc(:,:)

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




; plot = gsn_csm_contour_map_polar(wks,ta_nn(:,:),res)



wks = gsn_open_wks("png","panel")
cmap = read_colormap_file("gui_default")
plot = new(2,graphic)                          ; create a plot array

res                     = True
res@gsnDraw             = False           ; don't draw
res@gsnFrame            = False           ; don't advance frame
res@cnInfoLabelOn       = False           ; turn off cn info label
res@cnFillOn            = True            ; turn on color
res@cnFillPalette       = cmap(0:18,:)
res@lbLabelBarOn        = False
res@gsnPolar            = "NH"                   ;select norther hemisphere
res@mpMaxLatF           = 90          ; choose subregion
res@mpMinLatF           = 50
res@cnLinesOn           = False              ; turn off contour lines
res@mpFillOn            = False   ; turn off individual cb's
res@gsnTickMarksOn      = False           ;turn off all tick marks
res@mpGridAndLimbOn     = False              ; turn off lat/lon lines
res@mpPerimOn           = False
res@txFontHeightF     = 0.07
res@cnLevelSelectionMode =  "ManualLevels"
res@cnMinLevelValF       = 0.
res@cnMaxLevelValF       = 700.
res@cnLevelSpacingF      = 10.

; res@gsnCenterString        = "Variance of dRa at TOA (NN)"
; plot(0)=gsn_csm_contour_map_polar(wks,ta_k(:,:),res)
res@gsnCenterString        = "Variance of dRc at TOA (NN)"
; plot(1)=gsn_csm_contour_map_polar(wks,tc_k(:,:),res)


plot(0)=gsn_csm_contour_map_polar(wks,ta_nn(:,:),res)
plot(1)=gsn_csm_contour_map_polar(wks,tc_nn(:,:),res)
;************************************************
; create panel
;************************************************
resP                     = True                ; modify the panel plot
; resP@gsnPanelMainString = "A plot with a common label bar"

resP@gsnPanelLabelBar    = True                ; add common colorbar
resP@lbLabelFontHeightF  = 0.03              ; make labels smaller

gsn_panel(wks,plot,(/2,1/),resP)               ; now draw as one plot
end




; gsn_define_colormap(wks, "WhViBlGrYeOrReWh")   ; choose colormap  ; send graphics to PNG file
;
; plot_A = new(2,graphic)
; plot_B = new(2,graphic)
; ;*******************************************
; ; common resources
; ;*******************************************
; res                      = True
; res@cnFillOn             = True     ; turn on color fill
; res@cnFillPalette        = "gui_default"
; res@lbLabelBarOn         = False    ; no individual label bars
; res@cnLinesOn            = False    ; no contour lines
; res@mpFillOn             = False              ; turn off gray filled land
; res@gsnPolar             = "NH"               ; specify the hemisphere
; res@mpMinLatF            = 70
;
; res@tiMainString         = "Variance of dRa at TOA (NN)"
; plot_A(0)=gsn_csm_contour_map_polar(wks,ta_nn(:,:),res)
; res@tiMainString         = "Variance of dRc at TOA (NN)"
; plot_A(1)=gsn_csm_contour_map_polar(wks,tc_nn(:,:),res)
;
; res@tiMainString         = "Variance of dRa at TOA (kernel)"
; plot_B(0)=gsn_csm_contour_map_polar(wks,ta_k(:,:),res)
; res@tiMainString         = "Variance of dRc at TOA (kernel)"
; plot_B(1)=gsn_csm_contour_map_polar(wks,tc_k(:,:),res)
; ;*******************************************
; ;  create panel resources
; ;*******************************************
; pres1                  = True                ; mods desired
; pres1@gsnFrame         = False               ; save panel until both ready
; pres1@gsnPanelRight    = 0.5                 ; draw from left edge to center
; pres1@gsnPanelLabelBar = True                ; commom label bar
; pres1@lbLabelStride    = 2                   ; every other label
;
; pres2                  = True                ; mods desired
; pres2@gsnFrame         = False               ; save panel until both ready
; pres2@gsnPanelLeft     = 0.5		      ; draw from center to right edge
; pres2@gsnPanelLabelBar = True		      ; commom label bar
; pres2@lbLabelStride    = 4		      ; every 4th label
;
; gsn_panel(wks,plot_A,(/2,1/),pres1)          ; create first panel
; gsn_panel(wks,plot_B,(/2,1/),pres2)          ; create second panel
;
; frame(wks)                                   ; draw both on same page
; end
