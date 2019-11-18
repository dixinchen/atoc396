load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
  a = addfile("/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-17/var_dra_drc_nn.nc","r")
  b = addfile("/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-17/var_dra_drc_kernel.nc","r")

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

  ; tc_nn!0         = "lat"
  ; tc_nn&lat       = a->latitude(:)
  ;
  ; tc_nn!1         = "lon"
  ; tc_nn&lon       = a->longitude(:)
  ;
  ; tc_nn&lat@units="degree_north"
  ; tc_nn&lon@units="degree_east"


  ; tc_k!0         = "lat"
  ; tc_k&lat       = b->latitude(:)
  ;
  ; tc_k!1         = "lon"
  ; tc_k&lon       = b->longitude(:)
  ;
  ; tc_k&lat@units="degree_north"
  ; tc_k&lon@units="degree_east"



  wks = gsn_open_wks("png","color3_bl")
  gsn_define_colormap(wks, "WhViBlGrYeOrReWh")   ; choose colormap

   res                      = True                ; plot mods desired

   res@tiMainString         = "Variance of dRa at TOA (NN)" ; plot title
   ; res@gsnCenterString      = "Paleo_Ice"                     ; center title

   res@cnFillOn             = True               ; turns on the color
  ;  res@cnLinesOn            = False              ; turn off contour lines
   res@cnFillPalette        = "BlAqGrYeOrRe"             ; set color map
   res@cnLinesOn            = False              ; turn off contour lines
   res@mpFillOn             = False            ; turn off gray filled land
  ; ; select which colors out of the colormap to use. These are the colors of
  ; ; the contours only. They do not effect the missing value, continental fill
  ; ; color etc. The -1 is no color or "transparent"
  ;  res@cnFillColors         = (/-1,4,9,19,35,81,75,91,99/)
  ;  res@cnFillDrawOrder       = "PreDraw"         ; make sure fill map on top
  ;
  ;  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  ;  res@cnMinLevelValF       = min(ta_nn)                 ; set min contour level
  ;  res@cnMaxLevelValF       = max(ta_nn)         ; set max contour level

   res@gsnPolar             = "NH"               ; specify the hemisphere
   res@mpMinLatF            = 70                ; specify min lat

plot = gsn_csm_contour_map_polar(wks,ta_nn(:,:),res)
end
