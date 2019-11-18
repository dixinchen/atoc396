;*************************************************
; asr_4.ncl
;
; Concepts illustrated:
;   - Reading a multiple variables
;   - Associating grid coordinates with a variable
;     via the reserved attributes 'lat2d' and 'lon2d'
;   - Plotting all variables on a panel plot
;************************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;*******************************************
;read netCDF file
;*******************************************
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

  ;*******************************************
  ; Get dimensions
  ;*******************************************
  ;*******************************************
  ; Get human readable time
  ;*******************************************

  ;*******************************************
  ; Associate lat2d/lon2d with variables for georeferenced graphics
  ;*******************************************

  ta_nn!0         = "lat"
  ta_nn&lat       = a->latitude(:)

  ta_nn!1         = "lon"
  ta_nn&lon       = a->longitude(:)

  ta_nn&lat@units="degree_north"
  ta_nn&lon@units="degree_east"

  ; printVarSummary(ta_nn)





  ; Create plot(s)
  ; *******************************************

  ;************************************************
    wks = gsn_open_wks("png","color")               ; send graphics to PNG file

  ;***********************************************
  ; create array of named colors. Note that you
  ; should not include the background or foreground
  ; colors.
  ;***********************************************
    colors = (/"white","royal blue","light sky blue",\
  	     "powder blue","light sea green",\
                "pale green","wheat","brown","pink"/)

    res = True                                     ; plot mods desired

    res@tiMainString         = "Example of Using Named Colors" ; plot title
    res@gsnCenterString      = "Paleo_Ice"                     ; center title

    res@cnFillPalette        = colors
    res@cnFillOn             = True               ; turns on the color
    res@mpFillOn             = False              ; turns off continent gray
    res@cnLinesOn            = False              ; turn off contour lines

    ; res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
    ; res@cnMinLevelValF       = 1.                 ; set min contour level
    ; res@cnMaxLevelValF       = 8.                 ; set max contour level

    res@gsnPolar             = "NH"               ; specify the hemisphere
    res@mpMinLatF            = 70                 ; specify min lat

  ; note: since ice data is stored on a reduced grid with a data gap in
  ; latitude from -35 to +35 degrees, it is necessary to provide gsun
  ; with a sub-set of the data. Otherwise, an error will occur and the
  ; plot will not be correct.

    plot = gsn_csm_contour_map_polar(wks,ta_nn(:,:),res)
  end
