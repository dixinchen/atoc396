load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
  ; a = addfile("/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-17/var_dra_drc_nn.nc","r")
  ; b = addfile("/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-17/var_dra_drc_kernel.nc","r")
  a = addfile("/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-17/var_arctic_june_nn.nc","r")
  b = addfile("/Users/dchen/OneDrive - McGill University/myCourses/atoc_396/nov/nov-17/var_arctic_june_kernel.nc","r")

  tna = a->var_sfc_dra(:,:)
  tnc = a->var_sfc_drc(:,:)
  tnw = a->var_sfc_drw(:,:)
  tka = b->var_toa_dra(:,:)
  tkc = b->var_toa_drc(:,:)
  tkw = b->var_toa_drw(:,:)

  aa!0         = "lat"
  aa&lat       = a->latitude(:)
  aa!1         = "lon"
  aa&lon       = a->longitude(:)
  aa&lat@units="degree_north"
  aa&lon@units="degree_east"

  bb!0         = "lat"
  bb&lat       = b->latitude(:)
  bb!1         = "lon"
  bb&lon       = b->longitude(:)
  bb&lat@units="degree_north"
  bb&lon@units="degree_east"

  wks = gsn_open_wks("png","sfc_nn_drw")
  gsn_define_colormap(wks, "WhViBlGrYeOrReWh")   ; choose colormap

  res                      = True                ; plot mods desired
  ; res@tiMainString         = "Variance of dRa at TOA (NN)" ; plot title
  res@cnFillOn             = True               ; turns on the color
  res@cnFillPalette        = "BlAqGrYeOrReVi200"             ; set color map
  res@cnLinesOn            = False              ; turn off contour lines
  res@mpFillOn             = False            ; turn off gray filled land
  res@gsnPolar             = "NH"               ; specify the hemisphere
  res@mpMinLatF            = 70                ; specify min lat

  plot = gsn_csm_contour_map_polar(wks,aa(:,:),res)
end
