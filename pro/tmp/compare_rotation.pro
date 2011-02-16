PRO compare_rotation, cat, rot, shstruct, rotstruct

  simpctable
  !p.background=!white
  !p.color=!black

  sphoto_match, cat, rot, mc, mr
  combine_zlensum, cat[mc], 100., 20., 1000., 1.,sum,shstruct
  combine_zlensum, rot[mr], 100., 20., 1000., 1.,sum,rotstruct

  yrange=[-.001, .0025]
  erase & multiplot, [1,2]

  plot, [0], /nodata, xrange=[0,max(shstruct.meanr)], yrange=yrange, $
    ytitle='shear'
  oploterror, shstruct.meanr, shstruct.shear/shstruct.ssh, shstruct.shearerr/shstruct.ssh, $
             psym=5, color=!blue
  oploterror, rotstruct.meanr, rotstruct.shear/rotstruct.ssh, rotstruct.shearerr/rotstruct.ssh, $
             psym=4, color=!red
  oplot,[0,10000],[0,0]
  legend, ['original','new'], psym=[5,4], colors=[!blue,!red],/right
  multiplot
  plot, [0], /nodata, xrange=[0,max(shstruct.meanr)], yrange=yrange, $
    ytitle='orthoshear',xtitle='radius(kpc)'
  oploterror, shstruct.meanr, shstruct.ortho/shstruct.ssh, shstruct.orthoerr/shstruct.ssh, $
             psym=5, color=!blue
  oploterror, rotstruct.meanr, rotstruct.ortho/rotstruct.ssh, rotstruct.orthoerr/rotstruct.ssh, $
             psym=4, color=!red
  oplot,[0,10000],[0,0]
  multiplot,/reset

  key=get_kbrd(1)

  
  erase & multiplot, [1,2]
  error = sqrt( (shstruct.shearerr/shstruct.ssh)^2 + (rotstruct.shearerr/rotstruct.ssh)^2 )
  ploterror, shstruct.meanr, $
             (shstruct.shear/shstruct.ssh-rotstruct.shear/rotstruct.ssh), $
             error,  $
             psym=1, ytitle='orig-new'
  oplot,[0,10000],[0,0]
  multiplot
  error = sqrt( (shstruct.shearerr/shstruct.shear)^2 + (rotstruct.shearerr/rotstruct.shear)^2 )
  ploterror, shstruct.meanr, $
             (shstruct.shear/shstruct.ssh)/(rotstruct.shear/rotstruct.ssh), $
             error,  $
             psym=1, ytitle='orig/new',xtitle='radius(kpc)'
  oplot,[0,10000],[1,1]
  multiplot,/reset
return
end
