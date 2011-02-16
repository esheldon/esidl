FUNCTION denscont_sim_pq_limits, x

  return, [-sqrt(!rmax^2 - x^2), sqrt(!rmax^2 - x^2)]

END 


FUNCTION denscont_sim_sis2d, x, y

  ;; must have !cenx, !ceny, and !sigma_v defined

  ;; x,y in kpc
  radius = sqrt( (x-!cenx)^2 + (y-!ceny)^2 )

  ;; core of 1kpc

  term = 1./sqrt(1. + radius^2)
  return,  3.36e3*(!sigma_v/170.)^2*term ;Msolar/pc^2

END 

PRO denscont_sim_calcdensity, func, pqlim_func, $
                              Npts, Npts_circ, $
                              rmax, $
                              density_area, density_circ, density_contrast

  t1=systime(1)
  nrmax = n_elements(rmax)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; define variables
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; define points around unit circle
  xx = [arrscl(dindgen(Npts_circ), -1., 1.), $
        arrscl(dindgen(Npts_circ), -1., 1.)]
  yy = xx
  yy[0L:Npts_circ-1] = sqrt(1. - xx[0L:Npts_circ-1]^2)
  yy[Npts_circ:2L*Npts_circ-1] = -yy[0L:Npts_circ-1]

  ;; mean density in disk 
  density_area = dblarr(nrmax)
  ;; mean density around circle
  density_circ = density_area
  ;; density_contrast
  density_contrast = density_area


  FOR ir=0L, nrmax-1 DO BEGIN 

      defsysv, '!rmax', rmax[ir]
      ab_limits = [-!rmax, !rmax]

      ;; mean within radius !rmax: integrate
      area = !dpi*!rmax^2
      density_area[ir] = qgauss2d(func, ab_limits, pqlim_func,$
                                      Npts, /double)/area

      ;; mean around circle radius !rmax: generate x's and y's and
      ;; average
      density_circ[ir] = mean( denscont_sim_sis2d(!rmax*xx, !rmax*yy) )

      ;; density contrast
      density_contrast[ir] = density_area[ir] - density_circ[ir]

      print,'.',format='(a,$)'

  ENDFOR 
  print
  ptime,systime(1)-t1

END 

PRO density_cont_sim_plot, rmax, $
                  density_area_central, density_circ_central, density_contrast_central,$
                  density_area_neigh, density_circ_neigh, density_contrast_neigh,$
                  density_area, density_circ, density_contrast, loglog=loglog

  setup_mystuff

  ;; styles
  types = !csym.delta_cap+!csym.sigma_cap+'!D'+!csym.plus+'!N'
  types = types + ' !D'+['Total', 'Central', 'Neighbor']+'!N'
  types = [types, !csym.sigma_cap+' !DTotal!N     ']

  lines = [3, 2, 1, 0]

  colors = [!magenta, !red, !blue, !green]

  IF keyword_set(loglog) THEN BEGIN
      yrange=prange(density_contrast, density_contrast_central)
      yrange = [0.9*yrange[0], 10*yrange[1]]
      xlog=1 & ylog=1 
      aratio = 1.0
  ENDIF ELSE BEGIN 
      yrange=prange(density_contrast, density_contrast_central)
      yrange = [0.9*yrange[0], 1.1*yrange[1]]
      xlog=0 & ylog=0
      aratio = !gratio
  ENDELSE 

  xrange = [0.9*min(rmax), 1.1*max(rmax)]

  ytitle = '!S'+!csym.sigma_cap+'!R!A'+!csym.minus+'!N ('+!csym.ltequal+'R) '+$
    !csym.minus+' !S'+!csym.sigma_cap+'!R!A'+!csym.minus+$
    '!N (R) [M'+sunsymbol()+' pc!U'+!csym.minus+'2!N]'
  xtitle = 'R (kpc)'

  aplot, aratio,xlog=xlog,ylog=ylog,$
         rmax, density_contrast_central, $
         yrange=yrange,xrange=xrange,ytitle=ytitle,xtitle=xtitle,$
         line=lines[1], ystyle=1, xstyle=1
  oplot,rmax,density_contrast_central, color=!red, line=lines[1]

  oplot, rmax, density_contrast_neigh, color=!blue, line=lines[2]
  oplot, rmax, density_contrast, color=!magenta, line=lines[0]

  oplot, rmax, density_circ, color=!green, line=lines[0]
  
  add_labels, xtickv=[200,500,2000]

  IF !d.name EQ 'PS' THEN fac=0.7 ELSE fac=1.0
  legend, types, line=lines, color=colors, /right, charsize=!p.charsize*fac,$
          thick=replicate(!p.thick,4)

END 

PRO denscont_sim, central_sigmav, neighbor_sigmav, neighbor_pos, Npts, $
                  rmax, $
                  density_area_central, density_circ_central, density_contrast_central,$
                  density_area_neigh, density_circ_neigh, density_contrast_neigh,$
                  density_area, density_circ, density_contrast

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: denscont_sim, central_sigmav, neighbor_sigmav, neibhbor_dist, Npts , $'
      print,'     rmax, $'
      print,'     density_area_central, density_circ_central, density_contrast_central,$'
      print,'     density_area_neigh, density_circ_neigh, density_contrast_neigh,$'
      print,'     density_area, density_circ, density_contrast'
      print
      print,' neibhbor_dist in kpc, sigmav in km/s.  Note rmax is 2000kpc'
      return
  ENDIF 

  ;; Npts for integration 1500L works ok
  ;; improves some with 3000L

  print
  print,'Central object: sigma_v = '+ntostr(central_sigmav)+' km/s'
  print,'Position of neighbor: X = '+ntostr(neighbor_pos)+' kpc'
  print,'Neighbor object: sigma_v = '+ntostr(neighbor_sigmav)+' km/s'
  print,'Using Npts = '+ntostr(Npts)+' for integration'
  print

  ;; The mass model
  func = 'denscont_sim_sis2d'
  pqlim_func = 'denscont_sim_pq_limits'

  ;; number of points for mean around circle
  Npts_circ = Npts

  ;; number of radii
  nrmax = 100L
  
  ;; range of radii
  min_rmax = 25d                ;kpc
  max_rmax = 2000d

  rmax = arrscl( dindgen(nrmax), min_rmax, max_rmax)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; the central object
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  defsysv, '!sigma_v', double(central_sigmav)

  defsysv, '!cenx', 0d
  defsysv, '!ceny', 0d
  
  print,'Calculating density contrast from central object: sigma_v = '+$
        ntostr(!sigma_v)
  denscont_sim_calcdensity, func, pqlim_func, $
                            Npts, Npts_circ, $
                            rmax,$
                            density_area_central, density_circ_central, $
                            density_contrast_central

  ;; Now the neighbor
  defsysv, '!sigma_v', double(neighbor_sigmav)

  defsysv, '!cenx', -double(neighbor_pos)       ;kpc
  defsysv, '!ceny', 0d

  print,'Calculating density contrast from neighbor at center: ',!cenx,!ceny
  print,'Neighbor has sigma_v = '+ntostr(!sigma_v)
  denscont_sim_calcdensity, func, pqlim_func, $
                            Npts, Npts_circ, $
                            rmax,$
                            density_area_neigh, density_circ_neigh, $
                            density_contrast_neigh


  density_area = density_area_central + density_area_neigh
  density_circ = density_circ_central + density_circ_neigh

  density_contrast = density_area - density_circ

  density_cont_sim_plot, rmax, $
                  density_area_central, density_circ_central, density_contrast_central,$
                  density_area_neigh, density_circ_neigh, density_contrast_neigh,$
                  density_area, density_circ, density_contrast

  

END 
