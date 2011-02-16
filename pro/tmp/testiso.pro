PRO rotate_angle, theta

  w1=where(theta GT 0., n1)
  IF n1 NE 0 THEN theta[w1] = theta[w1] - 180.
  w2=where(theta LE 0., n2)
  IF n2 NE 0 THEN theta[w2] = 90. + theta[w2]

END 


PRO testiso, struct
  

  IF n_elements(struct) EQ 0 THEN BEGIN 

      fetch_dir, 756, 1, 1, dir, corrdir=corrdir
      read_corr, corrdir, tstruct, start=300,nf=100

      addsdsstag, tstruct, ['iso_a', 'iso_b', 'iso_phi'], struct
      delvarx, tstruct
  ENDIF 

  clr=2
  maxmag = 22.

  w=where(struct.e1[clr] NE 1.e10 AND struct.iso_a[clr] NE 0. AND $
          struct.petrocounts[clr] LT maxmag, nw)

  ;; iso_* are measured from different axis
  e1=struct[w].e1[clr]
  e2=struct[w].e2[clr]
  iso_b = struct[w].iso_b[clr]
  iso_a = struct[w].iso_a[clr]
  iso_phi = struct[w].iso_phi[clr]

  rotate_angle, iso_phi

;  rotate_e1e2,!pi/2.,e1,e2,e1_out,e2_out

  findabtheta, e1, e2, aratio, posangle
;  findabtheta, e1_out, e2_out, aratio, posangle
  posangle = posangle*180./!pi

  !p.multi=[0,1,2]

  simpctable
  w2=where(aratio NE -1000,nw2)
  w=w[w2]
  aplot, 1, aratio[w2], iso_b[w2]/iso_a[w2], psym=3, $
    xtitle='aratio (adaptive)', ytitle='aratio (iso_b/iso_a)'
  oplot,[-1000,1000], [-1000,1000], color=!red, thick=2
  aplot, 1, posangle[w2], iso_phi[w2], psym=3, $
    xtitle='pos angle (adaptive)', ytitle='pos angle (iso_phi)'
  oplot,[-1000,1000], [-1000,1000], color=!red,thick=2



  return
END 
