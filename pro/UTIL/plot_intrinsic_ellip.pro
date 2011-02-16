PRO plot_intrinsic_ellip, str

  ;; plots of ellipticities etc.

  run=756
  rerun=1
  nframes=300
  start=300

  IF n_elements(str) EQ 0 THEN BEGIN 

      FOR camcol=1, 6 DO BEGIN 
          fetch_dir, run, camcol, rerun, dir, corrdir=corrdir
          read_tsobj, corrdir, tmp,$
                      front='adatc', /all, taglist=['petrocounts', $
                                                    'reddening',   $
                                                    'r',           $
                                                    'e1',          $
                                                    'e2']

          r = tmp.petrocounts[2] - tmp.reddening[2]
          gals=where(r LT 18.0 AND $
                     r GT 16.0, ngal)

          tmp = tmp[gals]
          
          IF camcol GT 1 THEN BEGIN 
              concat_structs, temporary(str), temporary(tmp), tmp2
              str=temporary(tmp2)
          ENDIF ELSE str=tmp

      ENDFOR 
  ENDIF 

  
  r = str.petrocounts[2] - str.reddening[2]

  clr=1
  ggals=where(str.e1[clr] NE 1.e10 AND str.R[clr] LT 0.8 AND $
              str.R[clr] GT 0.0, nggal)

  clr=2
  rgals=where(str.e1[clr] NE 1.e10 AND str.R[clr] LT 0.8 AND $
              str.R[clr] GT 0.0, nrgal)

  clr=3
  igals=where(str.e1[clr] NE 1.e10 AND str.R[clr] LT 0.8 AND $
              str.R[clr] GT 0.0, nigal)

  gec = sqrt( str[ggals].e1[1]^2 + str[ggals].e2[1]^2 )
  rec = sqrt( str[rgals].e1[2]^2 + str[rgals].e2[2]^2 )
  iec = sqrt( str[igals].e1[3]^2 + str[igals].e2[3]^2 )
  plothist, rec, rxhistc, ryhistc, bin=0.02, xrange=[0,1], /noplot
  plothist, gec, gxhistc, gyhistc, bin=0.02, xrange=[0,1], /noplot
  plothist, iec, ixhistc, iyhistc, bin=0.02, xrange=[0,1], /noplot

  gyhistc = gyhistc/total(gyhistc)
  ryhistc = ryhistc/total(ryhistc)
  iyhistc = iyhistc/total(iyhistc)

  simpctable
  aplot, !gratio, gxhistc, 2.*!pi*gxhistc*gyhistc, psym=10, xrange=[0,1], $
         xtitle='e', ytitle='2 '+!tsym.pi+' eP(e)'

  oplot, gxhistc, 2.*!pi*gxhistc*gyhistc, psym=10, color=!green
  oplot, rxhistc, 2.*!pi*rxhistc*ryhistc, psym=10, color=!red
  oplot, ixhistc, 2.*!pi*ixhistc*iyhistc, psym=10, color=!magenta

  legend,[!colorsp[1],!colorsp[2],!colorsp[3]],line=[0,0,0],$
         colors=[!green,!red,!magenta],$
         /bottom,/center,box=0,thick=replicate(!p.thick,3)

END 
