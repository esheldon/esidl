PRO test_rho4, struct

  IF n_elements(struct) EQ 0 THEN BEGIN

      taglist = ['IXX','IYY','IXY','RHO4',$
                 'R', 'MOMERR','E1','E2']

      run=756
      camcol=3
      rerun=1
      start=100
      nf=100
      fetch_dir, run, camcol, rerun, dir, corrdir=corrdir

      read_tsobj, corrdir, str, start=start, nf=nf, taglist=taglist, $
        front='adatc'

  ENDIF 

  clr=2
  w=where(str.rho4[clr] NE 0., nw)

  wgal = where(str[w].r[clr] LT 0.8 AND str[w].momerr[clr] LT 0.64)
  
  corr1 = 
