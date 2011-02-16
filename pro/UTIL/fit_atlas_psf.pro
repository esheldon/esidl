PRO checkwin, wnum1, wnum2

  IF (!d.window NE wnum2) THEN BEGIN 
      WHILE !d.window NE -1 DO BEGIN 
          wdelete,!d.window
      ENDWHILE 
      wnum1 = !d.window+1
      wnum2 = !d.window+2
      window,wnum1,title='Atlas Image',xsize=500,ysize=400,xpos=0,ypos=1000
      window,wnum2,title='Fitted Gaussian',xsize=500,ysize=400,ypos=1000
  ENDIF

END 

PRO fit_atlas_psf, struct, clr

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: fit_atlas_psf, struct, clr'
      return
  ENDIF 


  scale = .4                    ;arcsec/pixel
  sky = 1000.                   ;sky in atlas images

  nn=n_elements(struct)

  oldwin = !d.window
  wnum1 = !d.window+1
  wnum2 = !d.window+2
  window,wnum1,title='Atlas Image',xsize=500,ysize=400,xpos=0,ypos=1000
  window,wnum2,title='Fitted Gaussian',xsize=500,ysize=400,xpos=1000,ypos=1000

  FOR i=0, nn-1 DO BEGIN 

      checkwin, wnum1, wnum2
      wset, wnum1

      run = struct[i].run
      camcol = struct[i].camcol
      rerun = struct[i].rerun

      print
      size=2.*sqrt(2.*alog(2))*sqrt( (struct[i].ixx[clr] + struct[i].iyy[clr])/2.)*scale
;      psfsize = struct[i].R[clr]*size
      print,'Adaptive Mom. Size: ',ntostr(size)
;      print,'Psf Size: ',ntostr(psfsize)

      fetch_dir,run, camcol, rerun, dir, atldir

      get_atlas, struct, i, dir=atldir, clr=clr, $
        imtot=imtot, row0=row0, col0=col0, /noprompt
      print_flags, struct, i, clr

      cenx = struct[i].colc[clr] - col0
      ceny = struct[i].rowc[clr] - row0

      checkwin, wnum2, wnum1
      wset,wnum2

      fitgauss2im, imtot, cenx, ceny, sky, scale=scale

      print,format='($, "Hit a key (q to quit) (p for previous): " )'
      key=get_kbrd(1)
      print

      IF (key EQ 'q') OR (key EQ 'Q') THEN BEGIN
          checkwin, wnum1, wnum2
          wdelete, wnum1
          wdelete, wnum2
         
          print
          return
      ENDIF 
      IF (key EQ 'p') OR (key EQ 'P') THEN BEGIN
          print
          i=i-2 > (-1)
      ENDIF 

  ENDFOR 
  checkwin, wnum1, wnum2
  wdelete, wnum1
  wdelete, wnum2


return
END 

      
