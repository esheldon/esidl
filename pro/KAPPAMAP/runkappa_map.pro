PRO runkappa_map, lenses, w, $
                  slength=slength, $
                  rfac=rfac, $
                  gridsize=gridsize, $
                  stepfac=stepfac, $
                  yoverx=yoverx, $
                  allign=allign, $
                  sum=sum, $
                  scat=scat, $
                  fgal=fgal, $
                  write=write, $
                  surface=surface, rotate=rotate, $
                  noprompt=noprompt, verbose=verbose, $
                  doprint=doprint,$
                  donoise=donoise, $
                  outdir=outdir,$
                  nocheck=nocheck, $
                  _extra=extra


  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: runkappa_map, lenses, w, '
      print,'     slength=slength, '
      print,'     rfac=rfac, '
      print,'     gridsize=gridsize,' 
      print,'     stepfac=stepfac, '
      print,'     yoverx=yoverx, '
      print,'     allign=allign, '
      print,'     sum=sum, '
      print,'     scat=scat,' 
      print,'     write=write, '
      print,'     surface=surface, rotate=rotate, '
      print,'     noprompt=noprompt, verbose=verbose, '
      print,'     doprint=doprint,'
      print,'     donoise=donoise, '
      print,'     outdir=outdir,'
      print,'     _extra=extr'
      return
  ENDIF 

  COMMON seed,seed
  time = systime(1)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Some parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF NOT keyword_set(sum) THEN sum=0
  IF NOT keyword_set(doprint) THEN doprint=0
  IF NOT keyword_set(nocheck) THEN nocheck=0

  run1 = 752
  run2 = 756
  clr  = 2

  IF n_elements(gridsize) EQ 0 THEN gridsize = 1000.
  IF n_elements(w) EQ 0 THEN w=lindgen( n_elements(lenses) )
  nw = n_elements(w)
  print
  print,'Using ',ntostr(nw),' Lenses'
  print

  IF sum THEN nw = 1
  FOR i=0, nw-1 DO BEGIN 
      check = 1
      
      gsize = gridsize
 
      IF doprint THEN print,'Lens # ',ntostr(i+1),'/',ntostr(nw)
      IF sum THEN send = lenses[w] ELSE send = lenses[ w[i] ]
      WHILE check DO BEGIN 
          kappa_map, run1, run2, clr, send, $
                  slength=slength, $
                  rfac=rfac, $
                  gridsize=gsize, $
                  stepfac=stepfac, $
                  yoverx=yoverx, $
                  allign=allign, $
                  scat=scat, $
                  fgal=fgal, $
                  write=write, $
                  surface=surface, rotate=rotate, $
                  noprompt=noprompt, verbose=verbose,$
                  check=check, $
                  /abs, $
                  donoise=donoise, $
                  outdir=outdir, $
                  _extra=extra

          IF nocheck THEN check = 0
          IF check THEN BEGIN
              print
              print,'Shrinking gridsize'
              print
              gsize = gsize*.9
              IF gsize LE 200. THEN BEGIN
                  print,'Skipping This Lens'
                  check = 0
              ENDIF 
          ENDIF 
      ENDWHILE 
  ENDFOR 
 

  print
  ptime,systime(1)-time

  return
END 
