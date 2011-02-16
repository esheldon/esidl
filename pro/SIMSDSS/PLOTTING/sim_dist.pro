PRO sim_dist, cat, phcat=phcat

  if N_params() eq 0 then begin
    print,'Syntax: sim_dist, cat, phcat=phcat'
    return
  endif


sizerange=[0,10.0]
magrange=[14.0,26.0]
cattitle = 'Simulated (All)'
phtitle = 'SDSS (All)'
;cattitle='Binned'
;phtitle = 'Rebinned'


IF n_elements(phcat) NE 0 THEN BEGIN
  fpc = 1
  wp = where(phcat.x2_ad NE 0.0)
ENDIF ELSE fpc = 0
w=where(cat.goodflag EQ 1)

sizebin=0.3
magbin = 0.3
peak=1
; find stars
ws = where(cat[w].stype EQ 2)
ws = w[ws]
; find galaxies
wg = where(cat[w].stype EQ 1)
wg = w[wg]

galsym = 3
starsym = 4
gsymsize= .4
ssymsize = .4


continue = 0
WHILE (NOT continue) DO BEGIN

  plot, cat[wg].mag_best, sqrt(cat[wg].x2_ad + cat[wg].y2_ad), psym=galsym,$
    symsize=gsymsize, xtitle = 'mag_best', ytitle = 'Sqrt( ixx + iyy )', $
    title=cattitle, yrange=sizerange, xrange=magrange
  oplot, cat[ws].mag_best, sqrt(cat[ws].x2_ad + cat[ws].y2_ad), psym=starsym, $
    symsize=ssymsize
  legend,['Galaxies','Stars'],psym=[galsym,starsym]

  key=get_kbrd(1)
  IF (key EQ 'q') THEN return

  IF fpc THEN BEGIN 
    plot, phcat[wp].mag_best, sqrt(phcat[wp].x2_ad + phcat[wp].y2_ad), psym=3,$
      xtitle = 'mag_best', ytitle = 'Sqrt( ixx + iyy )', title=phtitle, $
      yrange=sizerange, xrange=magrange
    key=get_kbrd(1)
    IF (key EQ 'q') THEN return
  ENDIF
  IF (key NE 'p') THEN continue = 1
ENDWHILE

continue = 0
WHILE (NOT continue) DO BEGIN

  xtitle = 'Sqrt( ixx + iyy )'
  plothist,sqrt(cat[w].x2_ad + cat[w].y2_ad),xhist,yhist,bin=sizebin,$
    xtitle=xtitle,title=cattitle, xrange=sizerange, peak=peak
  key=get_kbrd(20)
  IF (key EQ 'q') THEN return

  IF fpc THEN BEGIN
    xtitle = 'Sqrt( ixx + iyy )'
    plothist,sqrt(phcat[wp].x2_ad + phcat[wp].y2_ad),xhist,yhist,bin=sizebin,$
      xtitle=xtitle,title=phtitle, xrange=sizerange, peak=peak
    key=get_kbrd(20)
    IF (key EQ 'q') THEN return
  ENDIF
  IF (key NE 'p') THEN continue = 1
ENDWHILE

continue = 0
WHILE (NOT continue) DO BEGIN
  xtitle = 'mag_best'
  plothist,cat[w].mag_best,xhist,yhist,bin=magbin,$
    xtitle=xtitle,title=cattitle,xrange=magrange, peak=peak
  key=get_kbrd(20)
  IF (key EQ 'q') THEN return

  IF fpc THEN BEGIN
    plothist,phcat[wp].mag_best,xhist,yhist,bin=magbin,$
      xtitle=xtitle,title=phtitle,xrange=magrange, peak=peak
    key=get_kbrd(20)
    IF (key EQ 'q') THEN return
  ENDIF 
  IF (key NE 'p') THEN continue = 1
ENDWHILE



IF n_elements(phcat) NE 0 THEN BEGIN
  fpc = 1
  wp = where(phcat.r NE -10.0)
ENDIF ELSE fpc = 0

continue=0
WHILE (NOT continue) DO BEGIN
  xtitle = 'mag_best'
  ytitle = 'R'

  plot,cat[wg].mag_best,cat[wg].sr,psym=galsym,symsize=gsymsize, $
    xrange=magrange,yrange=[0,2],ytitle=ytitle,xtitle=xtitle,title=cattitle
  oplot,cat[ws].mag_best,cat[ws].sr,psym=starsym,symsize=ssymsize
  legend,['Galaxies','Stars'],psym=[galsym,starsym]
  key=get_kbrd(20)
  IF (key EQ 'q') THEN return

  IF fpc THEN BEGIN
    plot,phcat[wp].mag_best,phcat[wp].r,psym=3,xrange=magrange,yrange=[0,2],$
         ytitle=ytitle,xtitle=xtitle,title=phtitle
    key=get_kbrd(20)
    IF (key EQ 'q') THEN return
  ENDIF 
  IF (key NE 'p') THEN continue = 1
ENDWHILE



return
END















