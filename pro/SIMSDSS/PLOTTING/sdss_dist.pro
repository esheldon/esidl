PRO sdss_dist, str, nfields, galsize, galmag, starmag, factor=factor

IF (n_params() EQ 0) THEN BEGIN
  print,'-Syntax: sdss_dist,str, nfields, galsize, galmag, starmag, factor=factor'
  return
ENDIF

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; extraction routines throw out stuff.  This is a test
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

IF (NOT keyword_set(factor) ) THEN factor = 8.0
nfields=float(nfields)

w=where(str.x2_ad NE 0.0)

sizerange = [0,10.0]
magrange = [14.0,26.0]
starmagrange = [16.0,26.0]
galmagrange = magrange

sizebin=0.3
magbin = 0.3
galmagbin = 0.3
starmagbin = 0.2

starbreak = .9
galbreak = .2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Everything first
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

plot,str[w].mag_best,sqrt(str[w].x2_ad + str[w].y2_ad),psym=3,$
     xtitle='mag_best',ytitle='Sqrt( ixx + iyy )',title='Everything',$
     yrange=sizerange, xrange=magrange
key=get_kbrd(20)
IF (key EQ 'q') THEN return

xtitle = 'Sqrt( ixx + iyy )'
plothist,sqrt(str[w].x2_ad + str[w].y2_ad),xhist,yhist,bin=sizebin,$
  xtitle=xtitle,title='Everything',xrange=sizerange
key=get_kbrd(20)
IF (key EQ 'q') THEN return

xtitle='mag_best'
plothist,str[w].mag_best,xhist,yhist,bin=magbin,$
         xtitle=xtitle,title='Everything', xrange=magrange
key=get_kbrd(20)
IF (key EQ 'q') THEN return


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Look at the galaxies
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


wg = where(str[w].class_star LT galbreak)
wg = w[wg]

plot,str[wg].mag_best,sqrt(str[wg].x2_ad + str[wg].y2_ad),psym=3,$
     xtitle='mag_best',ytitle='Sqrt( ixx + iyy )',title='Galaxies',$
     yrange=sizerange, xrange=galmagrange
key=get_kbrd(20)
IF (key EQ 'q') THEN return

                      ;; Sizes
xtitle = 'Sqrt( ixx + iyy )'
plothist,sqrt(str[wg].x2_ad + str[wg].y2_ad),xhist,yhist,bin=sizebin,$
  xtitle=xtitle,title='Galaxies',xrange=sizerange

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; IMPORTANT:  Multiplying by about 8 to get actual number of galaxies
;;; seen out to 24
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

peak=factor*max(yhist)/nfields
plothist,sqrt(str[wg].x2_ad + str[wg].y2_ad),xhist,yhist,$
         bin=sizebin,peak=peak,title='Galaxies Per Field',xtitle=xtitle,$
         xrange=sizerange

;;; make nicer return values
low = where(yhist EQ 0.0,nlow)
IF (nlow NE 0) THEN BEGIN 
  good = where(yhist NE 0.0)
  mingood = min(yhist[good])
  yhist[low] = mingood
ENDIF

n=n_elements(yhist)
galsize = create_struct('yhist',fltarr(n),'xhist',fltarr(n))
galsize.xhist=xhist
galsize.yhist=yhist
key=get_kbrd(20)
IF (key EQ 'q') THEN return

                      ;; Magnitudes
xtitle='mag_best'
plothist,str[wg].mag_best,xhist,yhist,bin=galmagbin,$
         xtitle=xtitle,title='Galaxies',xrange=galmagrange

peak=factor*max(yhist)/nfields
plothist,str[wg].mag_best,xhist,yhist,bin=galmagbin,peak=peak,$
  title='Galaxies Per Field',xtitle=xtitle,xrange=galmagrange

;;; make nicer return values
low = where(yhist EQ 0.0,nlow)
IF (nlow NE 0) THEN BEGIN 
  good = where(yhist NE 0.0)
  mingood = min(yhist[good])
  yhist[low] = mingood
ENDIF

n=n_elements(yhist)
galmag = create_struct('maghist',fltarr(n),'mag',fltarr(n))
galmag.mag=xhist
galmag.maghist=yhist
key=get_kbrd(20)
IF (key EQ 'q') THEN return

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Look at the stars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ws = where(str[w].class_star GT starbreak AND $
           sqrt(str[w].x2_ad + str[w].y2_ad) LT 2.0 )
ws = w[ws]


xtitle='mag_best'
ytitle='sqrt(ixx + iyy)'
plot,str[ws].mag_best,sqrt(str[ws].x2_ad + str[ws].y2_ad),psym=3,$
  xtitle=xtitle,ytitle=ytitle, title='Stars'
key=get_kbrd(20)
IF (key EQ 'q') THEN return

xtitle='mag_best'
plothist,str[ws].mag_best,xhist,yhist,bin=starmagbin,xtitle=xtitle,$
  title='Stars',xrange=starmagrange

peak=factor*max(yhist)/nfields
plothist,str[ws].mag_best,xhist,yhist,bin=starmagbin,peak=peak,$
  title='Stars Per Field',xtitle=xtitle,xrange=starmagrange


;;; make nicer return values
low = where(yhist EQ 0.0,nlow)
IF (nlow NE 0) THEN BEGIN 
  good = where(yhist NE 0.0)
  mingood = min(yhist[good])
  yhist[low] = mingood
ENDIF

n=n_elements(yhist)
starmag = create_struct('maghist',fltarr(n),'mag',fltarr(n))
starmag.mag=xhist
starmag.maghist=yhist


return
END










