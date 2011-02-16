PRO compare, nframes, lupstr, adstr, ilup, iad

IF n_params() EQ 0 THEN BEGIN
  print,'-Syntax: compare, nframes, lupstr, adstr, ilup, iad'
  return
ENDIF


;; compare e1, e2, q, u by converting q,u to e1 and e2
pold = !p.multi
!p.multi = [0,1,2]

colors = ['u','g','r','i','z']
run=752
camcol = 3
clr=2

dir=run_dir(run)
tsobj_name, run,camcol, name
filec = dir+name

file = '~philf/run752/adat3c.fit'

make_corrected_tags, t

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Read in the catalogs if not input
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

IF n_elements(adstr) LT 10 THEN read_photo_col, file, adstr, $
                                nframes=nframes,taglist=t, str = 'bbllaahh'
IF n_elements(lupstr) LT 10 THEN read_photo_col, filec, lupstr, $
                                nframes=nframes, str='hheeyy'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Match up the structures if not input
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

IF n_elements(iad) LT 10 THEN matchphoto, adstr, lupstr, iad, ilup

ixx = adstr[iad].ixx[clr]
iyy = adstr[iad].iyy[clr]
ixy = adstr[iad].ixy[clr]

q=lupstr[ilup].q[clr]
u=lupstr[ilup].u[clr]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; calculate e1 and e2 from q and u
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

qu2e, q, u,  e1lup, e2lup

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; calculate raw e1 and e2 from adaptive moments
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

s = ixx + iyy
w=where(s NE 0.0 AND lupstr[ilup].fibercounts[clr] LT 23.0 $
       AND lupstr[ilup].fibercounts[clr] GT 18.0)

e1ad = (ixx[w] - iyy[w])/s[w]
e2ad = 2*ixy[w]/s[w]

wr = where(abs(e1ad) GT .01 AND abs(e2ad) GT .01 AND $
           abs(e1ad) LT  1.0 AND abs(e2ad) LT 1.0)
w=w[wr]

diffe1 = e1lup[w]-e1ad[wr]
diffe2 = e2lup[w]-e2ad[wr]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; calculate reletive difference
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

rele1 = diffe1/e1ad[wr]
rele2 = diffe2/e2ad[wr]
rele1med = median(rele1)
rele2med = median(rele2)
rele1sig = sqrt( (moment(rele1))[1] )
rele2sig = sqrt( (moment(rele2))[1] )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Bins for sigma vs magnitude plots
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

bins = [0., .5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
bins = bins+18.0

n=n_elements(bins)
mag = fltarr(n-1)
sigma1 = fltarr(n-1)
sigma2 = fltarr(n-1)

ytitle = 'Number'
xrange=[-1,1]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Calculate the sigma of distribution vs magnitude
; Plot the distributions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FOR i=0, n-2 DO BEGIN
  ww=where(lupstr[ ilup[w] ].fibercounts[clr] LT bins[i+1] AND $
           lupstr[ ilup[w] ].fibercounts[clr] GT bins[i] , nw)

  IF nw NE 0 THEN BEGIN
    mag[i] = median( lupstr[ ilup[w[ww]] ].fibercounts[clr]  )
    sigma1[i] = sqrt( ( moment( diffe1[ww] ) )[1] )
    sigma2[i] = sqrt( ( moment( diffe2[ww] ) )[1] )

;    title = colors[clr]+' Magnitude = '+ntostr(mag[i], 6)
;    xtitle='e1 lupton  -  e1 adaptive'
;    plothist, diffe1[ww], xhist, yhist, bin=.01,xrange=xrange,$
;              xtitle=xtitle,ytitle=ytitle, title=title
;    xyouts, .95*xrange[0],.95*max(yhist), 'sigma = '+ntostr(sigma1[i],6)
;    xtitle='e2 lupton  -  e2 adaptive'
;    plothist, diffe2[ww], xhist, yhist, bin=.01,xrange=xrange,$
;              xtitle=xtitle,ytitle=ytitle
;    xyouts, .95*xrange[0],.85*max(yhist), 'sigma = '+ntostr(sigma2[i],6)
;    key=get_kbrd(20)
    IF key EQ 'q' THEN BEGIN
      !p.multi=pold & return
    ENDIF 
  ENDIF 

ENDFOR 


continue = 1

WHILE continue DO BEGIN 
  continue = 0

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;Make histograms of differences
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  xrange=[-1,1]
  ytitle='Number'
  title = colors[clr] + ' Magnitude in [18.0,23.0]'
  xtitle='e1 lupton - e1 adaptive'
  plothist, diffe1, bin=.01,xrange=xrange,$
            ytitle=ytitle,xtitle=xtitle,title=title

  xtitle='e2 lupton - e2 adaptive'
  plothist, diffe2, bin=.01,xrange=xrange,ytitle=ytitle,xtitle=xtitle

  key = get_kbrd(20)
  IF key EQ 'q' THEN BEGIN
    !p.multi=pold & return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Make histograms of relative differences
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  xtitle='(e1 lupton - e1 adaptive)/e1 adaptive'
  plothist, rele1, bin=.01,xrange=xrange,ytitle=ytitle,xtitle=xtitle
  xyouts, -.75, 450, 'Median = '+ntostr(rele1med,6)
  ;xyouts,-.75, 400, 'Sigma = '+ntostr(rele1sig,6)

  xtitle='(e2 lupton - e2 adaptive)/e2 adaptive'
  plothist, rele2, bin=.01,xrange=xrange,ytitle=ytitle,xtitle=xtitle
  xyouts, -.75, 275, 'Median = '+ntostr(rele2med,6)
  ;xyouts,-.75, 250, 'Sigma = '+ntostr(rele2sig,6)

  key = get_kbrd(20)
  IF key EQ 'q' THEN BEGIN
    !p.multi=pold & return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;  Make plots of spread in differences vs magnitude
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  xtitle='Median Mag'
  ytitle = 'stdev( e1lup - e1ad )'
  yrange=[0,.5]
  plot, mag, sigma1, psym=7, xtitle=xtitle,ytitle=ytitle,yrange=yrange

  ytitle = 'stdev( e2lup - e2ad )'
  plot, mag, sigma2, psym=7, xtitle=xtitle,ytitle=ytitle,yrange=yrange

  key = get_kbrd(20)
  IF key EQ 'p' THEN continue = 1

ENDWHILE 
!p.multi = pold

return
END









