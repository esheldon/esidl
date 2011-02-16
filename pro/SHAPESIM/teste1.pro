PRO teste1, me1, fwhm, wme1, err, measerr


e1=.2
n=20
start=0.5
step=.1
fwhm=fltarr(n) + start
wt = fltarr(n)
me1 = fltarr(n)
err = fltarr(n)
measerr = fltarr(n)

nframes=5
galtype = 0
aratio=0.816497
theta=0.0
s2n = 100
FOR i=0, n-1 DO BEGIN

  fwhm[i] = start + i*step
  print,'HERE IS FWHM: ',fwhm[i]
  IF i NE 0 AND i MOD 10 EQ  0 THEN key = get_kbrd(20)
  testmom, nframes, galtype, aratio,theta,cat,s2n=s2n,fwhm=fwhm[i]
  

    ;; find stars
  w=where(cat.e1_ad NE -10.0, nw)

  IF nw GT 1 THEN BEGIN
    c=cat[w]
    me1[i] = median(c.e1_ad)
    measerr[i] = sqrt( (moment(c.e1_ad))[1] )
  
    ;; do weighted averages
    wt = 1.0/c.e1_aderr^2
    wtsum = total(wt)
    wme1[i] = total( wt*c.e1_ad )/wtsum
    err[i] = sqrt(1.0/wtsum)
  ENDIF
ENDFOR

plot, fwhm, me1, xtitle='fwhm (arcsec)',ytitle='e1 measured',$
         title = 'Shape bias',psym=4
oploterr, fwhm, me1, measerr
oplot,[0,10], [e1,e1]
key = get_kbrd(20)

plot, fwhm, wme1, xtitle='fwhm (arcsec)',ytitle='e1 measured',$
         title = 'Shape bias',psym=4
oploterr, fwhm, wme1, err
oplot,[0,10], [e1,e1]
key = get_kbrd(20)

plot, fwhm, wme1, xtitle='fwhm (arcsec)',ytitle='e1 measured',$
         title = 'Shape bias',psym=4
oploterr, fwhm, wme1, measerr
oplot,[0,10], [e1,e1]

return
END


