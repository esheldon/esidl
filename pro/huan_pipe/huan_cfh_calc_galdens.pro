PRO huan_cfh_calc_galdens, type, mag, fwhm, dens2d, star_fwhm, gooddens

  ;; This is doing things with teh fwhm only.  This pretty much sucks, don't
  ;; use it

  indir = '~/Huan/fwhm_only'
  
  fwhm_file23_7 = indir + 'matchcat12k.fwhm.0920.IAB.23.7'
  fwhm_file24_1 = indir + 'matchcat12k.fwhm.0920.IAB.24.1'

  ctfw_file23_7 = indir + 'matchcat12k.ctfw.0920.IAB.23.7'
  ctfw_file24_1 = indir + 'matchcat12k.ctfw.0920.IAB.24.1'


;  readcol, fwhm_file23_7, fwhm, num, dens, cumdens
  
;  w = where(fwhm GT 0.7, nw)
;  print,total(dens[w])

  IF type EQ 1 THEN file = ctfw_file23_7 ELSE file = ctfw_file24_1

  ;; fwhm are the galaxy sizes
  readcol, file, $
           mag_binnum, I_AB_low, I_AB_high, $
           fwhm_low, fwhm_high, $
           num, dens, cumdens, /silent

  rmd = rem_dup(mag_binnum)
  nmag = n_elements(rmd)
  mag = (I_AB_low[rmd] + I_AB_high[rmd])/2.0

  w=where(mag_binnum EQ 1, nfwhm)
  fwhm = (fwhm_low[w] + fwhm_high[w])/2.0
  

  dens2d = fltarr(nmag, nfwhm)

  FOR imag = 0L, nmag-1 DO BEGIN 
      w=where(mag_binnum EQ imag+1, nw)
      IF nw NE nfwhm THEN message,'What?'

      dens2d[imag, *] = dens[w]
  ENDFOR 

  xrange = [min(mag),max(mag)]
  yrange=[min(fwhm),max(fwhm)]

  tvim2, dens2d, xrange=xrange,yrange=yrange, $
         xtitle='I_AB', ytitle='FWHM'

  ;; Now figure out which galaxies are useful for lensing given a range
  ;; of PSf fwhm


  star_fwhm = [0.6, 0.7, 0.8, 0.9, 1.0]
;  star_fwhm = fwhm
  w=where(star_fwhm GE 0.6 AND star_fwhm LT 1.0)
  star_fwhm = star_fwhm[w]
  nstar_fwhm = n_elements(star_fwhm)

  gooddens = fltarr(nstar_fwhm)

  FOR i=0L, nstar_fwhm-1 DO BEGIN 

      r = star_fwhm[i]^2/fwhm^2

      w=where(r LT 0.8)

      gooddens[i] = total(dens2d[*, w])

  ENDFOR 

  colprint, star_fwhm, gooddens

;  plot, star_fwhm, gooddens, $
;        xtitle='PSF_FWHM', ytitle='density (#/arcmin!U2!N)'



END 
