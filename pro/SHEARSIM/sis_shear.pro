PRO sis_shear, zlens, zsource, sigma, sdec, sra, cat, cen, shear, $
               incat=incat, $
               incenter=incenter, $
               radmin=radmin, $
               radmax=radmax, $
               density=density, $
               h=h,$
               omegamat=omegamat, $
               error=error, $
               plot=plot

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;       
; PURPOSE:
;	
;
; CALLING SEQUENCE:
;      
;                 
;
; INPUTS: 
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
;       
; OUTPUTS: 
;
; OPTIONAL OUTPUTS:
;
; CALLED ROUTINES:
; 
; PROCEDURE: 
;	
;	
;
; REVISION HISTORY:
;	
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF N_params() EQ 0 THEN BEGIN 
     print,'-Syntax: sis_shear, zlens, zsource, sigma, sdec, sra'
     print,'     [, cat, cen, '
     print,'     incat=incat, '
     print,'     incenter=incenter, '
     print,'     radmin=radmin,'
     print,'     radmax=radmax,'
     print,'     density=density, '
     print,'     h=h,'
     print,'     omegamat=omegamat,'
     print,'     plot=plot ]'
     print,'- sdec, sra in degrees'
     print,'Use doc_library,"sis_shear"  for more help.'  
     return
  ENDIF 

  IF n_elements(h) EQ 0 THEN h = .7
  IF n_elements(omegamat) EQ 0 THEN omegamat=1.0
  IF NOT keyword_set(plot) THEN plot=0
  IF NOT keyword_set(error) THEN error=0

  IF n_elements(density) EQ 0 THEN density = 7200. ;gal per square degree

  ;; if set, give objects intrinsic shape
  ;; average intrinsic shape
  e_int = double(.32)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Find angular diameter distances
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  DS = angdist(zsource, h=h, omegamat=omegamat)
  DLS = angdist(zsource, zlens, h=h, omegamat=omegamat)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Factor in front of shear equation
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  shfact = DLS/DS*(sigma/220.)^2*1.4  ;; 1.4 arcsec

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; If not input, make the catalog
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(incat) EQ 0 THEN BEGIN 
      sra = double(sra)
      sdec = double(sdec)
      ramin  = double(0.)
      decmin = double(0.)

      ;; Assuming sx, sy in degrees
      area = sdec*sra
      ngal = long( density*area )

      ;; Default ra, dec position
      ramax = ramin + sra
      decmax = decmin + sdec
      cen = [ mean([decmax,decmin]), mean([ramax, ramin]) ]

      typ='tmp'+ntostr(long(systime(1)))
      default = double(0.)
      str = create_struct(name=typ, $
                          'ra', default, $
                          'dec', default, $
                          'e1', 0.,$
                          'e2', 0., $
                          'uncert',0. )

      cat = replicate(str, ngal)
      shear = replicate(str, ngal)
      ra = arrscl( randomu(seed, ngal), ramin, ramax, arrmin=0., arrmax=1. )
      dec = arrscl( randomu(seed, ngal), decmin, decmax, arrmin=0., arrmax=1. )
      cat.ra = ra
      cat.dec = dec
      shear.ra = ra
      shear.dec = dec
  ENDIF ELSE BEGIN 
      ngal = n_elements(incat)
      cat = incat
      shear = incat
      dec = cat.dec
      ra  = cat.ra
      IF n_elements(incenter) EQ 0 THEN BEGIN
          cen=[median(dec), median(ra)]
      ENDIF ELSE cen = incenter
  ENDELSE 
      
  x = (dec - cen[0])*3600.     ; Relative positions in arcseconds
  y = (ra  - cen[1])*3600. 
  
  R = sqrt( x^2 + y^2 )

  IF n_elements(radmax) NE 0 THEN $
    w=where( R LT radmax, nw ) $
  ELSE BEGIN
      w=lindgen(ngal)
      nw=ngal
  ENDELSE 

  IF n_elements(radmin) NE 0 THEN BEGIN
      w2=where( R[w] GT radmin, nw )
      w=w[w2]
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Shear the catalog
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  tn = 5*nw
  continue=1
  WHILE continue DO BEGIN 
      continue=0
      e=abs( e_int*randomn(seed,tn)+.55 )
      theta=arrscl(randomu(seed,tn),-!pi/2.,!pi/2.,$
                   arrmin=0., arrmax=1.)
      e1 = e*cos(2*theta)
      e2 = e*sin(2*theta)
      ;; can remove this later, as is seen in 
      ;;real dist
      ;;g=where( e1^2 + e2^2 LE 1., ng)
      ng=tn
      g=lindgen(tn)
      IF ng LT nw THEN continue=1 ELSE BEGIN 
          cat[w].e1 = e1[g(0:nw-1)]
          cat[w].e2 = e2[g(0:nw-1)]
          cat[w].uncert = .22
      ENDELSE 
  ENDWHILE

  fac = shfact/R[w]^3

  shear.uncert = 1.
  shear[w].e1 = fac*( y[w]^2 - x[w]^2 )/2. < .8
  shear[w].e2 = fac*( -2*x[w]*y[w] )/2. < .8

  cat[w].e1 = cat[w].e1 + shear[w].e1*2. ; for sis, kappa=shear=e/2.
  cat[w].e2 = cat[w].e2 + shear[w].e2*2.

  IF error THEN BEGIN 
      print,'Adding Measurement Error'
      file='/sdss3/usrdevel/esheldon/idl.lib/SHEARSIM/e_uncert_hist.fit'
      uhist=mrdfits(file,1,/silent)
      genrand, uhist.hist, uhist.x, nw, uncert, /double
      ff = ( lindgen(nw) MOD 2 )*2 - 1 ;-1 and 1
      s=sort(randomu(seed,nw))
      s2 = sort(randomu(seed,nw))

      cat[w].e1 = cat[w].e1 + ff[s]*uncert[s]
      cat[w].e2 = cat[w].e2 + ff[s2]*uncert[s]
      cat[w].uncert = uncert[s]
  ENDIF 

  IF plot THEN plot, cat[w].dec, cat[w].ra, psym=3

  return 
END 
