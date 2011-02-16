pro sdss_pipe,filelist,dir,parfile,whyflag=whyflag,fwhm=fwhm
; NAME:
;       SDSS_PIPE
; PURPOSE:
;	Process all images described in file using sextractor 
;	
;
; CALLING SEQUENCE:
;       deep_pipe, filelist,dir,parfile
;
; INPUTS:
;       filename an array of strings with filenames (no path)
;
; OPTIONAL OUTPUT ARRAYS:
;
; INPUT KEYWORD PARAMETERS:
;
; PROCEDURE:
;	
;
; REVISION HISTORY:
;  Tim McKay	UM	6/69/98
;  David Johnston  UM      7/12/98
;      tuning parameters
;  Erin Scott Sheldon UM  rewrote for analysing sdss fields and simulations
;-
; On_error,2              ;Return to caller

if N_params() eq 0 then begin
        print,'Syntax - sdss_pipe, filelist,dir,parfile,'
	print,'whyflag=whyflag'
        return
endif
 
if n_elements(parfile) eq 0 then begin
	parfile='~/idl.lib/sdss.par'
endif

print,'parameters from ',parfile
numb=n_elements(filelist)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; setup the parameters for SExtractor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

sdss_setup,ps
ps.parameters_name=parfile
ps.starnnw_name=$
  '~/SExtractor/mysextractor2013/config/default.nnw'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Check on fwhm keyword
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

nfwhm = n_elements(fwhm)
IF keyword_set(fwhm) THEN BEGIN
  seeingfwhm=fwhm
  IF nfwhm NE 1 THEN BEGIN
    IF nfwhm NE numb THEN BEGIN
      print,'Number in fwhm not equal to number of input files.  Using first'
      seeingfwhm=fwhm[0]
    ENDIF 
  ENDIF 
ENDIF ELSE BEGIN
  print,'Using default seeing fwhm: ', strtrim(string(ps.seeing_fwhm),2)
  seeingfwhm = ps.seeing_fwhm
ENDELSE
nfwhm = n_elements(seeingfwhm)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Loop through the input files (images)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

for i=0,numb-1 do begin 	
  name=filelist(i)

  namearray=str_sep(name,'.')
  catfile=namearray(0)+"_sobj.fits"
  catfile=dir+catfile

  print
  print, "Will process frame:",dir+name
  print, "Will be writing objects to file:",catfile
  print

  IF nfwhm EQ 1 THEN ps.seeing_fwhm = seeingfwhm $
  ELSE ps.seeing_fwhm = seeingfwhm[i]
  print,'Using seeing fwhm:  ',ps.seeing_fwhm
  ps.catalog_name=catfile
  ps.checkimage_name = dir+'check.fits'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Call SExtractor
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  call_sex,ps,dir+name

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Add adaptive moments to the catalog
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  im=mrdfits(dir+name,0,hdr,/silent)

  ;; scale image if needed
  bzero=fxpar(hdr,'BZERO')
  bscale=fxpar(hdr,'BSCALE')
  if bscale ne 0 then begin
    im=im*bscale+bzero
  endif
  cat=mrdfits(catfile,1,hdr,/silent)
	
  ;put x,y in IDL notation (not SExtractor)
  cat.x_image=cat.x_image-1.0
  cat.y_image=cat.y_image-1.0

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;  Call adaptive moment wrapper
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ad_mom,cat,im,ixx,iyy,ixy,err,$
    numiter,wcenx,wceny,whyflag,rho4,$
    sky=cat.background

  nn=n_elements(cat)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; define a new strucure that will contain adaptive moment info
  ;; stuff beginning with 's' is simulated parameters that will be added
  ;; later
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  str=cat(0)
  strf2=create_struct(str,$
                      'sfield',-10,$
                      'stype',-10 ,$
                      'smag',-10.0,$
                      'saratio',-10.0,$
                      'sposangle',-10.0,$
                      'se1',-10.0,$
                      'se2',-10.0,$
                      'sr',-10.0,$
                      'e1',0.0,$
                      'e2',0.0,$
                      'e1_ad',0.0,$
                      'e2_ad',0.0,$
                      'e_ad',0.0,$
                      'x2_ad',0.0,$
                      'y2_ad',0.0,$
                      'xy_ad',0.0,$
                      'rho4',0.0, $
                      'r',-10.0,$
                      'uncert_ad',0.0,$
                      'psf_fwhm',0.0, $
                      's2n',0.0, $
                      'goodflag',0,$
                      'typeflag',-10)	
		
  str2=replicate(strf2,nn)
  copy_struct,cat,str2

  xx=cat.x2_image
  yy=cat.y2_image
  xy=cat.xy_image

  e1=replicate(-10.0,nn)         ;bad are set to -10.0
  e2=e1
  e=e1
  e1_ad=e1
  e2_ad=e1
  e_ad=e1
  r=e1

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Calculate ellipticities
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  t=xx+yy
  w=where(t gt .1)
  e1(w)=(xx(w)-yy(w))/t(w)
  e2(w)=2.0*xy(w)/t(w)
  e(w)=sqrt(e1(w)^2+e2(w)^2)
	
  t=ixx+iyy
  w=where(t gt .1, nw)
  IF (nw NE 0) THEN begin
    e1_ad(w)=(ixx(w)-iyy(w))/t(w)
    e2_ad(w)=2.0*ixy(w)/t(w)
    e_ad(w)=sqrt(e1_ad(w)^2+e2_ad(w)^2)
  ENDIF

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;  Calculate the polarization r  NOTE  This may not work for really 
  ;;;;;;;;;  Bad seeing                 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  good=where(e1_ad NE -10.0 AND rho4 GT .1,ngood)
  IF (ngood NE 0) THEN BEGIN
    stars = where(cat[good].class_star GT .9 AND $
                  sqrt(ixx[good]+iyy[good]) LT 2.0 AND $
                  sqrt(ixx[good]+iyy[good]) GT 1.2, nstars)
    IF nstars NE 0 THEN BEGIN
      stars=good[stars]
      psfrho4 = median(rho4[stars])
      psfsize = median( ixx[stars] + iyy[stars] )
      ;; Use gary's correction factor
      psfsize = psfsize*(4.0/psfrho4 - 1)
      sizes = (ixx[good] + iyy[good])*(4.0/rho4[good] -1)
      r[good] = psfsize/sizes
    ENDIF
  ENDIF

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Calculate signal to noise
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  f2n=cat.flux_aper
  s_two_n = f2n/sqrt(f2n+16.0*!pi*cat.background)
  s_two_n = s_two_n*sqrt(ps.gain)	;gain

  str2.e1 = e1
  str2.e2 = e2
  str2.x2_ad = ixx
  str2.y2_ad = iyy
  str2.xy_ad = ixy
  str2.rho4 = rho4
  str2.e1_ad = e1_ad
  str2.e2_ad = e2_ad
  str2.e_ad = e_ad
  str2.uncert_ad = err	
  str2.s2n = s_two_n
  str2.r = r

  mwrfits,str2,catfile,/create
endfor

return
end








