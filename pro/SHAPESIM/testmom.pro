pro testmom,numfr,galtype,aratio,theta,tot, galaxy,image,$
                smear=smear, start=start, galsize=galsize,fwhm=fwhm,s2n=s2n

;test the moment errors

if n_params() eq 0 then begin
   print,'-syntax  testmom, numfr, galtype, aratio, theta, tot, galaxy ,image, smear=smear, start=start'
   return
endif

if keyword_set(start) then begin
	print,'Starting at frame number: ',start
	stop = start+numfr-1
endif else begin
	start=1
	stop = numfr
endelse

;;;;;; Some image parameters
sky = 186.0
gain = 1.0

IF NOT keyword_set(galsize) THEN galsize = 100L

;;; input fwhm is in arcseconds
;;; use rzero in pixels 0.4arcsec/pixel
IF keyword_set(fwhm) THEN rzero = 0.42466091*fwhm/.4 ELSE rzero = 3.0
arrsize=4*galsize
image = fltarr(arrsize,arrsize)

print,'-----------------------------------------------'
;;;;;;;;;;Determine user input and make appropriate galaxy;;;;;;;;;;;;;
case 1 of

  (galtype eq 0): BEGIN 

    ;;;;;;;;;;;;GAUSSIAN GALAXIES;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    print, 'USING GAUSSIAN GALAXIES'
    sigma=rzero
    makegauss,galaxy,[galsize,galsize],sigma,counts=counts,aratio=aratio,$
      theta=theta
                  END 
  (galtype eq 1): BEGIN 
    
    ;;;;;;;;;;DEVAUCOULEURS GALAXIES;;;;;;;;;;;;;;;;;;;;;;;;;;;
    print, 'USING DEVAUCOULEURS GALAXIES'
    make_devauc,galaxy,[galsize,galsize],rzero,aratio=aratio,theta=theta
  END
  (galtype EQ 2): BEGIN
    print,'USING EXPONENTIAL DISKS'
    make_exp, edisk,[galsize,galsize],rzero, aratio=aratio, theta=theta
  END 
    
  (galtype eq 3): begin
      
    ;;;;;;;;;;Devaucouleurs and Exponential Disk;;;;;;;;;;;;;;
    print,'USING DEVAUC. AND EXP. DISK'
    make_devauc,galaxy,[galsize,galsize],rzero,aratio=aratio,theta=theta
    
    ;;;Use B/T of 0.4, assume light same in both pieces;;;;
    find_rd,0.4,rd
    make_exp, edisk,[galsize,galsize],rd*rzero, aratio=aratio, theta=theta

    ;;;;;;both;;;;;;;;;;;
    galaxy = galaxy + edisk
  END 
  ELSE: BEGIN 
    print,'galtype = ',galtype,' not valid'
    return
  END 
ENDCASE 

;;;;;;;;;;;;Smear up the image if requested;;;;;;;;;;;;;;;;;;;;;;;
if keyword_set(smear) then begin
	print,'Smearing Image'
        sigma = 2.0
        makegauss, psf, [3*sigma,3*sigma],sigma
	galaxy = convol(galaxy,psf,/edge_truncate)
endif
print,'-----------------------------------------------'

;;;;;;;;;;;;;;;;;build up statistics, do many frames;;;;;;;;;;;;;;

for i=start,stop do BEGIN
  
  print, 'rzero = ',rzero
  print, 'Aratio = ',aratio
  
  ;;; increment the s2n if not input
  IF NOT keyword_set(s2n) THEN s2n = float(i)
  print,'S/N: ',s2n

  counts = findflux(s2n, sky, gain=gain)
  gal = galaxy * counts
  for k = 0,3 do begin
    for l = 0,3 do begin
      image[k*galsize:(k+1)*galsize -1,l*galsize:(l+1)*galsize -1] = gal
    endfor
  endfor
  image = image+sky
  
  add_noise,image,gain=gain
  sdss_extract,image,cat
  if i eq start then begin
    tot=cat
  endif else begin
    concat_structs,tot,cat,temp
    tot=temp
  endelse
endfor


return
end
	

	

