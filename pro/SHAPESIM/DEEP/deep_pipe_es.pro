pro deep_pipe_es,filelist,write_dir,read_dir,parfile,noback=noback,$
drobs=drobs,nosmear=nosmear
; NAME:
;       DEEP_PIPE
; PURPOSE:
;	Process all images described in file using sextractor 
;	Tuned for the deeprange project data
;
; CALLING SEQUENCE:
;       deep_pipe, filelist,write_dir,read_dir,parfile
;
; INPUTS:
;       filename an array of strings with filenames (no path)
;
; OPTIONAL OUTPUT ARRAYS:
;
; INPUT KEYWORD PARAMETERS:
;
; PROCEDURE:
;	the file 
;
; REVISION HISTORY:
;	Tim McKay	UM	6/69/98
;	David Johnston  UM      7/12/98
;	tuning parameters
;-
 On_error,2              ;Return to caller

if N_params() eq 0 then begin
        print,'Syntax - deep_pipe_es, filelist,write_dir,read_dir,parfile,'
	print,'noback=noback,drobs=drobs'
        return
endif
 
if n_elements(parfile) eq 0 then begin
;	parfile='/sdss/products/idltools/sdss/sim/deep.par'
	parfile='/usr/users/esheldon/idl.lib/deep_es.par'
endif

flagdir='/sdss/data1/lensing/deeprange/flags/'
halodir='/sdss/data2/lensing/deeprange/flags/'

print,'parameters from ',parfile

if n_elements(drobs) eq 0 then read_drobs,drobs,make=1
	
;deep_setup,ps
deep_setup_es,ps

numb=n_elements(filelist)
fwhm_list=fltarr(numb)
zero_pt=fltarr(numb)

assoc_drobs,filelist,drobs,drindex
wg=where(drindex ne -1,wgif)
wb=where(drindex eq -1,wbif)

if wgif gt 0 then gindex=drindex(wg)	;where info comes from drobs
if wbif gt 0 then bindex=drindex(wb)
	;where they must be set to default values

if wgif gt 0 then begin
	fwhm_list(wg)=drobs(gindex).fwhm
	zero_pt(wg)=31.0-drobs(gindex).zpt
endif
if wbif gt 0 then begin
	fwhm_list(wb)=1.3	
	zero_pt(wb)=31.0
	;nothing in drobs regarding this file -default
endif


for i=0,numb-1 do begin 	
	name=filelist(i)

	ps.parameters_name=parfile
	ps.starnnw_name=$
;	'/sdss/products/sextractor/sextractor1.2b10b/config/default.nnw'
	'/usr/users/esheldon/SExtractor/mysextractor2013/config/default.nnw'
	ps.gain='8.2'
	ps.pixel_scale='0.47'
	ps.seeing_fwhm=string(fwhm_list(i))
	ps.filter='N'
	ps.clean='Y'
	ps.deblend_mincont='.0001'
	ps.deblend_nthresh='32'
	ps.clean_param='1.5'
	ps.detect_thresh='2.2'
	ps.analysis_thresh='2.7'
	ps.satur_level='20000.0'
     	ps.back_size='20'
	ps.back_filtersize='3'
	ps.backphoto_type='GLOBAL'
	ps.memory_pixstack='1000000'
	ps.memory_bufsize='512'
	ps.mag_zeropoint=string(zero_pt(i))
	ps.verbose_type='normal'
	ps.phot_apertures='4'

	namearray=str_sep(name,'.')
	catfile=namearray(0)+"_sobj.fits"
	catfile=write_dir+catfile
	;name=read_dir+name

	if keyword_set(noback) then begin
		;make -background files
		ps.checkimage_type='-background'
		checkname=namearray(0)+'_noback.fits'
		checkname=write_dir+checkname
		ps.checkimage_name=checkname	
	endif	

	print
	print, "Will process corrected frame:",read_dir+name
	print, "Will be writing objects to file:",catfile
	print
	ps.catalog_name=catfile
;	deep_extract,ps,read_dir+name
	deep_extract_es,ps,read_dir+name		   
	    
        if keyword_set(noback) then begin
		check=mrdfits(checkname,0,hdr)
		check=round(check-10000.0)
		check=fix(check)	
		;pack it into integer
		writefits,checkname,check
	endif	

	;add adaptive moments
	im=mrdfits(read_dir+name,0,hdr)
	bzero=fxpar(hdr,'BZERO')
	bscale=fxpar(hdr,'BSCALE')
	if bscale ne 0 then begin
		im=im*bscale+bzero
	endif
	cat=mrdfits(catfile,1,hdr)
	
	;put x,y in IDL notation (not SExtractor)
	cat.x_image=cat.x_image-1.0
	cat.y_image=cat.y_image-1.0
	ad_mom,cat,im,ixx,iyy,ixy,err,$
	   numiter,wcenx,wceny,whyflag,$
	   sky=cat.background

	nn=n_elements(cat)
	;now add elliptcities and ad_moms
	str=cat(0)
	strf2=create_struct(str,$
	'e1',0.0,$
	'e2',0.0,$
	'x2_ad',0.0,$
	'y2_ad',0.0,$
	'xy_ad',0.0,$
	'e1_ad',0.0,$
	'e2_ad',0.0,$
	'e_ad',0.0,$
	'uncert_ad',0.0,$
	's2n',0.0,$
	'flag_sat',0,$
	'flag_halo',0,$
	'fwhm_post',0.0,$
	'ra',0.0,$
	'dec',0.0,$
	'p1',0.0,$
	'p2',0.0,$
	'out_focus',0,$
	'wild1',0.0,$
	'wild2',0.0)	
		
	str2=replicate(strf2,nn)
	copy_struct,cat,str2
	
	flagname=flagdir+strmid(name,0,4)+'_flag.fits'
	haloname=halodir+strmid(name,0,4)+'_halo.fits'
	
	if exist(flagname) and exist(haloname) then begin
		flagim=mrdfits(flagname,0,fhdr)
		haloim=mrdfits(haloname,0,hhdr)
		remove_sat,flagim,haloim,str2,10
	endif else print,'flag or halo file not found'

	xx=cat.x2_image
	yy=cat.y2_image
	xy=cat.xy_image

	t=xx+yy
	w=where(t gt .1)
	e1=replicate(-9.9,nn)	;bad are set to -9.9
	e2=e1
	e=e1
	e1_ad=e1
	e2_ad=e1
	e_ad=e1

	e1(w)=(xx(w)-yy(w))/t(w)
	e2(w)=2.0*xy(w)/t(w)
	e(w)=sqrt(e1(w)^2+e2(w)^2)
	
	t=ixx+iyy
	w=where(t gt .1)
	if (n_elements(w) eq 1) then begin
		if (w[0] ne -1) then begin
	  	  e1_ad(w)=(ixx(w)-iyy(w))/t(w)
	  	  e2_ad(w)=2.0*ixy(w)/t(w)
	  	  e_ad(w)=sqrt(e1_ad(w)^2+e2_ad(w)^2)
        	endif else begin
          	  e1_ad = -10.0
          	  e1_ad = -10.0
          	  e_ad  = -10.0
		endelse
	endif else begin
		e1_ad(w)=(ixx(w)-iyy(w))/t(w)
		e2_ad(w)=2.0*ixy(w)/t(w)
		e_ad(w)=sqrt(e1_ad(w)^2+e2_ad(w)^2)
	endelse

	f2n=cat.flux_aper
	s_two_n=f2n/sqrt(f2n+16.0*!pi*cat.background)
	s_two_n=s_two_n*sqrt(8.2)	;gain

	str2.e1=e1
	str2.e2=e2
	str2.x2_ad=ixx
	str2.y2_ad=iyy
	str2.xy_ad=ixy
	str2.e1_ad=e1_ad
	str2.e2_ad=e2_ad
	str2.e_ad=e_ad
	str2.uncert_ad=err	
	str2.s2n=s_two_n    			

	;this next part determins the stellar polarizations p1,p2
	;from the smear map , also marks the out of focus points
	;where stars are not found by the classifier
	;most of these areas are out of focus but not all of them

	if keyword_set(nosmear) eq 0 then begin
	
		smear_map,str2,smear,5,5,imsize=[2048,1988],/ad,$
		goodnum=4,class=.75,magmax=23.0,magmin=18.5,x2max=3	
	
		polarizations,smear,str2.x_image,str2.y_image,p1,p2,out_focus
	
		str2.p1=p1
		str2.p2=p2
		str2.out_focus=out_focus
	endif
	

	mwrfits,str2,catfile,/create
endfor

return
end







