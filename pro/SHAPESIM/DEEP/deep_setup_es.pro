pro deep_setup_es, param_struct
;+
; NAME:
;       REXTRACT_SETUP
; PURPOSE:
;	Set up a parameter structure for sextractor operation
;
; CALLING SEQUENCE:
;       rextract_setup
;
; INPUTS:
;       
; OUTPUTS:
;	param_struct: a structure containing all the parameters needed
;		      to run sextractor
;
; OPTIONAL OUTPUT ARRAYS:
;
; INPUT KEYWORD PARAMETERS:
; 
; PROCEDURE: This just reads in the formatted information from the CV file
;	
;
; REVISION HISTORY:
;	Tim McKay	UM	1/8/98
;	Tim McKay	UM	3/7/98  
;		Added check for environment variables 
;			EXTRACT_CONFIG, and EXTRACT_PAR
;-
 On_error,2              ;Return to caller

 if N_params() ne 1 then begin
        print,'Syntax - deep_setup_es, param_struct
        return
 endif
 
 config_dir=getenv('EXTRACT_CONFIG')
 if (config_dir eq "") then begin
;	config_dir='/sdss/products/sextractor/sextractor1.2b10b/config'
	config_dir='/usr/users/esheldon/SExtractor/mysextractor2013/config'
 endif
 par_dir=getenv('EXTRACT_PAR')
 if (par_dir eq "") then begin
;	par_dir='/sdss/products/idltools/rotse/rotse_idl/pipeline'
	par_dir='/usr/users/esheldon/idl.lib'
 endif

 param_struct = { $
	CATALOG_NAME: 'test.fts', $
	CATALOG_TYPE: 'FITS_1.0', $
	DETECT_TYPE: "CCD", $
	DETECT_IMAGE: 'SAME', $
	FLAG_IMAGE: 'flag.fits', $
	DETECT_MINAREA: '5', $
	DETECT_THRESH: '1.0', $
	ANALYSIS_THRESH: '1.2', $
	FILTER: 'Y', $
	FILTER_NAME: config_dir+"/gauss_1.5_3x3.conv", $
	DEBLEND_NTHRESH: '32', $
	DEBLEND_MINCONT: '0.00001', $
	CLEAN: 'N', $
	CLEAN_PARAM: '1.0', $
	BLANK: 'Y', $
	PHOT_APERTURES: '5', $
	PARAMETERS_NAME: par_dir+"/rotse.par", $
	SATUR_LEVEL: '16000.0', $
	MAG_ZEROPOINT: '30.0', $
	MAG_GAMMA: '4.0', $
	GAIN: '5.25', $
	PIXEL_SCALE: '14.4', $
	SEEING_FWHM: '20.0', $
	STARNNW_NAME: 'default.nnw',$
	BACK_SIZE: '32', $
	BACK_FILTERSIZE: '3', $
	BACKPHOTO_TYPE: 'GLOBAL', $
;	CHECKIMAGE_TYPE: 'MINIBACKGROUND', $
	CHECKIMAGE_TYPE: 'APERTURES', $
	CHECKIMAGE_NAME: '/sdss4/data1/esheldon/check.fits', $
	MEMORY_OBJSTACK: '30000', $
	MEMORY_PIXSTACK: '2000000', $
	MEMORY_BUFSIZE: '512', $
	SCAN_ISOAPRATIO: '0.6', $
	VERBOSE_TYPE: 'NORMAL'}

 return

 end



