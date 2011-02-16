pro sdss_setup, param_struct
;+
; NAME:
;       SDSS_SETUP
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
; On_error,2              ;Return to caller

 if N_params() ne 1 then begin
        print,'Syntax - sdss_setup, param_struct'
        return
 endif
 
 config_dir=getenv('EXTRACT_CONFIG')
 IF (config_dir eq "") THEN $
   config_dir='~/SExtractor/mysextractor2013/config'

 par_dir=getenv('EXTRACT_PAR')
 IF (par_dir eq "") THEN par_dir='~/idl.lib/'

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;; Define SExtractor parameters to be input by hand
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 param_struct = { $
	ANALYSIS_THRESH: '2.7', $
	BACK_SIZE: '20', $
	BACK_FILTERSIZE: '3', $
	BACKPHOTO_TYPE: 'GLOBAL', $
;	BLANK: 'Y', $
	CATALOG_NAME: 'test.fits', $
	CATALOG_TYPE: 'FITS_1.0', $
	CHECKIMAGE_TYPE: 'APERTURES', $
	CHECKIMAGE_NAME: 'check.fits', $
	CLEAN: 'Y', $
	CLEAN_PARAM: '1.5', $
	DEBLEND_NTHRESH: '32', $
	DEBLEND_MINCONT: '0.0001', $
	DETECT_TYPE: "CCD", $
;	DETECT_IMAGE: 'SAME', $
	DETECT_MINAREA: '5', $
	DETECT_THRESH: '2.2', $
	FILTER: 'N', $
	FILTER_NAME: config_dir+"/gauss_1.5_3x3.conv", $
	FLAG_IMAGE: 'flag.fits', $
	GAIN: '5.25', $
	MAG_ZEROPOINT: '30.0', $
	MAG_GAMMA: '4.0', $
	MEMORY_OBJSTACK: '30000', $
	MEMORY_PIXSTACK: '1000000', $
	MEMORY_BUFSIZE: '3000', $
	PARAMETERS_NAME: par_dir+"sdss.par", $
	PHOT_APERTURES: '5', $
	PIXEL_SCALE: '.4', $
	SATUR_LEVEL: '20000.0', $
;	SCAN_ISOAPRATIO: '0.6', $
	SEEING_FWHM: '1.5', $
	STARNNW_NAME: 'default.nnw',$
	VERBOSE_TYPE: 'NORMAL'}

 return

 end



