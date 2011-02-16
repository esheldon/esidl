pro esrextract_setup, param_struct
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
;-
 On_error,2              ;Return to caller

 if N_params() ne 1 then begin
        print,'Syntax - esrextractor_setup, param_struct
        return
 endif
 
 param_struct = { $
	CATALOG_NAME: 'cat.fits', $
	CATALOG_TYPE: 'FITS_1.0', $
	DETECT_TYPE: "CCD", $
	DETECT_IMAGE: 'SAME', $
	FLAG_IMAGE: 'flag.fits', $
	FLAG_TYPE: 'OR',$
        DETECT_MINAREA: '5', $
	DETECT_THRESH: '1.0', $
	ANALYSIS_THRESH: '2.5', $
	FILTER: 'Y', $
	FILTER_NAME: $
   '/sdss/products/sextractor/sextractor1.2b10b/config/gauss_1.5_3x3.conv', $
	DEBLEND_NTHRESH: '32', $
	DEBLEND_MINCONT: '0.001', $
	CLEAN: 'N', $
	CLEAN_PARAM: '1.0', $
	BLANK: 'Y', $
	PHOT_APERTURES: '5', $
	PARAMETERS_NAME: $
	   '/sdss/products/idltools/sdss/sim/deep.par', $ 
	SATUR_LEVEL: '16000.0', $
	MAG_ZEROPOINT: '20.0', $
	MAG_GAMMA: '4.0', $
	GAIN: '18.0', $
	PIXEL_SCALE: '14.4', $
	SEEING_FWHM: '20.0', $
	STARNNW_NAME: 'default.nnw',$
	BACK_SIZE: '16', $
	BACK_FILTERSIZE: '3', $
	BACKPHOTO_TYPE: 'GLOBAL', $
	CHECKIMAGE_TYPE: 'APERTURES', $
	CHECKIMAGE_NAME: 'check.fits', $
	MEMORY_OBJSTACK: '30000', $
	MEMORY_PIXSTACK: '2000000', $
	MEMORY_BUFSIZE: '512', $
	SCAN_ISOAPRATIO: '0.6', $
	VERBOSE_TYPE: 'NORMAL'}

 return

 end


