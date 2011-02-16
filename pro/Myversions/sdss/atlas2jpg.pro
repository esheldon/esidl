;+
; NAME:
;   atlas2jpg
; PURPOSE:
;   Create a JPG of atlas
; CALLING SEQUENCE:
;   atlas_jpg, obj [, outdir=, cutout=, scales=, nonlinearity=, satvalue=, /all]
; INPUTS:
;   obj - [N] structure with RUN, CAMCOL, FIELD, ID, RERUN
; OPTIONAL INPUTS:
;   cutout - size of image to make (default 300)
;   scales - scaling values for JPG (default [20., 20., 20.]
;   nonlinearity - nonlinearity for JPG (default 3.)
;   satvalue - saturation for JPG (default 30.)
;   rebin - expand image into rebinned pixels (default 1)
; OPTIONAL KEYWORDS:
;   /all - if set, make an extra image with ALL field objects
; COMMENTS:
;   Makes output file in current directory of form:
;     atlas-RUN6-RERUN-CAMCOL-FIELD4-ID4.jpg
;   If /all is set makes an extra output file with all objects called:
;     atlas-RUN6-RERUN-CAMCOL-FIELD4-ID4-all.jpg
;   The rebin input allows you to make a nice big pixelized
;    image of a small object. It should be an integer >= 1
; REVISION HISTORY:
;   07-Jan-2009  Written by Mike Blanton, NYU
;   2009-01-28 Copied to altas2jpg, added more file flexibility and
;		sdssidl style naming
;-
;------------------------------------------------------------------------------
;;
function _atlas2jpg_rebin, image, rebin

	nx=(size(image,/dim))[0]
	ny=(size(image,/dim))[1]
	return, rebin(image, nx*rebin, ny*rebin, /sample)

end
;;
pro atlas2jpg, obj, outdir=outdir, files=files, afiles=afiles, all=all, cutout=cutout, scales=scales, nonlinearity=nonlinearity, satvalue=satvalue, rebin=rebin, sheldon=sheldon

	if(NOT keyword_set(scales)) then scales=[13.,13.,13.]
	if(NOT keyword_set(satvalue)) then satvalue=30.
	if(NOT keyword_set(nonlinearity)) then nonlinearity=3.
	if(NOT keyword_set(cutout)) then cutout=300
	if(NOT keyword_set(rebin)) then rebin=1

	delvarx, files, afiles

	for i=0L, n_elements(obj)-1L do begin
		afilename = sdss_objname(obj[i], prefix='atlas-', suffix='.jpg')
		ffilename = sdss_objname(obj[i], prefix='atlas-', suffix='-all.jpg')

		if n_elements(outdir) ne 0 then begin
			afilename = filepath(root=outdir, afilename)
			ffilename = filepath(root=outdir, ffilename)
		endif

		add_arrval, afilename, files
		add_arrval, ffilename, afiles
    
		gim=fpbin_to_frame(obj[i].run, obj[i].camcol, obj[i]. field, $
			obj[i].id, rerun=obj[i].rerun, filter=1, $
			/calibrate, /register, /literal, /center, $
			cutout=cutout, sheldon=sheldon)
		rim=fpbin_to_frame(obj[i].run, obj[i].camcol, obj[i]. field, $
			obj[i].id, rerun=obj[i].rerun, filter=2, $
			/calibrate, /register, /literal, /center, $
			cutout=cutout, sheldon=sheldon)
		iim=fpbin_to_frame(obj[i].run, obj[i].camcol, obj[i]. field, $
			obj[i].id, rerun=obj[i].rerun, filter=3, $
			/calibrate, /register, /literal, /center, $
			cutout=cutout, sheldon=sheldon)

		gim= _atlas2jpg_rebin(gim, rebin)
		rim= _atlas2jpg_rebin(rim, rebin)
		iim= _atlas2jpg_rebin(iim, rebin)
		djs_rgb_make, iim, rim, gim, name=afilename, $
			scales=scales, nonlinearity=nonlinearity, $
			satvalue=satvalue, quality=100.

		if(keyword_set(all)) then begin
			gim=fpbin_to_frame(obj[i].run, obj[i].camcol, obj[i]. field, $
				obj[i].id, rerun=obj[i].rerun, filter=1, $
				/calibrate, /register, /literal, /center, $
				cutout=cutout, /all, sheldon=sheldon)
			rim=fpbin_to_frame(obj[i].run, obj[i].camcol, obj[i]. field, $
				obj[i].id, rerun=obj[i].rerun, filter=2, $
				/calibrate, /register, /literal, /center, $
				cutout=cutout, /all, sheldon=sheldon)
			iim=fpbin_to_frame(obj[i].run, obj[i].camcol, obj[i]. field, $
				obj[i].id, rerun=obj[i].rerun, filter=3, $
				/calibrate, /register, /literal, /center, $
				cutout=cutout, /all, sheldon=sheldon)
			gim= _atlas2jpg_rebin(gim, rebin)
			rim= _atlas2jpg_rebin(rim, rebin)
			iim= _atlas2jpg_rebin(iim, rebin)
			djs_rgb_make, iim, rim, gim, name=ffilename, $
				scales=scales, nonlinearity=nonlinearity, $
				satvalue=satvalue, quality=100.
		endif

	endfor

end
