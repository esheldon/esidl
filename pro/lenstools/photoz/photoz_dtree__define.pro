function photoz_dtree::init
	return,1
end 


function photoz_dtree::match_struct, base_struct, n

	arr1 = replicate(-9999.0, 5)
	arr2 = fltarr(5)
	arr3 = replicate(9999.0, 5)
	newstruct={$
		matchid: -9999L, $
		modelflux:arr1, $
		modelflux_ivar: arr2, $
		modelmag: arr1, $
		modelmag_err: arr3 $
	}
	outst = create_struct(base_struct, newstruct)

	if n_elements(n) ne 0 then begin
		outst=replicate(outst, n)	
	endif
	return, outst
end

pro photoz_dtree::match2training, phot=phot, htmid_phot=htmid_phot, revphot=revphot
	sw=obj_new('sweeps')

	photfile=file_basename( sw->output_file('pzgal','m01',/gather) )
	if n_elements(phot) eq 0 then begin
		print,'reading photometric catalog'
		phot=sw->read_gather('pzgal','m01')
	endif

	tf='/home/users/esheldon/photoz/dtree/training/specz_ellip_dr6_slaq_deep_vagc7_11_2008.fit'

	matchf='/home/users/esheldon/photoz/dtree/training/specz_ellip_dr6_slaq_deep_vagc7_11_2008_match.fit'

	t=mrdfits(tf,1)

	outst=self->match_struct(t[0], n_elements(t))
	struct_assign, t, outst, /nozero

	print,'Matching to photometric catalog'
	angle = 2d/3600d*!dpi/180d
	htm_match, t.ra, double(t.dec), phot.ra, phot.dec, angle, $
		htmid2=htmid_phot, htmrev2=revphot, $
		mt, mphot, maxmatch=1

	help,mphot

	outst[mt].matchid = mphot
	outst[mt].modelflux = phot[mphot].modelflux_dered
	outst[mt].modelflux_ivar = phot[mphot].modelflux_dered_ivar

	outst[mt].modelmag = nmgy2mag(outst[mt].modelflux, $
		ivar=outst[mt].modelflux_ivar, $
		err=err)
	outst[mt].modelmag_err=err

	h=['']
	sxaddpar, h, 'photfile', photfile
	hdr={photfile:photfile}
	print,'Outputting matched catalog: ',matchf
	mwrfits, outst, matchf, h, /create
	;write_idlstruct, outst, matchf, hdr=hdr, /csv
end


pro photoz_dtree__define
  struct = {$
             photoz_dtree, $
             type:'' $
           }
end
