pro test_distort_invert
	; this file has a SIP header
	file='~/data/astrometry/newheader.fits'

	hdr=headfits(file)	

	extast, hdr, astr

	;hprint, hdr

	naxis1 = sxpar(hdr, 'naxis1')
	naxis2 = sxpar(hdr, 'naxis2')
	new_astr = hogg_ab2apbp(astr, [naxis1,naxis2])

end
