function sdss_solarmags, band_shift=band_shift

    ; Only call blanton's code if band_shift is present and not zero
    common sdss_solarmags_common, solarmags0, old_band_shift, old_solarmags

    ; keep this for when no band_shift is sent
    if n_elements(solarmags0) eq 0 then solarmags0=k_solar_magnitudes()

    ; Just return the value unshifted
    if n_elements(band_shift) eq 0 then begin
        old_band_shift=0.0
        old_solarmags=solarmags0
        return, solarmags0
    endif

    ; We need to do some parameter checking now. First, make sure the 
    ; defaults are set
    if n_elements(old_band_shift) eq 0 then begin
        old_band_shift=band_shift
        old_solarmags=k_solar_magnitudes(band_shift=band_shift)
    endif

    ; Now see if we need to calculate it again (slow!).  We only need to
    ; calculate if band_shift has changed
    if band_shift eq old_band_shift then begin
        return, old_solarmags
    endif else begin
        solar_mags = k_solar_magnitudes(band_shift=band_shift)
        old_band_shift = band_shift
        old_solarmags = solar_mags
        return, solar_mags
    endelse
    
end
