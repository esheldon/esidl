;+
; NAME:
;   sdss_catmatch
;
; PURPOSE:
;   Match RA,DEC positions to a set of external catalogs.
;
; CALLING SEQUENCE:
;   outdat = sdss_catmatch( [ra, dec, matchdist=, catalog= ] )
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   ra         - Right ascension(s) [deg]
;   dec        - Declination(s) [deg]
;   matchdist  - Matching distance for all catalogs; default to 3./3600 deg.
;   catalog    - External catalog names to match against; default to none;
;                setting with /CATALOG is equivalent to ['2MASS','FIRST','USNO']
;                The (case-insensitive) options are:
;                  FIRST
;                  2MASS
;                  Tycho
;                  USNO
;                This can be an array of catalog names, or a single
;                concatentated string such as '2MASS FIRST USNO'.
;
; OUTPUTS:
;   outdat     - Structure with one element per coordinate with match info;
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   If RA,DEC are not set, then return an empty structure.
;
; EXAMPLES:
;
; BUGS:
;   This routine is **not** well-suited for searching for a list
;   of objects spread out around the sky.  This is because it reads
;   all objects within a circular region on the sky that covers the list.
;   For example, if you have a list of objects at all RAs and/or DECs,
;   then the entire 2MASS catalog will be read -- this will run out
;   of memory.
;
; PROCEDURES CALLED:
;   copy_struct_inx
;   djs_angle_match()
;   first_coverage()
;   first_read()
;   tmass_read()
;   tycho_read()
;   usno_read()
;
; REVISION HISTORY:
;   17-Jul-2003  Written by D. Schlegel and D. Finkbeiner, Princeton.
;-
;------------------------------------------------------------------------------
function es_catmatch, ra, dec, matchdist=matchdist, $
 catalog=catalog1

   if (NOT keyword_set(catalog1)) then return, 0
   if (size(catalog1,/tname) NE 'STRING') then $
    catalog = ['2MASS','FIRST','USNO','RASS'] $
   else $
    catalog = strupcase(catalog1)

   q_tmass = total(strmatch(catalog,'*2MASS*')) NE 0
   q_first = total(strmatch(catalog,'*FIRST*')) NE 0
   q_tycho = total(strmatch(catalog,'*TYCHO*')) NE 0
   q_usno = total(strmatch(catalog,'*USNO*')) NE 0

   q_rass = total(strmatch(catalog,'*RASS*')) NE 0

   if (NOT keyword_set(matchdist)) then matchdist = 3./3600

   rass_dat = { $
		   rass_cat:       ' ',   $
           rass_ra:			0.d, $
           rass_dec:		0.d, $
		   rass_pos_err:	0,   $
		   rass_cps:		0.,  $
		   rass_cps_err:	0.,  $
		   rass_exp:        0,   $
		   rass_hr1:       0.,   $
		   rass_hr1_err:   0.,   $
		   rass_hr2:       0.,   $
		   rass_hr2_err:   0.,   $
           rass_matchdist:  -1.d, $
           rass_nmatch:  0 }

   first_dat = { $
           first_ra:      0.d, $
           first_dec:     0.d, $
           first_warning: ' ', $
           first_fint:    0., $
           first_fpeak:   0., $
           first_rms:     0., $
           first_maj:     0., $
           first_min:     0., $
           first_pa:      0., $
           first_fmaj:    0., $
           first_fmin:    0., $
           first_fpa:     0., $
           first_skyrms:  0., $
           first_matchdist:   -1., $
           first_nmatch:  0 }

   tmass_dat = { $
           tmass_ra:      0.d, $
           tmass_dec:     0.d, $
           tmass_err_maj: 0., $
           tmass_err_min: 0., $
           tmass_err_ang: 0., $
           tmass_j:       0.,$
           tmass_j_ivar:  0., $
           tmass_h:       0.,$
           tmass_h_ivar:  0., $
           tmass_k:       0., $
           tmass_k_ivar:  0., $
           tmass_ph_qual: '   ', $
           tmass_rd_flg:  0, $
           tmass_bl_flg:  0, $
           tmass_cc_flg:  '   ', $
           tmass_gal_contam: 0b, $
           tmass_mp_flg: 0b, $
           tmass_jdate: 0d, $
           tmass_dup_src: 0b, $
           tmass_matchdist:   -1., $
           tmass_nmatch:  0 }

   tyc_dat = { $
           tyc_ra:       0.d, $
           tyc_dec:      0.d, $
           tyc_pmra:     0. , $
           tyc_pmdec:    0. , $
           tyc_btmag:    0. , $
           tyc_vtmag:    0. , $
           tyc_matchdist:    -1. , $
           tyc_nmatch:   0 }

   usno_dat = { $
           usno_ra:      0.d, $
           usno_dec:     0.d, $
           usno_rmag:    0. , $
           usno_bmag:    0. , $
           usno_matchdist:   -1. , $
           usno_nmatch:  0 }

   ;----------
   ; Make a blank output structure

   if (q_first) then outdat = first_dat

   if (q_rass) then outdat = keyword_set(outdat) ? $
    create_struct(outdat,rass_dat) : rass_dat

   if (q_tmass) then outdat = keyword_set(outdat) ? $
    create_struct(outdat,tmass_dat) : tmass_dat
   if (q_tycho) then outdat = keyword_set(outdat) ? $
    create_struct(outdat,tyc_dat) : tyc_dat
   if (q_usno) then outdat = keyword_set(outdat) ? $
    create_struct(outdat,usno_dat) : usno_dat
   if (n_params() LT 2) then return, outdat
   outdat = replicate(outdat, n_elements(ra))

   ;----------
   ; Determine a central RA,DEC and search radius

   decmin = min(dec, max=decmax)
   if (decmin LT 0 AND decmax GT 0) then cosdec = 1. $
    else cosdec = cos( min( abs( [decmin,decmax] ) ) / !radeg)
   ramin = min(ra, max=ramax)
   rarange = ramax - ramin

   ; See if we can get a smaller RA range by shifting RAs by 180 deg
   ramin2 = min(ra + (ra LT 180)*360, max=ramax2)
   rarange2 = ramax2 - ramin2
   if (rarange2 LT rarange) then begin
      ramin = ramin2
      ramax = ramax2
      rarange = rarange2
   endif

   searchrad = 0.5*(decmax - decmin) > (0.5*rarange * cosdec)
   searchrad = searchrad + matchdist
   racen = 0.5 * (ramin + ramax)
   cirrange, racen
   deccen = 0.5 * (decmin + decmax)

   ;----------
   ; Match against FIRST

   if (keyword_set(q_first)) then begin
      thisdat = first_read(racen, deccen, searchrad)
      if (keyword_set(thisdat)) then begin
         nmatch = djs_angle_match(ra, dec, $
          thisdat.first_ra, thisdat.first_dec, dtheta=matchdist, $
          mcount=mcount, mindx=mindx, mdist=mdist)
         if (nmatch GT 0) then begin
            indx1 = where(mindx NE -1)
            indx2 = mindx[indx1]
            copy_struct_inx, thisdat, outdat, index_from=indx2, index_to=indx1
            outdat.first_nmatch = mcount
            outdat.first_matchdist = mdist
         endif
      endif

      ; Read the FIRST coverage maps
      outdat.first_skyrms = first_coverage(ra, dec)
   endif

   ;----------
   ; Match against RASS

   if (keyword_set(q_rass)) then begin
      thisdat = rass_read('all',/rassnames)
      if (keyword_set(thisdat)) then begin
         nmatch = djs_angle_match(ra, dec, $
          thisdat.rass_ra, thisdat.rass_dec, dtheta=matchdist, $
          mcount=mcount, mindx=mindx, mdist=mdist)
         if (nmatch GT 0) then begin
            indx1 = where(mindx NE -1)
            indx2 = mindx[indx1]
            copy_struct_inx, thisdat, outdat, index_from=indx2, index_to=indx1
            outdat.rass_nmatch = mcount
            outdat.rass_matchdist = mdist
         endif
      endif
   endif


   ;----------
   ; Match against 2-MASS

   if (keyword_set(q_tmass)) then begin
      thisdat = tmass_read(racen, deccen, searchrad)
      if (keyword_set(thisdat)) then begin
         nmatch = djs_angle_match(ra, dec, $
          thisdat.tmass_ra, thisdat.tmass_dec, dtheta=matchdist, $
          mcount=mcount, mindx=mindx, mdist=mdist)
         if (nmatch GT 0) then begin
            indx1 = where(mindx NE -1)
            indx2 = mindx[indx1]
            copy_struct_inx, thisdat, outdat, index_from=indx2, index_to=indx1
            outdat.tmass_nmatch = mcount
            outdat.tmass_matchdist = mdist
         endif
      endif
   endif

   ;----------
   ; Match against Tycho-2

   if (keyword_set(q_tycho)) then begin
      thisdat = tycho_read(racen=racen, deccen=deccen, radius=searchrad)
      if (keyword_set(thisdat)) then begin
         nmatch = djs_angle_match(ra, dec, $
          thisdat.ramdeg, thisdat.demdeg, dtheta=matchdist, $
          mcount=mcount, mindx=mindx, mdist=mdist)
         if (nmatch GT 0) then begin
            indx1 = where(mindx NE -1)
            indx2 = mindx[indx1]
            outdat[indx1].tyc_ra = thisdat[indx2].ramdeg
            outdat[indx1].tyc_dec = thisdat[indx2].demdeg
            outdat[indx1].tyc_pmra = thisdat[indx2].pmra
            outdat[indx1].tyc_pmdec = thisdat[indx2].pmde
            outdat.tyc_nmatch = mcount
            outdat.tyc_matchdist = mdist
         endif
      endif
   endif

   ;----------
   ; Match against USNO

   if (keyword_set(q_usno)) then begin
      thisdat = usno_read(racen, deccen, searchrad)
      if (keyword_set(thisdat)) then begin
         nmatch = djs_angle_match(ra, dec, $
          thisdat.ra, thisdat.dec, dtheta=matchdist, $
          mcount=mcount, mindx=mindx, mdist=mdist)
         if (nmatch GT 0) then begin
            indx1 = where(mindx NE -1)
            indx2 = mindx[indx1]
;            copy_struct_inx, thisdat, outdat, index_from=indx2, index_to=indx1
            outdat[indx1].usno_ra = thisdat[indx2].ra
            outdat[indx1].usno_dec = thisdat[indx2].dec
            outdat[indx1].usno_rmag = thisdat[indx2].rmag
            outdat[indx1].usno_bmag = thisdat[indx2].bmag
            outdat.usno_nmatch = mcount
            outdat.usno_matchdist = mdist
         endif
      endif
   endif

   return, outdat
end
;------------------------------------------------------------------------------
