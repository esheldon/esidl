PRO test_htm_find_neighbors, lcat, scat, depth, rmax, ind1, ind2,$
                             outfile=outfile,convert=convert, remove=remove,$
                             cheops=cheops

  ;; rmax in Mpc

  IF n_elements(lcat) EQ 0 THEN BEGIN 
      IF keyword_set(cheops) THEN BEGIN 
          lcat=mrdfits('/net/cheops1/data0/esheldon/testdata/stripe10_classgal_spec.fit',1)
      ENDIF ELSE BEGIN 
          lcat=mrdfits('/sdss7/data1/esheldon/stripe10_classgal_spec.fit',1)
      ENDELSE 
  ENDIF 
  IF n_elements(scat) EQ 0 THEN BEGIN 
      IF keyword_set(cheops) THEN BEGIN 
          scat=mrdfits('/net/cheops1/data0/esheldon/testdata/stripe10_srcgal_r_radec.fit',1)
      ENDIF ELSE BEGIN 
          scat=mrdfits('/sdss7/data1/esheldon/stripe10_srcgal_r_radec.fit',1)
      ENDELSE 
  ENDIF 

  
  w=where(lcat.z1d LT 0.4 AND lcat.z1d GT 0.025, nw)
  DL = angdist_lambda(lcat[w].z1d, h=1.0, omegamat=0.3) ;Mpc
  srad = double(rmax/DL)                ;radians

;  nw=n_elements(lcat)
;  w=lindgen(nw)
;  srad = replicate(600d/3600d*!dpi/180d, nw)

  ;htmlookupradec, scat.ra, scat.dec, depth, leafids

  htm_find_neighbors, depth, lcat[w].ra, lcat[w].dec, scat.ra, scat.dec, srad, $
   ind1, ind2,outfile=outfile, convert=convert, remove=remove, leafids=leafids

;  srad = srad*180d/!dpi
;  find_neighbors_radec, lcat[w].ra, lcat[w].dec, scat.ra, scat.dec, srad, ind1, ind2

END 
