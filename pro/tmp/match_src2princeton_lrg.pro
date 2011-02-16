PRO match_src2princeton_lrg_convert

  print
  print,'Converting file: '+file+' --> '+outfile
  print,'stopping: are you sure you want to overwrite?'
stop
  file = '~/Rachel/LRGcatalog.tbl'
  outfile = '~/Rachel/LRGcatalog_fix.dat'

  nlines = numlines(file)
  
  openr, lun, file, /get_lun
  openw, olun, outfile, /get_lun

  ra = dblarr(nlines)
  dec = dblarr(nlines)

  e1 = fltarr(nlines)
  e2 = e1
  e1e1var = e1
  e1e2var = e1
  e2e2var = e1

  c0 = e1
  c1 = e1
  c2 = e1
  c3 = e1
  c4 = e1

  e0err = e1
  c1err = e1
  c2err = e1
  c3err = e1
  c4err = e1

  f1 = lonarr(nlines)
  f2 = f1
  f3 = f1

  r_i = e1
  photoz_z = e1
  photoz_z_lrg = e1

  line = ''

  tt = systime(1)
  FOR i=0L, nlines-1 DO BEGIN 

      ;;print,'-------------------------------------------------------'
      readf, lun, line
      ;;print,line

      tmp = strsplit(line, /extract)

      hex2dec,tmp[17],f1, /quiet
      tmp[17] = ntostr(f1)

      hex2dec,tmp[18],f2, /quiet
      tmp[18] = ntostr(f2)

      hex2dec,tmp[19],f3, /quiet
      tmp[19] = ntostr(f3)

      outline = strjoin(tmp, ' ')
      printf,olun,outline

  ENDFOR 
  ptime,systime(1)-tt
  free_lun, lun
  free_lun, olun


END 

PRO match_src2princeton_lrg, struct, leafids, scat, match1, match2

  file = '~/Rachel/LRGcatalog_fix.dat'

  IF NOT fexist(file) THEN match_src2princeton_lrg_convert
  outfile = '~/Rachel/LRGscat.fit'

  IF n_elements(struct) EQ 0 THEN BEGIN 
      s=create_struct('ra',0d,'dec',0d,$
                      'e1',0.,'e2',0.,$
                      'e1e1var',0.,'e1e2var',0.,'e2e2var',0.,$
                      'counts_model',fltarr(5), 'counts_modelerr', fltarr(5), $
                      'flags',0L,'flags2',0L,'flags3',0L,$
                      'r_i',0.,$
                      'photoz_z',0.,$
                      'photoz_z_lrg',0.)

      read_struct, file, s, struct

      depth = 10
      htmLookupRadec, struct.ra, struct.dec, depth, leafids
      
  ENDIF 

  IF n_elements(scat) EQ 0 THEN BEGIN 
      ;; no stripe 9 in their cat
      stripes = [10,11,12,13,14,15,$
                 27,28,29,30,31,32,33,34,35,36,37,$
                 76,82,86]

      columns = ['clambda','ceta','e1_recorr','e2_recorr',$
                 'e1e1err','e1e2err','e2e2err','photoz_z','photoz_zerr']
      get_scat, stripes, [1,2,3], scat, $
        /nzcuts, /pixelMaskFlags, /rlrgMask, $
        columns=columns
  ENDIF 

  IF n_elements(match1) EQ 0 AND n_elements(match2) EQ 0 THEN BEGIN 
      csurvey2eq, scat.clambda, scat.ceta, ra, dec
      
      w=where(abs(scat.e1_recorr) lt 2 and abs(scat.e2_recorr) lt 2)

      tol = 3d/3600d
      tt = systime(1)
      spherematch, $
        struct.ra, struct.dec, ra[w], dec[w], $
        tol, match1, match2, dis12

      match2 = w[match2]

      ptime, systime(1)-tt
  ENDIF 



  nmatch = n_elements(match1)
  s = create_struct('leafid', 0L, $
                    scat[0], $
                    'photoz_z_lrg', 0.0)
  tcat = replicate(s, nmatch)

  tcat.leafid = leafIDs[match1]

  tcat.clambda = scat[match2].clambda
  tcat.ceta = scat[match2].ceta

  tcat.e1_recorr = scat[match2].e1_recorr
  tcat.e2_recorr = scat[match2].e2_recorr

  tcat.e1e1err = scat[match2].e1e1err
  tcat.e1e2err = scat[match2].e1e2err
  tcat.e2e2err = scat[match2].e2e2err

  tcat.photoz_z = scat[match2].photoz_z
  tcat.photoz_zerr = scat[match2].photoz_zerr

  tcat.photoz_z_lrg = struct[match1].photoz_z_lrg

  sl = sort(tcat.leafid)
  tcat = tcat[sl]

  minid = min(tcat.leafid, max=maxid)
  leaf_hist = histogram(tcat.leafid, min=minid, max=maxid, rev=rev)

  hdr = ['END']
  sxAddPar, hdr0, 'depth', depth, ' HTM depth used to calculate leafids'
  print
  print,'Writing file: ',outfile
  mwrfits2, tcat, outfile, /create, /destroy, hdr0=hdr0
  print,'Appending the reverse indices'
  mwrfits, rev, outfile

END 
