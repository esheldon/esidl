PRO match_voronoi2lensout, stripe

  ;; this replaces lensum with new version containing the
  ;; voronoi density. Should copy old version to _old file

  stripestr = ntostr(long(stripe))
  get_tsgals, stripe, tsgals, /spec

  indir = '/sdss5/data0/lensout/stripe'+stripestr+'/'
  
  ;;;;;;;;;;;;;;;;;
  ;; match g'
  ;;;;;;;;;;;;;;;;;

  infile=indir+'zgal_gal_stripe'+stripestr+'_g_lensum_N1.fit'
  print,'infile: ',infile
  hdr0=headfits(infile)
  lensum=mrdfits(infile,1,hdr1)
  
  photo_match, tsgals.run, tsgals.rerun,tsgals.camcol,tsgals.field,tsgals.id,$
               lensum.run,lensum.rerun,lensum.camcol,lensum.field,lensum.id,$
               matchts,matchlensum

  help,matchlensum
  add_tag, lensum, 'voronoi_dens', -1d, newlensum
  newlensum[matchlensum].voronoi_dens = tsgals[matchts].voronoi_dens
  print,'Are you sure you want to replace file: ',infile,' (y/n)?'
  key=get_kbrd(1)
  IF strlowcase(key) EQ 'y' THEN BEGIN 
      FXHCLEAN, hdr1
      mwrfits2, newlensum, infile, hdr1, /create, hdr0=hdr0
  ENDIF 

  setzero,lensum, newlensum

  ;;;;;;;;;;;;;;;;;
  ;; match r'
  ;;;;;;;;;;;;;;;;;

  infile=indir+'zgal_gal_stripe'+stripestr+'_r_lensum_N1.fit'
  print,'infile: ',infile
  hdr0=headfits(infile)
  lensum=mrdfits(infile,1,hdr1)
  
  photo_match, tsgals.run, tsgals.rerun,tsgals.camcol,tsgals.field,tsgals.id,$
               lensum.run,lensum.rerun,lensum.camcol,lensum.field,lensum.id,$
               matchts,matchlensum

  help,matchlensum
  add_tag, lensum, 'voronoi_dens', -1d, newlensum
  newlensum[matchlensum].voronoi_dens = tsgals[matchts].voronoi_dens
  print,'Are you sure you want to replace file: ',infile,' (y/n)?'
  key=get_kbrd(1)
  IF strlowcase(key) EQ 'y' THEN BEGIN 
      FXHCLEAN, hdr1
      mwrfits2, newlensum, infile, hdr1, /create, hdr0=hdr0
  ENDIF 

  setzero,lensum, newlensum

  ;;;;;;;;;;;;;;;;;
  ;; match i'
  ;;;;;;;;;;;;;;;;;

  infile=indir+'zgal_gal_stripe'+stripestr+'_i_lensum_N1.fit'
  print,'infile: ',infile
  hdr0=headfits(infile)
  lensum=mrdfits(infile,1,hdr1)
  
  photo_match, tsgals.run, tsgals.rerun,tsgals.camcol,tsgals.field,tsgals.id,$
               lensum.run,lensum.rerun,lensum.camcol,lensum.field,lensum.id,$
               matchts,matchlensum

  help,matchlensum
  add_tag, lensum, 'voronoi_dens', -1d, newlensum
  newlensum[matchlensum].voronoi_dens = tsgals[matchts].voronoi_dens
  print,'Are you sure you want to replace file: ',infile,' (y/n)?'
  key=get_kbrd(1)
  IF strlowcase(key) EQ 'y' THEN BEGIN 
      FXHCLEAN, hdr1
      mwrfits2, newlensum, infile, hdr1, /create, hdr0=hdr0
  ENDIF 

  setzero,lensum,tsgals, newlensum


END 
