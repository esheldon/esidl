PRO s_g_contam, run, scontam, lcontam, maxf=maxf, read=read

  IF n_params() EQ 0 THEN BEGIN 
      print,'-Syntax: s_g_contam, run, scontam, lcontam, maxf=maxf'
      return
  ENDIF 

  IF NOT keyword_set(read) THEN read=0

  mincol = 1
  maxcol = 6
  ncolumn = maxcol-mincol+1

  red = 2
  runstr = ntostr(long(run))
  outdir = '/sdss4/data1/esheldon/CORRECTED/tmp/'
  indir = '/sdss3/usrdevel/philf/run'+runstr+'/' 

  lcontam = fltarr(ncolumn)
  scontam = lcontam

  IF read THEN BEGIN 
      tl = ['OBJC_TYPE','OBJC_FLAGS','FLAGS','PETROCOUNTS',$
            'PETRORAD','RA','DEC','E1','E2','R','MOMERR']

      IF run EQ 756 THEN BEGIN & start = 187 & nframes = 607 & END 
      IF run EQ 752 THEN BEGIN & start = 2   & nframes = 607 & END 

      ii = 0
      FOR col=mincol, maxcol DO BEGIN 

          tmp=0
          cstr = ntostr(col)
          typ = 'blah'+ntostr(long(systime(1)))+cstr
          file = indir + 'adat'+cstr+'c.fit'
          lensfile = outdir + 'run'+runstr+'_lensgal'+cstr+'.fit'
          srcfile = outdir + 'run'+runstr+'_srcgal'+cstr+'.fit'
          
          IF n_elements(maxf) EQ 0 THEN fits_info, file, /silent, n_ext=maxf

          read_photo_col, file, tmp, $
                      taglist=tl, $
                      start=start, nframes=nframes, $
                      maxf=maxf, $
                      struct_typ = typ

          ;; just select by red for now (as phil does)
      
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Source Galaxies
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          max_mag = 22.
          min_mag = 18.
          extract_lensgalaxies, tmp, red, src_ind, $
                            max_mag=max_mag, min_mag=min_mag,/plot

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Lens Galaxies
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          max_mag = 18.
          min_mag = 16.
          extract_lensgalaxies, tmp, red, lens_ind, $
                            max_mag=max_mag, min_mag=min_mag,/plot

          mwrfits, tmp[src_ind], srcfile, /create
          mwrfits, tmp[lens_ind], lensfile, /create
          
          ii = ii+1

      ENDFOR 
  ENDIF ELSE BEGIN 
      ii = 0
      FOR col=mincol, maxcol DO BEGIN 

          cstr = ntostr(col)
          lensfile = outdir + 'run'+runstr+'_lensgal'+cstr+'.fit'
          srcfile = outdir + 'run'+runstr+'_srcgal'+cstr+'.fit'

          lenses=0
          sources=0

          lenses = mrdfits(lensfile, 1,/silent) & nlens = n_elements(lenses)
          sources = mrdfits(srcfile, 1,/silent) & nsrc  = n_elements(sources)

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Lens Galaxies
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          max_mag = 18.
          min_mag = 16.
          plot,lenses.petrocounts(red), lenses.petrorad(red), psym=3, $
               ytitle='petrorad',xtitle='petrocounts',yrange=[0,10]
          print,format='($, "Petrorad Cutoff ")'
          read, pcut
          print,format='($, "Mag Cutoff ")'
          read, mcut
          w = where(lenses.petrorad[red] LE pcut AND $
                    lenses.petrocounts[red] LE mcut, nw)

          w2 = where(lenses.petrorad[red] LE pcut AND $
                     lenses.petrocounts[red] GE mcut, nw2)
          nw2 = nw2/2.

          ;; Extrapolate beyond the mag cut
          napprox = nw + nw2
          lcontam[ii] = napprox/nlens
          print,'Contamination ~ ',lcontam[ii]

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Source Galaxies
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          max_mag = 22.
          min_mag = 18.
          plot,sources.petrocounts(red), sources.petrorad(red), psym=3, $
               ytitle='petrorad',xtitle='petrocounts', yrange=[0,10]
          print,format='($, "Petrorad Cutoff ")'
          read, pcut
          print,format='($, "Mag Cutoff ")'
          read, mcut
          w = where(sources.petrorad[red] LE pcut AND $
                    sources.petrocounts[red] LE mcut, nw)

          w2 = where(sources.petrorad[red] LE pcut AND $
                     sources.petrocounts[red] GE mcut, nw2)
          nw2 = nw2/2.
          ;; Extrapolate beyond the mag cut
          napprox = nw + nw2
          scontam[ii] = napprox/nsrc
          print,'Contamination ~ ',scontam[ii]
          ii = ii+1
      ENDFOR 
  ENDELSE 

  print,'Mean lens contamination: ',total(lcontam)/ncolumn
  print,'Mean source contamination: ',total(scontam)/ncolumn

  return
END 
