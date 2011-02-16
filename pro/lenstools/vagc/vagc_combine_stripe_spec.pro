;; Don't need separate lss files, can just stuff all into database
;; and select what we want

PRO vagc_combine_stripe_spec, stripe, outdir=outdir, lrg=lrg, southern=southern, lss=lss, letter=letter, post=post, sample=sample

  on_error,2

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: vagc_combine_stripe_spec, stripe, outdir=outdir, /lrg, /southern, /lss, letter=letter, post=post, sample=sample'
      return
  ENDIF 

  ;; set up the system variables
  sdssidl_setup, /silent
  setup_mystuff

  time = systime(1)

  IF keyword_set(lss) AND keyword_set(lrg) THEN BEGIN 
      message,'The lss sample is MAIN only'
  ENDIF 
  IF keyword_set(lss) AND keyword_set(southern) THEN BEGIN 
      message,'The lss sample is MAIN only'
  ENDIF 


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; directories
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  byrun_dir = $
    sdssidl_config('spec_dir') + 'blanton/gal_collated/byrun_matched/'
  IF n_elements(outdir) EQ 0 THEN BEGIN 
      outdir = vagc_lensinput_dir(lss=lss)
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Output file names
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  catname = vagc_catname(lss=lss, letter=letter, post=post, sample=sample)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; look in the collate directory and get runs
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  files = findfile(byrun_dir, count=nf)
  w=where( strmatch(files,'*'+catname+'*'), nmatch)
  IF nmatch EQ 0 THEN BEGIN 
      message,'No '+catname+' files found'
  ENDIF 

  ;; Append for special samples

  IF keyword_set(lrg) THEN BEGIN 
      catname = 'lrg_'+catname
  ENDIF ELSE IF keyword_set(southern) THEN BEGIN 
      catname = 'southern_'+catname 
  ENDIF
  sname='stripe'+stripe2string(stripe)+'_'+catname+'.fit'
  sname = outdir+sname


  nf = nmatch
  files = files[w]
  runs = lonarr(nf)
 
  FOR i=0L, nf-1 DO BEGIN 
      t=( strsplit(files[i],'run',/extract) )[0]
      runstr = ( strsplit(t,'-',/extract) )[0]
      runs[i] = long(runstr)
  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Get stripes from run.par file
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  COMMON vagc_combine_block, runstruct

  IF n_elements(runstruct) EQ 0 THEN read_runlist, runstruct

  ;; a few won't match,but there are only 169 galaxies in those
  ;; as of 11-May-2004
  match, runs, runstruct.run, mruns, mstruct

  stripes = runstruct[mstruct].stripe
  strips = runstruct[mstruct].strip

  ;; Will not do this for lss?  Maybe they won't be in lss
  IF stripe EQ 10 THEN BEGIN 
      w=where(stripes EQ 10 OR stripes EQ 100, nrun)
  ENDIF ELSE BEGIN 
      w=where(stripes EQ stripe, nrun)
  ENDELSE 

  IF nrun EQ 0 THEN BEGIN 
      print
      print,'No runs for stripe '+ntostr(stripe)+'. Exiting'
      return
  ENDIF 

  runs = runs[mruns[w]]
  files = files[mruns[w]]
  stripes = stripes[w]
  strips = strips[w]
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Print some stuff
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  stripestrip = ntostr(stripes)+'-'+strips
  s=sort(stripestrip)
  
  print
  print,'---------------------------------------------------------'
  print,'Combining Stripe: '+ntostr(stripe)
  print
  print,'        Runs   Strips'
  forprint,runs[s],'    '+stripestrip[s]+'    '+files[s]
  print
  print,'---------------------------------------------------------'

  files = byrun_dir + files

  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set up header stuff
  ;;;;;;;;;;;;;;;;;;;;;;;;

  run_str = ntostr(long(runs))

  run_string = ''
  FOR irun=0L, nrun-1 DO BEGIN 
      run_string = run_string + run_str[irun]+' '
  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; read in all the runs
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  mrdfits_multi, files, struct
  print,'---------------------------------------------------------'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; choose galaxies
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ntot = n_elements(struct)
  noriginal = ntot

;  GALAXY_RED = 2L^5
;  GALAXY = 2L^6
;  MINZ_LRG = 0.15

;  IF keyword_set(lss) THEN BEGIN 

;      ngal = n_elements(struct)
;      wgal = lindgen(ngal)

;  ENDIF ELSE IF keyword_set(lrg) THEN BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; select LRG sample
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;      wgal = where(struct.z GT MINZ_LRG, ngal)

;      make_tsflag_struct, tsf
;      tsf.GALAXY_RED = 'Y'
;      tsf.SOUTHERN_SURVEY = 'N'
;      tsflag_select, struct, tsf, wgal, ngal, input_index=wgal

;      IF ngal EQ 0 THEN BEGIN
;          print
;          print,'No LRG galaxies in stripe: '+ntostr(stripe)
;          return
;      ENDIF ELSE BEGIN 
;          IF ngal LT ntot THEN struct = temporary(struct[wgal])
;          nremove = ntot-ngal
;          print
;          print,ntostr(nremove),$
;            ' Non-LRG or SOUTHERN galaxies removed from sample'
;          ntot = ngal
;      ENDELSE 

;  ENDIF ELSE IF keyword_set(southern) THEN BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; select SOUTHERN sample
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;      make_tsflag_struct, tsf
;      tsf.SOUTHERN_SURVEY = 'Y'
;      tsflag_select, struct, tsf, wgal, ngal

;      IF ngal EQ 0 THEN BEGIN
;          print
;          print,'No SOUTHERN galaxies in stripe: '+ntostr(stripe)
;          return
;      ENDIF ELSE BEGIN 
;          IF ngal LT ntot THEN struct = temporary(struct[wgal])
;          nremove = ntot-ngal
;          print
;          print,ntostr(nremove),$
;            ' Non-SOUTHERN galaxies removed from sample'
;          ntot = ngal
;      ENDELSE 


;  ENDIF ELSE BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; select MAIN sample
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;      make_tsflag_struct, tsf
;      tsf.GALAXY = 'Y'
;      tsf.SOUTHERN_SURVEY = 'N'
;      tsflag_select, struct, tsf, wgal, ngal

;      IF ngal EQ 0 THEN BEGIN
;          print
;          print,'No galaxies in stripe: '+ntostr(stripe)
;          return
;      ENDIF ELSE BEGIN 
;          IF ngal LT ntot THEN struct = temporary(struct[wgal])
;          nremove = ntot-ngal
;          print
;          print,ntostr(nremove),' Non-main galaxies removed from sample'
;          ntot = ngal
;      ENDELSE 

;  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; remove close pairs, preferentially keeping those that have matched
  ;; to the photo outputs
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  s = sort(struct.clambda)
;  struct = temporary(struct[s])

;  rmclose_lameta, struct.clambda, struct.ceta, keep, bad

;  nremove = ntot - n_elements(keep)
;  print,ntostr(nremove),' galaxies matched in overlap region'

;  struct = temporary(struct[keep])
;  ntot = n_elements(struct)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Output the file
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        
  print
  print,'Number in stripe: '+ntostr(ntot);+'/'+ntostr(noriginal)
  print
  print,'Writing combined file: ',sname

  ;; info in both headers
  outhdr = ['END']
  SXADDPAR, outhdr, 'RUN', run_string
  SXADDPAR, outhdr, 'STRIPE', ntostr(stripe)

  mwrfits2, struct, sname, outhdr, /create, /destroy


  print,'---------------------------------------------------------'
  print

  ptime,systime(1)-time

  return
END 

