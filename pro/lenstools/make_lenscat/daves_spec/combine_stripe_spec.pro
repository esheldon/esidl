;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    COMBINE_STRIPE
;       
; PURPOSE:
;    Combine shape files from different runs (but the same stripe) into
;    one big catalog.
;
; CALLING SEQUENCE:
;    combine_stripe, runs, reruns, clr, indir=indir
;
; INPUTS: 
;    runs: an integer vector containing the runs.
;    reruns: an integer vector with reruns.
;    clr: the bandpass to use in integer form (g=1, r=2, i=3)
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;    /lrg: create lrg sample
;       
; OUTPUTS: 
;    One big catalog for whole stripe.
;   
; CALLED ROUTINES:
;    SDSSIDL_SETUP
;    MRDFITS
;    SXADDPAR
;    RMCLOSE_LAMETA
;    MWRFITS
;
; PROCEDURE: 
;    
;	
;
; REVISION HISTORY:
;    1-June-2000  Erin Scott Sheldon UofMichigan
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


PRO combine_stripe_spec_addabs, struct

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Reddening and K-corrections
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Reddening correction'
  petro = struct.petrocounts - struct.reddening
  model = struct.counts_model - struct.reddening
  print,'K-correcting and calculating luminosity'
  kcorrect, petro, struct.petrocountserr, struct.z,$
            kmag, /sdssfix, kcorrectz=0.0
  kcorrect, model, struct.counts_modelerr, struct.z,$
            kmag_model, /sdssfix, kcorrectz=0.0
  

  kcorr = petro - kmag
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; absolute mags and luminosities
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  calc_sdss_absmag, kmag, struct.z, absmag, lum_solar
  calc_sdss_absmag, kmag_model, struct.z, absmag_model, lum_solar_model

  arrval = replicate(-9999., 5)

  ;; k-correction
  IF NOT tag_exist(struct[0], 'kcorr') THEN BEGIN 
      add_tag, temporary(struct), 'kcorr', arrval, struct
  ENDIF ELSE BEGIN 
      struct.kcorr = -9999.
  ENDELSE 

  ;; reddening and k-corrected counts
  IF NOT tag_exist(struct[0], 'abscounts') THEN BEGIN 
      add_tag, temporary(struct), 'abscounts', arrval, struct
  ENDIF ELSE BEGIN 
      struct.abscounts = -9999.
  ENDELSE 

  ;; absolute magnitude
  IF NOT tag_exist(struct[0], 'absmag') THEN BEGIN 
      add_tag, temporary(struct), 'absmag', arrval, struct
  ENDIF ELSE BEGIN 
      struct.absmag = -9999.
  ENDELSE 
  IF NOT tag_exist(struct[0], 'absmag_model') THEN BEGIN 
      add_tag, temporary(struct), 'absmag_model', arrval, struct
  ENDIF ELSE BEGIN 
      struct.absmag = -9999.
  ENDELSE 

  ;; luminosity in units of 10^10 L_{\odot}
  IF NOT tag_exist(struct[0], 'lum') THEN BEGIN 
      add_tag, temporary(struct), 'lum', arrval, struct
  ENDIF ELSE BEGIN 
      struct.lum = -9999.
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; copy in the good ones
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  FOR clr=0,4 DO BEGIN 
      w=where(struct.petrocounts[clr] NE 0.0 AND $
              kmag[clr,*] NE 1000., nw)
      IF nw NE 0 THEN BEGIN 
          struct[w].kcorr[clr]     = reform(kcorr[clr,w])
          struct[w].abscounts[clr] = reform(kmag[clr,w])
          struct[w].absmag[clr]    = reform(absmag[clr,w])
          struct[w].absmag_model[clr]    = reform(absmag_model[clr,w])
          struct[w].lum[clr]       = reform(lum_solar[clr,w])
      ENDIF 
      
  ENDFOR 
        
END 

PRO combine_stripe_spec, stripe, lrg=lrg


  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: combine_stripe_spec, stripe, /lrg'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Some parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  time = systime(1)
  colors=['u','g','r','i','z']

  ;; get the runs for this stripe
  read_stripe_index, sindex

  w=where(sindex.stripe EQ stripe, nst)
  IF nst EQ 0 THEN message,'No runs in stripe: '+ntostr(stripe)

  runs = sindex[w].run
  stripes = sindex[w].stripe
  strips = sindex[w].strip

  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set up runs, arrays
  ;;;;;;;;;;;;;;;;;;;;;;;;

  nrun = n_elements(runs)
  run_str = ntostr(long(runs))
  neach = lonarr(nrun)

  ;; set up the system variables
  sdssidl_setup
  setup_mystuff

  sdss_shapecorr_dir = sdssidl_config('shapecorr_dir')

  IF n_elements(indir) EQ 0 THEN $
    indir = SDSS_SHAPECORR_DIR+'spec_index/'
  IF n_elements(outdir) EQ 0 THEN $
    outdir = SDSS_SHAPECORR_DIR+'combined/'

  ;;;;;;;;;;;;;;;;;;;;;
  ;; Set up header
  ;;;;;;;;;;;;;;;;;;;;;

  run_string = ''
  FOR irun=0L, nrun-1 DO BEGIN 
      run_string = run_string + run_str[irun]+' '
  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; read in all the runs
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nnm = 'spec'
  nend = '.fit'

  IF keyword_set(lrg) THEN nnmout = 'lrg_spec' ELSE nnmout = 'spec'
  sname='stripe'+stripe2string(stripe)+'_'+nnmout+nend
  sname = outdir+sname

  print
  print,'---------------------------------------------------------'
  print,'Combining Stripe: '+ntostr(stripe)
  print
  print,'        Runs   Strips'
  colprint,runs,'     '+strips
  print
  print,'---------------------------------------------------------'

  FOR irun = 0L, nrun-1 DO BEGIN 
      
      name = indir + 'run'+ run_str[irun] +'_'+ nnm +nend
      print,'Reading ',name
          
      ;; read in one of the runs
      tstr = mrdfits(name, 1)
      sz_st = size(tstr, /structure)
      IF sz_st.type_name NE 'STRUCT' THEN BEGIN
          message,'Error reading run: '+run_str[irun]
      ENDIF 
      neach[irun] = sz_st.n_elements

      IF irun EQ 0 THEN BEGIN
          struct = tstr
      ENDIF ELSE BEGIN 
          concat_structs, temporary(struct), temporary(tstr), tmp_struct
          struct = temporary(tmp_struct)
      ENDELSE 

  ENDFOR 
  print,'---------------------------------------------------------'

  print
  colprint,'Stripe    Strip    Run       Number'
  colprint,ntostr(stripes),'         '+strips,'       '+run_str,neach

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; choose galaxies
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ntot = n_elements(struct)
  noriginal = ntot

  ;; note only 0.2% of SPEC_GALAXY objects
  ;; have no ECLASS measured

  SPEC_GALAXY = 2
  NOECLASS = 0.0
  GALAXY_RED = 2L^5
  GALAXY = 2L^6
  MINZ_LRG = 0.15

  IF keyword_set(lrg) THEN BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; select LRG sample
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      wgal = where(struct.spec_cln EQ SPEC_GALAXY AND $
                   struct.eclass NE NOECLASS AND $
                   ( (struct.primtarget AND GALAXY_RED) NE 0) AND $
                   struct.z GT MINZ_LRG, ngal)

      IF ngal EQ 0 THEN BEGIN
          print
          print,'No LRG galaxies in stripe: '+ntostr(stripe)
          return
      ENDIF ELSE BEGIN 
          IF ngal LT ntot THEN struct = temporary(struct[wgal])
          nremove = ntot-ngal
          print
          print,ntostr(nremove),' Non LRG/spec galaxies removed from sample'
          ntot = ngal
      ENDELSE 

  ENDIF ELSE BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; select MAIN sample
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      wgal = where(struct.spec_cln EQ SPEC_GALAXY AND $
                   struct.eclass NE NOECLASS AND $
                   ( (struct.primtarget AND GALAXY) NE 0), ngal)

      IF ngal EQ 0 THEN BEGIN
          print
          print,'No galaxies in stripe: '+ntostr(stripe)
          return
      ENDIF ELSE BEGIN 
          IF ngal LT ntot THEN struct = temporary(struct[wgal])
          nremove = ntot-ngal
          print
          print,ntostr(nremove),' Non spec galaxies removed from sample'
          ntot = ngal
      ENDELSE 

      ;; Don't use SOUTHERN_SURVEY stuff
      make_tsflag_struct,ts
      ts.SOUTHERN_SURVEY = 'N'
      tsflag_select, struct, ts, wgal, ngal
      
      IF ngal EQ 0 THEN BEGIN
          print
          print,'No Not-SOUTHERN_SURVEY galaxies in stripe: '+ntostr(stripe)
          return
      ENDIF ELSE BEGIN
          IF ngal LT ntot THEN struct = temporary(struct[wgal])
          nremove = ntot-ngal
          print
          print,ntostr(nremove),' SOUTHERN_SURVEY removed from sample'
          ntot = ngal
      ENDELSE 

  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; remove close pairs, preferentially keeping those that have matched
  ;; to the photo outputs
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  s = sort(struct.clambda)
  struct = temporary(struct[s])

  rmclose_lameta_matched, struct, keep, bad

  nremove = ntot - n_elements(keep)
  print,ntostr(nremove),' galaxies matched in overlap region'

  struct = temporary(struct[keep])
  ntot = n_elements(struct)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; add absolute mags if not there
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  combine_stripe_spec_addabs, struct
   
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Output the file
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        
  print
  print,'Number in stripe: '+ntostr(ntot)+'/'+ntostr(noriginal)
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
