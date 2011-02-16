; Order to run things
;    ::make_byrun - via scripts created by byrun_scripts
;    ::combine_stripe - via combine_allstripes
;    ::stuff
;    ::scat_write
;    ::make_htm_revind 
;
;

FUNCTION make_scat::init

  self->reset
  self.parstruct = ptr_new(/alloc)
  return,1

END 

PRO make_scat::reset

  self->free

;  self.parstruct = ptr_new(/alloc)

  self.cat = ptr_new(/alloc)
  self.keep = ptr_new(/alloc)

  self.band_keep[0] = ptr_new(/alloc)
  self.band_keep[1] = ptr_new(/alloc)
  self.band_keep[2] = ptr_new(/alloc)
  self.band_keep[3] = ptr_new(/alloc)
  self.band_keep[4] = ptr_new(/alloc)

  self.band_nkeep = 0

  self.final_keep = ptr_new(/alloc)
  self.final_nkeep = 0

  self.output_cat = ptr_new(/alloc)

END 

PRO make_scat::set_parstruct, parstruct
  ptr_free, self.parstruct
  self.parstruct = ptr_new(parstruct)
END 

FUNCTION make_scat::write_status
  return,self.write_status
END 
PRO make_scat::set_write_status, status
  self.write_status = status
END 

FUNCTION make_scat::status
  return,self.status
END 
PRO make_scat::set_status, status
  self.status = status
END 

PRO make_scat::help

  help,*self.cat
  help,*self.keep

END 


FUNCTION make_scat::parstruct_default, run=run, rerun=rerun, clrs=clrs

  IF n_elements(clrs) EQ 0 THEN clrs = [2,3]
  nclrs = n_elements(clrs)

  IF n_elements(run) EQ 0 THEN run = -1L
  IF n_elements(rerun) EQ 0 THEN rerun = -1

  parStruct = { $
                run: long(run[0]), $
                rerun: fix(rerun[0]), $
                camcol: 0, $
                clrs: clrs, $
                nclrs: nclrs, $
                hardrcut: 0.8, $
                hardprobcut: 0.8, $
                maxmag: 22.0, $
                minmag: 18.0, $
                edef: -9999.0, $
                e_errdef: 9999.0, $
                rotdef: 9999.0, $
                smdef: 9999.0, $
                seeingdef: 9999.0, $
                err_cut: 0.4*sqrt(nclrs) $
              }

  return,parstruct
END 

FUNCTION make_scat::getrunphotoz, run, rerun, status=status

  minid = ntostr( sdss_photoid(run, rerun, 1, 1, 0) )
  maxid = ntostr( sdss_photoid(run, rerun, 6, 5000, 5000) )
      
  query = $
    'SELECT photoid, photoz_z, photoz_zerr4 as photoz_zerr, photoz_zwarning '+$
    'FROM zphot '+$
    ' WHERE photoid BETWEEN '+minid+' AND '+maxid
      
  print,query
  pg=obj_new('postgres')
  zst = pg->query(query,status=status)
  
  return,zst

END 

PRO make_scat::main_cuts

  print
  print,'  Applying main cuts'
  ;; get parameter struct
  ps = *self.parstruct

  ;; temporarily set to standard variable
  cat = temporary( *self.cat )

  ncat = n_elements(cat)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Status cut
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; We will just cut stuff outside of its own stripe.  The resolution
  ;; within a stripe we will do ourselves.
  keep = $
    where( (cat.status AND sdss_flag('status', 'ok_stripe') ) NE 0, nkeep)

  pstr = ntostr( float(ncat-nkeep)/ncat*100., 5, /round)+'%'
  print,'    Threw out '+ntostr(ncat-nkeep)+' in status cut ('+pstr+')'


  IF nkeep NE 0 THEN BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Magnitude cut
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      w = where( cat[keep].cmodelmag_ext LE ps.maxmag AND $
                 cat[keep].cmodelmag_ext GE ps.minmag, nw)

      pstr = ntostr( float(nkeep-nw)/ncat*100., 5, /round)+'%'
      print,'    Threw out '+ntostr(nkeep-nw)+' in magnitude cut ('+pstr+')'
      nkeep = nw
      
      IF nkeep NE 0 THEN BEGIN 

          keep = keep[w]
          
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; flag cuts on each bandpass.  Each bandpass must pass.
          ;; This includes corrselect flags
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          nclr = n_elements(ps.clrs)
          FOR iclr=0L, nclr-1 DO BEGIN 
              IF keep[0] NE -1 THEN BEGIN 
                  w = where( (cat[keep].flags[iclr] AND sdss_flag('object1','deblended_as_psf') ) EQ 0 AND $
                             (cat[keep].corrselect_flags AND corrselect_flag(ps.clrs[iclr]) ) NE 0, nw)
                  IF nw EQ 0 THEN BEGIN 
                      keep = -1 
                  ENDIF ELSE BEGIN 
                      keep = keep[w]
                  ENDELSE 
              ENDIF 
          ENDFOR 
          pstr = ntostr( float(nkeep-nw)/ncat*100., 5, /round)+'%'
          print,'    Threw out '+ntostr(nkeep-nw)+$
            ' in flag/corrselect cuts ('+pstr+')'
          
          nkeep = nw
          
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Keep things with high galaxy probability
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          IF keep[0] NE -1 THEN BEGIN 
              
              w = where(cat[keep].objc_prob_gal GE ps.hardprobcut AND $
                        cat[keep].objc_prob_gal LE 1.0, nw)
              
              IF nw EQ 0 THEN BEGIN 
                  keep = -1
              ENDIF ELSE BEGIN 
                  keep = keep[w]
              ENDELSE 
              pstr = ntostr( float(nkeep-nw)/ncat*100., 5, /round)+'%'
              print,$
                '    Threw out '+ntostr(nkeep-nw)+' in probgal cuts ('+pstr+')'
              
              nkeep = nw
              
          ENDIF ;; passed flag cuts?
      ENDIF ;; passed magnitude cut?
  ENDIF ;; passed status cut?

  self.cat = ptr_new(cat, /no_copy)
  self.keep = ptr_new(keep, /no_copy)
  self.nkeep = nkeep

END 

PRO make_scat::plot

  cat = temporary( *self.cat )
  keep = temporary( *self.keep )

  plotrand, cat.clambda, cat.ceta, psym=3, frac=0.01
  plotrand, cat[keep].clambda, cat[keep].ceta, psym=3, color=!green, /over,$
    frac=0.01
  self.cat = ptr_new( cat, /no_copy )
  self.keep = ptr_new( keep, /no_copy )
END 

PRO make_scat::rcuts

  print
  print,'  Applying rcuts'
  ;; par struct
  ps = *self.parstruct

  ;; Make cuts in each bandpass
  cat = temporary( *self.cat )
  ncat = n_elements(cat)
  keep = temporary( *self.keep )
  nkeep = n_elements(keep)
  clrs = ps.clrs

  FOR iclr=0L, ps.nclrs-1 DO BEGIN 

      cstr = !colors[clrs[iclr]]

      bkeep = where(cat[keep].m_r_h[iclr] LT ps.hardrcut AND $
                    cat[keep].m_r_h[iclr] GT 0.0, nw)
      self.band_nkeep[iclr] = nw
      IF nw GT 0 THEN BEGIN 
          bkeep = keep[bkeep]
      ENDIF ELSE BEGIN 
          bkeep = -1L
      ENDELSE 

      pstr = ntostr( float(nkeep-nw)/ncat*100., 5, /round)+'%'
      print,'    Threw out '+ntostr(nkeep-nw)+' in '+cstr+' rcuts ('+pstr+')'


      self.band_keep[iclr] = ptr_new(bkeep, /no_copy)      
  ENDFOR 

  self.cat = ptr_new(cat, /no_copy)
  self.keep = ptr_new(keep, /no_copy)
  
END 

PRO make_scat::combine_shapes

  print
  print,'  Combining shapes'
  ;; parameter struct
  ps = *self.parstruct

  ;; Get these temporarily 
  cat = temporary( *self.cat )
  clrs = ps.clrs
  nclrs = ps.nclrs
  keep = temporary( *self.keep )

  ntotal = n_elements(cat)
  e1 = reform( replicate(ps.edef, nclrs, ntotal) )
  e2 = e1
  e1_recorr = e1
  e2_recorr = e2

  e1e1err = reform( replicate(ps.e_errdef, nclrs,  ntotal) )
  e1e2err = e1e1err
  e2e2err = e1e1err

  psfe1 = cat.m_e1_psf
  psfe2 = cat.m_e2_psf

  smear = reform( replicate(ps.smdef, nclrs, ntotal) )
  seeing = reform( replicate(ps.seeingdef, nclrs, ntotal) )
  corr = smear

  FOR iclr=0L, nclrs-1 DO BEGIN 

      cstr = !colors[clrs[iclr]]

      IF self.band_nkeep[iclr] GT 0 THEN BEGIN 

          ckeep = *self.band_keep[iclr]
          ;;help,ckeep
          e1[iclr, ckeep] = cat[ckeep].m_e1_corr_h[iclr]
          e2[iclr, ckeep] = cat[ckeep].m_e2_corr_h[iclr]

          e1e1err[iclr, ckeep] = cat[ckeep].m_e1e1err[iclr]
          e1e2err[iclr, ckeep] = cat[ckeep].m_e1e2err[iclr]
          e2e2err[iclr, ckeep] = cat[ckeep].m_e2e2err[iclr]
          smear[iclr, ckeep] = cat[ckeep].m_r_h[iclr]
          seeing[iclr, ckeep] = cat[ckeep].seeing[iclr]
          corr[iclr, ckeep] = 1./(1.-smear[iclr, ckeep])

          ;; correct the errors
          e1e1err[iclr, ckeep] = e1e1err[iclr, ckeep]*corr[iclr, ckeep]
          e1e2err[iclr, ckeep] = e1e2err[iclr, ckeep]*corr[iclr, ckeep]
          e2e2err[iclr, ckeep] = e2e2err[iclr, ckeep]*corr[iclr, ckeep]

          ;; Apply correction for linear egal vs. epsf
          print,'    corrected for egal vs. epsf: '+cstr
          fitstruct = egal_vs_epsf_read(ps.clrs[iclr])
          correct_eslope, fitStruct, $
            smear[iclr, ckeep], $
            e1[iclr, ckeep], e2[iclr, ckeep], $
            psfe1[iclr, ckeep], psfe2[iclr, ckeep], $
            e1c, e2c

          e1[iclr, ckeep] = e1c
          e2[iclr, ckeep] = e2c

      ENDIF 
  ENDFOR 

  ;; Apply rotations to the shapes
  print,'    Applying rotations'
  rotation = cat.rotation
  self->applyrot, rotation, nclrs, e1, e2, $
    e1e1err, e1e2err, e2e2err

  ;; Average the shapes
  print,'    Averaging shapes'
  combine_ellip_cove1e2, e1, e2, e1e1err, e1e2err, e2e2err, smear, seeing, $
    new_e1, new_e2, $
    new_e1e1err, new_e1e2err, new_e2e2err,$
    new_smear, new_seeing, $
    combine_flags, nave, gind, verbose=0, $
    err_cut=ps.err_cut, rcut=!hardrcut

  ;; gind now replaces keep
;  keep = gind

  IF gind[0] EQ -1 THEN BEGIN 
      print,'No objects bassed final cuts in combine_ellip'
      self.final_nkeep = 0
      self.cat = ptr_new( cat, /no_copy )
      self.keep = ptr_new( keep, /no_copy )
  ENDIF 

  nkeep = n_elements(gind)
  pstr = ntostr( float(nkeep)/ntotal*100., 5, /round)+'%'
  print,'    Keeping '+ntostr(nkeep)+' galaxies ('+pstr+')'

  ;; Remove unwanted measurements
  new_e1 = new_e1[gind]
  new_e2 = new_e2[gind]
  new_e1e1err = new_e1e1err[gind]
  new_e1e2err = new_e1e2err[gind]
  new_e2e2err = new_e2e2err[gind]
  new_smear = new_smear[gind]
  new_seeing = new_seeing[gind]
  combine_flags = combine_flags[gind]
  nave = nave[gind]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Make output structure
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  keep = gind
  nkeep = n_elements(keep)
  struct = self->byrun_structdef(nkeep)


  ;; Copy in
  struct.photoid = cat[keep].photoid
  
  struct.cmodelmag_r = cat[keep].cmodelmag_ext

  struct.e1 = temporary( new_e1 )
  struct.e2 = temporary( new_e2 )
  struct.e1e1err = temporary( new_e1e1err )
  struct.e1e2err = temporary( new_e1e2err )
  struct.e2e2err = temporary( new_e2e2err )
  
  struct.r = temporary( new_smear )
  struct.seeing = temporary( new_seeing )
  struct.combine_flags = temporary( combine_flags )
  struct.nave = temporary( nave )

  struct.clambda = cat[keep].clambda
  struct.ceta = cat[keep].ceta

  ;; put struct into self
  self.output_cat = ptr_new(struct, /no_copy)

  ;; Place these back into self
  self.cat = ptr_new( cat, /no_copy )
  self.final_keep = ptr_new( keep, /no_copy )
  self.final_nkeep = nkeep

END 

PRO make_scat::applyrot, rotation, $
        nclrs, e1, e2, e1e1err, e1e2err, e2e2err


  FOR clr=0,nclrs-1 DO BEGIN 

      w=where(e1e1err[clr,*] GT 0.0,nw)

      IF nw NE 0 THEN BEGIN 
          rotate_e1e2error, rotation[clr,w], e1[clr,w], e2[clr,w], $
                            e1e1err[clr,w],e1e2err[clr,w],e2e2err[clr,w],$
                            e1out,e2out,e1e1errout,e1e2errout,e2e2errout

          e1[clr,w] = e1out
          e2[clr,w] = e2out
          e1e1err[clr,w] = e1e1errout
          e1e2err[clr,w] = e1e2errout
          e2e2err[clr,w] = e2e2errout
          
      ENDIF 
  ENDFOR 

END 

FUNCTION make_scat::byrun_structdef, num

  ;; We are keeping clambda, ceta here for matching later
  struct = {$
             photoid: 0LL,$
             cmodelmag_r: 0.0, $
             e1: 0.0, $
             e2: 0.0, $
             e1e1err: 0.0, $
             e1e2err: 0.0, $
             e2e2err: 0.0, $
             r: 0.0, $
             seeing: 0.0, $
             combine_flags: 0L,$
             nave: 0b, $
             clambda: 0d, $
             ceta: 0d $
           }

  IF n_elements(num) NE 0 THEN BEGIN 
      struct = replicate(struct, num)
  ENDIF 

  return,struct
END 



FUNCTION make_scat::byrun_dir
  dir = self->source_dir()
  dir = concat_dir( dir, 'byrun' )
  return,dir
END 

FUNCTION make_scat::byrun_file, run, rerun

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: file = obj->byrun_file(run, rerun)'
      return,''
  ENDIF 

  sf = obj_new('sdss_files')
  file = 'scat-'+sf->run2string(run)+'-'+ntostr(rerun)+'.st'
  
  dir = self->byrun_dir()
  obj_destroy, sf

  return,concat_dir(dir,file)

END 


PRO make_scat::write_byrun_file

  ps = *self.parstruct
  camcol = ps.camcol
  struct = temporary( *self.output_cat )
  file = self->byrun_file(ps.run, ps.rerun)

  IF NOT fexist(file) THEN BEGIN 

      hdrStruct = {$
                    run: ps.run, $
                    rerun: ps.rerun, $
                    clrs: ps.clrs, $
                    hardrcut: ps.hardrcut, $
                    hardprobcut: ps.hardprobcut, $
                    maxmag: ps.maxmag, $
                    minmag: ps.minmag $
                  }

      print
      print,'  Writing to file: ',file
      write_idlstruct, struct, file, hdrStruct=hdrStruct, error=status
  ENDIF ELSE BEGIN 
      print
      print,'  Appending to file: ',file
      write_idlstruct, struct, file, error=status, /append
  ENDELSE 

  self->set_write_status, status
  self.output_cat = ptr_new( struct, /no_copy )

END 

FUNCTION make_scat::read_byrun_file, run, rerun

  file = self->byrun_file(run, rerun) 
  print
  print,'Reading file: ',file
  struct = read_idlstruct(file)
  return,struct

END 

PRO make_scat::make_byrun, run, rerun=rerun

  ;; Failure
  self->set_status, 1

  IF n_elements(rerun) EQ 0 THEN rerun = sdss_rerun(run)

  ;; read in what we need
  clrs = [2,3]

  ;; define parstruct
  parstruct = self->parstruct_default(run=run, rerun=rerun, clrs=clrs)
  ;; set global parstruct
  self->set_parstruct, parstruct

  rstr = ntostr(parstruct.run)
  rrstr = ntostr(parstruct.rerun)

  byrun_file = self->byrun_file(parStruct.run, parStruct.rerun)
  file_delete, byrun_file, /quiet

  ntotal = 0LL

  print,'-----------------------------------------------------------'
  print,'Processing Run = '+rstr+' Rerun = '+rrstr
  
  FOR camcol=1,6 DO BEGIN 

      colstr = ntostr(camcol)
      print
      print,'  ---------------------------------------------------------'
      print,'  Run = '+rstr+' Rerun = '+rrstr+' Camcol = '+colstr

      (*self.parstruct).camcol = camcol

      query = $
        'SELECT a.photoid, '+$
        'ae.status, '+$
        'ae.flags[2:3], '+$
        'a.corrselect_flags, '+$
        'a.objc_prob_gal, '+$
        'a.rotation[2:3], '+$
        'a.cmodelmag_ext[2], '+$
        'a.clambda, a.ceta, '+$
        'a.m_r_h[2:3], '+$
        'a.m_e1_corr_h[2:3], a.m_e2_corr_h[2:3], '+$
        'ae.m_e1e1err[2:3], ae.m_e1e2err[2:3], ae.m_e2e2err[2:3], '+$
        'ae.m_e1_psf[2:3], ae.m_e2_psf[2:3], '+$
        'a.seeing[2:3] '+$
        'FROM adatc AS a, adatc_extra AS ae '+$
        'WHERE a.run = '+rstr + $
        ' AND a.rerun = '+rrstr + $
        ' AND a.camcol = '+colstr + $
        ' AND a.photoid = ae.photoid'
      
;      print,query
      
      print
      print,'  Retrieving objects from database'
	  pg=obj_new('postgres')
      cat = pg->query(query)

      self.cat = ptr_new( cat, /no_copy )

      ;;  main cuts: status, magnitude, flags
      self->main_cuts

      IF self.nkeep GT 0 THEN BEGIN 
          ;; smear polarizability cuts
          self->rcuts

          IF total_int( self.band_nkeep ) GT 0 THEN BEGIN 
              ;; combine shapes, will apply additional cuts (err_cut)
              self->combine_shapes
          
              IF self.final_nkeep GT 0 THEN BEGIN 
                  ntotal = ntotal + self.final_nkeep
                  ;; Write to the file

                  self->write_byrun_file
                  IF self->write_status() NE 0 THEN BEGIN 
                      self->reset
                      print,'  * Write error'
                      return
                  ENDIF 
              ENDIF ELSE print,'  * No objects written for this camcol'
          ENDIF ELSE print,'  * No objects written for this camcol'
      ENDIF ELSE print,'  * No objects written for this camcol'
      
      self->reset

  ENDFOR 

  print
  print,'Wrote a total of '+ntostr(ntotal)+' objects to file: ',byrun_file
  print,'Done'
  self->set_status, 0

END 

PRO make_scat::runs, runs, reruns, nruns, $
             w_run_status=w_run_status, new=new

  ;; 5042 just has empty files, how stupid can you get
  w_run_status=where(!run_status.tsobj_photo_v GT 5.4 AND $
                     !run_status.run NE 94 AND $
                     !run_status.run NE 125 AND $
                     !run_status.run NE 1755 AND $
                     !run_status.run NE 3322 AND $
                     !run_status.run NE 5042, nruns)
  
  runs = !run_status[w_run_status].run
  reruns = !run_status[w_run_status].rerun

END 

FUNCTION make_scat::stripes, new=new

  self->runs, runs, reruns, w_run_status=w_run_status, new=new

  ;; get the unique stripes
  stripes = !run_status[w_run_status].stripe

  stripes = stripes[ rem_dup(stripes) ]
  
  return,stripes

END 



PRO make_scat::byrun_scripts, scriptnum

  IF 1 THEN BEGIN 

      self->runs, runs, reruns, nruns, /new
      eachdef = nruns/4
      
      w = lindgen(nruns)
      CASE scriptnum OF
          1: w = w[0:eachdef-1]
          2: w = w[eachdef:2*eachdef-1]
          3: w = w[2*eachdef:3*eachdef-1]
          4: w = w[3*eachdef:nruns-1]
          ELSE: message,'unknown script number: '+ntostr(scriptnum)
      ENDCASE 
      nruns = n_elements(w)
      runs = runs[w]
      reruns = reruns[w]

  ENDIF ELSE BEGIN 

      CASE scriptnum OF 
          1: runs = [1122, 1239, 1453, 1889, 1897, 2075, 2076, 2131, 2140, $
                     2207, 2305]
          2: runs = [2327, 2328, 2335, 2574, 2888, 2961, 3185]
          3: runs = [3478, 3513, 3605, 3911, 3927, 3958, 3965, 3967, 3972]
          4: runs = [3997, 4003, 4380, 4388, 4517, 4520, 4540, 4549, 4555, $
                     4571, 4600, 4633, 4671, 5045]
          ELSE: message,'unknown script number: '+ntostr(scriptnum)
      ENDCASE 
      reruns = sdss_rerun(runs)
      nruns = n_elements(runs)

  ENDELSE 

  dir = '/net/cheops2/home/esheldon/idlscripts/make_scat'

  file = concat_dir(dir,'make_byrun'+ntostr(scriptnum)+'.sh')

  openw, lun, file, /get_lun

  printf,lun,'#!/bin/bash'
  printf,lun
  FOR i=0L, nruns-1 DO BEGIN 

      run = runs[i]
      rerun = reruns[i]

      rstr = ntostr(run)
      rrstr = ntostr(rerun)

      printf,lun, 'idl<<EOF'
      printf,lun, '  run='+rstr
      printf,lun, '  rerun='+rrstr
      printf,lun, "  ss = obj_new('make_scat')"
      printf,lun, '  ss->make_byrun, run, rerun=rerun'
      printf,lun, '  if ss->status() ne 0 then exit, status=45'
      printf,lun, 'EOF'
      printf,lun, 'status=$?'
      printf,lun, 'if [ $status -eq 45 ]'
      printf,lun, 'then'
      printf,lun, '    echo Error in `basename $0` for '+rstr+' '+rrstr
      printf,lun, '    err="Write error for '+rstr+' '+rrstr+'"'
      printf,lun, '    echo "$err" | mail esheldon@cfcp.uchicago.edu -s "$err"'
      printf,lun, 'fi'
      printf,lun
  ENDFOR 
  printf,lun,'dt=`date`'
  printf,lun,'message="Finished making scat: '+ntostr(scriptnum)+'"'
  printf,lun,'echo "$message:  $dt" | mail esheldon@cfcp.uchicago.edu -s "$message"'

  free_lun, lun

  spawn,['chmod','755',file],/noshell

END 
















FUNCTION make_scat::bystripe_dir
  dir = self->source_dir()
  dir = concat_dir(dir, 'bystripe')
  return,dir
END 

FUNCTION make_scat::bystripe_file, stripe

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: file = obj->bystripe_file(stripe)'
      return,''
  ENDIF 

  sf = obj_new('sdss_files')
  file = 'scat-'+sf->stripe2string(stripe)+'.st'
  
  dir = self->bystripe_dir()
  obj_destroy, sf

  return,concat_dir(dir,file)

END 

FUNCTION make_scat::read_bystripe_file, stripe, hdrstruct=hdrstruct

  file = self->bystripe_file(stripe) 
  print
  print,'Reading file: ',file
  struct = read_idlstruct(file, hdrstruct=hdrstruct)
  return,struct

END 

FUNCTION make_scat::bystripe_structdef, num

  struct1 = self->byrun_structdef()

  struct2 = {$
              runcombine_flags: 0L, $
              nrunave: 0b, $
              photoz_z: -9999.0, $
              photoz_zerr: 9999.0, $
              photoz_zwarning: 9999, $
              maskflags: 0 $
            }

  struct = create_struct(struct1, struct2)

  IF n_elements(num) NE 0 THEN BEGIN 
      struct = replicate(struct, num)
  ENDIF 


  return,struct
END 


PRO make_scat::stripe_runs, stripe, runs, reruns, nruns

  ;; 5042 just has empty files, how stupid can you get
  w=where(!run_status.tsobj_photo_v GT 5.4 AND $
          !run_status.stripe EQ stripe AND $
          !run_status.run NE 94 AND $
          !run_status.run NE 125 AND $
          !run_status.run NE 1755 AND $
          !run_status.run NE 3322 AND $
          !run_status.run NE 5042, nruns)

  runs = !run_status[w].run
  reruns = !run_status[w].rerun

END 

PRO make_scat::apply_pixel_mask, struct, princeton=princeton

  nStruct = n_elements(struct)
  maskFlags = intarr(nStruct)
  
  IF keyword_set(princeton) THEN BEGIN 
      basic_mask = sdssidl_config('PIXEL_MASK_PRINCETON_BASIC')
  ENDIF ELSE BEGIN 
      basic_mask = sdssidl_config('PIXEL_MASK_BASIC')
      simple_mask = sdssidl_config('PIXEL_MASK_SIMPLE')
      bound_mask = sdssidl_config('PIXEL_MASK_BOUND')
  ENDELSE 

  IF n_elements(basic_mask) NE 0 THEN BEGIN 
      apply_pixel_mask, struct.clambda, struct.ceta, $
        basicMasked, basicUnmasked, maskfile=basic_mask
      IF basicMasked[0] NE -1 THEN BEGIN 
          maskFlags[basicMasked] = maskFlags[basicMasked] + esheldon_config('flags_masked_basic')
      ENDIF 
  ENDIF 

  IF n_elements(simple_mask) NE 0 THEN BEGIN 
      apply_pixel_mask, struct.clambda, struct.ceta, $
        simpleMasked, simpleUnmasked, maskfile=simple_mask
      IF simpleMasked[0] NE -1 THEN BEGIN 
          maskFlags[simpleMasked] = maskFlags[simpleMasked] + !FLAGS_MASKED_SIMPLE
      ENDIF 
  ENDIF 
 
  IF n_elements(bound_mask) NE 0 THEN BEGIN 
      apply_pixel_mask, struct.clambda, struct.ceta, $
        boundMasked, boundUnmasked, maskfile = bound_mask
      IF boundMasked[0] NE -1 THEN BEGIN 
          maskFlags[boundMasked] = maskFlags[boundMasked] + !FLAGS_MASKED_BOUND
      ENDIF 
  ENDIF 

  IF n_elements(combined_mask) NE 0 THEN BEGIN 
      apply_pixel_mask, struct.clambda, struct.ceta, $
        combinedMasked, combinedUnmasked, maskfile = combined_mask
      IF combinedMasked[0] NE -1 THEN BEGIN 
          maskFlags[combinedMasked] = $
            maskFlags[combinedMasked] + !FLAGS_MASKED_COMBINED
      ENDIF 
  ENDIF 


  struct.maskFlags = maskFlags

END 




PRO make_scat::combine_stripe, stripe

  stripe = stripe[0]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Read in all the runs
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,'-----------------------------------------------------------'
  print,'Combining stripe: '+ntostr(stripe)
  self->stripe_runs, stripe, runs, reruns, nruns
  files = self->byrun_file(runs, reruns)

  hdr = read_idlheader(files[0])

  struct = read_idlstruct_multi(files)
  ntot = n_elements(struct)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Create the new structure and copy in our the info we want
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  new_struct = self->bystripe_structdef(ntot)

  struct_assign, struct, new_struct, /verbose, /nozero

  ;; free old struct and rename this one to struct
  struct = 0
  struct = temporary(new_struct)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Get the photozs for each run and copy them in
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nmatch_tot = 0LL
  FOR i=0L, nruns-1 DO BEGIN 

      zst = self->getrunphotoz(runs[i], reruns[i], status=status)

      IF status EQ 0 THEN BEGIN 
          ;; match up to struct
          match, struct.photoid, zst.photoid, ms, mz, /sort
          
          IF mz[0] NE -1 THEN BEGIN 
              nmatch = n_elements(mz)
              struct[ms].photoz_z        = zst[mz].photoz_z
              struct[ms].photoz_zerr     = zst[mz].photoz_zerr
              struct[ms].photoz_zwarning = zst[mz].photoz_zwarning
          ENDIF ELSE nmatch = 0
          
          nmatch_tot = nmatch_tot + nmatch
          print,'Matched '+ntostr(nmatch)+' photozs for run '+ntostr(runs[i])
          
          zst = 0
      ENDIF 
  ENDFOR 
  print,'Matched a total of '+ntostr(nmatch_tot)+' photozs'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Stripe 82 is special.  We must deal with the many very low seeing
  ;; runs.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF stripe EQ 82 THEN BEGIN 

      minseeing = 0.9
      maxseeing = 1.5
      print
      print,'Stripe 82: '
      print,'  Removing objects with seeing < '+ntostr(minseeing)+$
        ' or seeing > '+ntostr(maxseeing)
      w = where(struct.seeing GE minseeing AND $
                struct.seeing LT maxseeing, nw, comp=comp, ncomp=ncomp)

      print,ntostr(ncomp)+' to be removed'

      tfile = '/net/cheops1/data6/db/tmp/tmp.bin'
      print
      print,'writing to temporary file: ',tfile
      openw, lun, tfile, /get_lun
      FOR i=0L, nw-1 DO writeu, lun, struct[w[i]]
      free_lun, lun


      tstruct = struct[0]
      zero_struct, tstruct
      struct = 0

      struct = replicate(tstruct, nw)
      
      openr, lun, tfile, /get_lun
      readu, lun, struct

      free_lun, lun

      file_delete, tfile, /quiet


  ENDIF 



  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Sort by lambda and combine.  We do it element-by-element to
  ;; conserve memory in stripe 82.  This only works for scalars
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,'Sorting by clambda for matching'
  s = sort(struct.clambda)
  tags = tag_names(struct)
  ntags = n_elements(tags)
  FOR i=0L, ntags-1 DO BEGIN 
      print,'Sorting tag: '+tags[i]
      struct.(i) = struct[s].(i)
  ENDFOR 
  

  rmclose_lameta_meane, struct, keep, out, $
    new_e1, new_e2, $
    new_e1e1err, new_e1e2err, new_e2e2err,$
    new_rsmear, new_seeing, $
    new_photoz_z, new_photoz_zerr, new_photoz_zwarning, $
    new_cflags, new_nave, $
    nkeep=nkeep

  struct = struct[keep]

  struct.e1 = new_e1
  struct.e2 = new_e2
  struct.e1e1err = new_e1e1err
  struct.e1e2err = new_e1e2err
  struct.e2e2err = new_e2e2err

  struct.r = new_rsmear
  struct.seeing = new_seeing

  struct.photoz_z = new_photoz_z
  struct.photoz_zerr = new_photoz_zerr
  struct.photoz_zwarning = new_photoz_zwarning

  struct.runcombine_flags = new_cflags
  struct.nrunave = new_nave

  nremove = ntot - nkeep
  print,ntostr(nremove),' galaxies matched in overlap region'
  ntot = n_elements(keep)
  
      
  print,'Final number of objects: '+ntostr(ntot)
  ;; Make the new header
  new_hdr = {$
            stripe: stripe, $
            runs: runs, $
            reruns: reruns, $
            clrs: fix(hdr.clrs), $
            hardrcut: hdr.hardrcut, $
            hardprobcut: hdr.hardprobcut, $
            minmag: hdr.minmag, $
            maxmag: hdr.maxmag  $
            }
      
  file = self->bystripe_file(stripe)
  print
  print,'Writing to file: ',file
  write_idlstruct, struct, file, hdrStruct=new_hdr
  
  struct = 0
END 

PRO make_scat::combine_allstripes, stripes=stripes
  
  IF n_elements(stripes) EQ 0 THEN stripes = self->stripes()
;  stripes = [17,18,19,20,21,22,23]

;  w = where(stripes GT 23)
;  stripes = stripes[w]

  nstripes = n_elements(stripes)
  FOR i=0L, nstripes-1 DO self->combine_stripe, stripes[i]

END 




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Final Catalogs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION make_scat::source_dir
  dir = concat_dir( esheldon_config('lensinput_dir'),'srcgal/chicago')
  dir = expand_tilde(dir)
  return,dir
END 

FUNCTION make_scat::source_file, sample
  IF n_elements(sample) EQ 0 THEN BEGIN 
      message,'You must enter sample number'
  ENDIF 
  CASE sample OF
      1: BEGIN 
          file = 'scat.st'
          dir = self->source_dir()
          return, concat_dir(dir,file)
      END 
      ELSE: message,'Unknown source sample: '+strn(sample)
  ENDCASE 
END 

PRO make_scat::source_write, struct, sample, append=append, hdr=hdr
  file = self->source_file(sample)
  print
  print,'Writing source file: ',file
  write_idlstruct, struct, file, append=append, hdr=hdr
END 
FUNCTION make_scat::source_get, sample
  file = self->source_file(sample)
  print
  print,'Reading source file: ',file
  struct = read_idlstruct(file)
  return,struct
END 


FUNCTION make_scat::stripes_use, nstripe

  stripes = self->stripes()

  w = where(stripes LE 39 OR $
            stripes EQ 76 OR $
            stripes EQ 82 OR $
            stripes EQ 86, nstripe)
  stripes = stripes[w]

  return,stripes

END 


FUNCTION make_scat::source_structdef, num

  struct = {photoid: 0LL, $
            run: 0, $
            rerun: 0, $
            camcol: 0, $
            field: 0, $
            id:0, $
            stripe: 0, $
            clambda:0d, $
            ceta:0d, $
            e1:0.0, $
            e2:0.0, $
            e1e1err:0.0, $
            e1e2err:0.0, $
            e2e2err:0.0, $
            seeing: 0.0, $
            r:0.0, $
            combine_flags: 0L, $
            nave: 0b, $
            nrunave: 0b, $
            cmodelmag_r:0.0, $
            photoz_z:0.0, $
            photoz_zerr:0.0, $
            maskflags: 0, $
            htm_index: 0L $
           }

  IF n_elements(num) NE 0 THEN BEGIN 
      retstruct = replicate(struct, num)
      return,retstruct
  ENDIF ELSE BEGIN 
      return,struct
  ENDELSE 

END 



;; Combine all stripes into one big cat
PRO make_scat::combine

  htm_depth = 10

  ;; sample 1
  sample = 1
  outfile = self->source_file(sample)
  print
  print,'Will write to file: ',outfile

  outst = self->source_structdef()

  stripes = self->stripes_use(nstripe)
  sstr = ntostr(stripes)

  FOR i=0L, nstripe-1 DO BEGIN 

      print,'--------------------------------------------------------'
      st = self->read_bystripe_file(stripes[i])
      ntot = n_elements(st)

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Only use objects with good photozs
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      print
      print,'Removing bad photoz'
      keep = where(st.photoz_zwarning EQ 0, nkeep)
      IF nkeep NE 0 THEN BEGIN 

          nremove = ntot - nkeep
          st = st[keep]
          print,'Removed '+ntostr(nremove)+' in photoz cuts'



          ;; Only keep those which pass the basic maskflags
          ;; Apply the mask and save the maskflags
          print
          print,'Getting maskflags'
          self->apply_pixel_mask, st

          keep2 = where((st.maskflags AND !FLAGS_MASKED_BOUND) EQ 0, nkeep2)

          IF nkeep2 NE 0 THEN BEGIN 
              
              nremove = nkeep - nkeep2
              st = st[keep2]
              print,'Removed '+ntostr(nremove)+' in mask'

              print,'Keeping '+ntostr(nkeep2)+'/'+ntostr(ntot)

              out_struct = replicate(outst, nkeep2)
          
              copy_struct, st, out_struct
              st = 0
              out_struct.stripe = stripes[i]
          
              ;; htm
              csurvey2eq, out_struct.clambda, out_struct.ceta, ra, dec
          
              out_struct.htm_index = htm_index(ra, dec, htm_depth)
          
              IF i EQ 0 THEN BEGIN 
                  print,'Writing to file: ',outfile
                  write_idlstruct, out_struct, outfile
              ENDIF ELSE BEGIN 
                  print,'Appending to file: ',outfile
                  write_idlstruct, out_struct, outfile, /append
              ENDELSE 

              out_struct = 0
          ENDIF ELSE BEGIN 
              print,'No objects in stripe '+sstr[i]+' passed mask cuts!'
          ENDELSE 
      ENDIF ELSE BEGIN 
          print,'No objects in stripe '+sstr[i]+' passed photoz cuts!'
      ENDELSE 

  ENDFOR 
  


END 

;; Combine all stripes into one big postgres table
PRO make_scat::stuff


  htm_depth = 10

  table = 'scat'
  print,'Will stuff to table: ',table

  dir = self->source_dir()

  outst = self->source_structdef()

  stripes = self->stripes_use(nstripe)
  sstr = ntostr(stripes)

  FOR i=0L, nstripe-1 DO BEGIN 

      print,'--------------------------------------------------------'
      st = self->read_bystripe_file(stripes[i])
      ntot = n_elements(st)

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Only use objects with good photozs
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      print
      print,'Removing bad photoz'
      keep = where(st.photoz_zwarning EQ 0, nkeep)
      IF nkeep NE 0 THEN BEGIN 

          nremove = ntot - nkeep
          st = st[keep]
          print,'Removed '+ntostr(nremove)+' in photoz cuts'



          ;; Only keep those which pass the basic maskflags
          ;; Apply the mask and save the maskflags
          print
          print,'Getting maskflags'
          self->apply_pixel_mask, st

          keep2 = where((st.maskflags AND !FLAGS_MASKED_BOUND) EQ 0, nkeep2)

          IF nkeep2 NE 0 THEN BEGIN 
              
              nremove = nkeep - nkeep2
              st = st[keep2]
              print,'Removed '+ntostr(nremove)+' in mask'

              print,'Keeping '+ntostr(nkeep2)+'/'+ntostr(ntot)

              out_struct = replicate(outst, nkeep2)
          
              copy_struct, st, out_struct
              st = 0
              out_struct.stripe = stripes[i]
          
              ;; htm
              csurvey2eq, out_struct.clambda, out_struct.ceta, ra, dec
          
              out_struct.htm_index = htm_index(ra, dec, htm_depth)

              self->struct2table, out_struct, table, status=stuff_status,$
                tmpdir=dir, primary_key='photoid'
          
              out_struct = 0
          ENDIF ELSE BEGIN 
              print,'No objects in stripe '+sstr[i]+' passed mask cuts!'
          ENDELSE 
      ENDIF ELSE BEGIN 
          print,'No objects in stripe '+sstr[i]+' passed photoz cuts!'
      ENDELSE 

  ENDFOR 
  

END 


PRO make_scat::create_index
  queries = $
    [ $
      'create index scat_stripe_index on scat (stripe)',$
      $
      'create index scat_cmodelmag_r_index on scat (cmodelmag_r)', $
      $
      'create index scat_photoz_zerr_index on scat (photoz_zerr)', $
      'create index scat_photoz_z_index on scat (photoz_z)', $
      $
      'create index scat_r_index on scat (r)', $
      'create index scat_seeing_index on scat (seeing)', $
      $
      'create index scat_maskflags_index on scat (maskflags)', $
      $
      'analyze scat' $
    ]

  nq = n_elements(queries)
  FOR i=0L, nq-1 DO BEGIN 

      print,queries[i]
      self->postgres::query, queries[i], status=status

      IF status NE self->postgres::status_val('no_result') THEN BEGIN 
          message,'Failed'
      ENDIF 

  ENDFOR 

END 

FUNCTION make_scat::where_clause, sample

  IF n_elements(sample) EQ 0 THEN $
    message,'You must send a sample number'

  CASE sample OF
      1: clause=''
      2: clause='( r >  0.585 )'  ; well resolved
      3: clause='( r <= 0.585 )'  ; less well resolved
      4: clause='( seeing < 1.5 )'
      6: clause='( r > 0.33333 )'
      7: clause='( r > 0.66666 )' ; this is actuallly the 2/3 point
      8: clause='( r > 0.4 )'
      ELSE: message,'Unsupported sample '+ntostr(sample)
  ENDCASE 
  return,clause

end

PRO make_scat::scat_write

  sample = 1
  outfile = self->source_file(sample)
  print
  print,'Will write to file: ',outfile

  query_front = $
    'SELECT '+$
    'stripe, '+$
    'e1, e2, '+$
    'e1e1err, e1e2err, e2e2err, '+$
    'clambda, ceta, '+$
    'photoz_z, photoz_zerr, '+$
    'htm_index FROM scat WHERE stripe = '

  ;; Do this in two chunks because of memory limitations

  stripes = self->stripes_use(nstripe)

  pg=obj_new('postgres')
  FOR i=0L, nstripe-1 DO BEGIN 
      stripe = stripes[i]

      query = query_front + ntostr(stripe)

      print
      print,query

      scat = pg->query(query, status=status)
      IF status NE self->postgres::status_val('success') THEN BEGIN 
          message,'Could not retrieve source info'
      ENDIF 
      
      IF i GT 0 THEN append = 1

      print,'Writing to file: ',outfile
      write_idlstruct, scat, outfile, append=append
      scat = 0

  ENDFOR 


END 




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Make htm reverse indices file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION make_scat::htmrev_file, sample
  file = self->source_file(sample)
  file = repstr(file, '.st', '-htmrev.bin')
  return,file
END 

PRO make_scat::make_htm_revind

  query = 'SELECT htm_index from scat'
  print,query
  pg=obj_new('postgres')
  st  = pg->query(query)
  htm_index = st.htm_index
  st = 0

  print
  print,'Creating reverse indices'
  
  minid = min(htm_index, max=maxid)
  h = histogram( htm_index-minid, min=0, rev=rev )
  nrev = n_elements(rev)
  
  ;; sample 1
  rev_file = self->htmrev_file(1)
  
  print
  print,'Writing to rev file: ',rev_file
  openw, lun, rev_file, /get_lun
  writeu, lun, nrev
  
  writeu, lun, minid
  writeu, lun, maxid
  
  writeu, lun, rev
  
  free_lun, lun
  
  rev = 0


END 





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Some plots
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO make_scat::plot_lameta, struct, w=w, stripe82=stripe82, _extra=_extra

  xtitle = !csym.lambda+'!Dc!N'
  ytitle = !csym.eta+'!Dc!N'

  IF keyword_set(stripe82) THEN BEGIN 
      w=where(struct.clambda GE -60.0 AND struct.clambda LE 60.0)

      xrange = [-60.0, 80.0]
      yrange = [144, 151]

      plotrand, struct[w].clambda, struct[w].ceta, $
        fracuse=0.01, psym=3, $
        /ynozero, xrange=xrange, yrange=yrange, xstyle=1, ystyle=1, $
        xtitle=xtitle, ytitle=ytitle, $
        _extra=_extra


  ENDIF ELSE BEGIN 

      plotrand, struct.clambda, struct.ceta, $
        xtitle=xtitle, ytitle=ytitle, $
        fracuse=0.01, psym=3, /ynozero, xstyle=1, $
        _extra=_extra
      
  ENDELSE 

END 

PRO make_scat::plot_masked, region, struct, _extra=_extra

  ;; This "region" is just northern and southern

  ;; stripes 61,62 are not in Ryan's window

  plot_dir = '~/plots/make_scat'
  pngfile = concat_dir(plot_dir, 'masks_region'+ntostr(region)+'.png')
  title = 'region '+ntostr(region)
  stripes = self->stripes()
  CASE region OF
      1: BEGIN 
          res = [750, 800]
          w=where(stripes LT 50)
      END 
      2: BEGIN 
          res = [860, 570]
          w=where(stripes GT 70)
      END 
      ELSE: message,'No such region: '+ntostr(region)
  ENDCASE 
  
  stripes = stripes[w]
  IF n_elements(struct) EQ 0 THEN BEGIN 

      files = self->bystripe_file(stripes)
      struct = $
        read_idlstruct_multi(files,columns=['clambda','ceta','maskflags'])

  ENDIF 

  setupplot, 'Z'
  simpctable, rmap, gmap, bmap
  device, set_resolution = res

  self->plot_lameta, struct, _extra=_extra, title=title

  colors = [!p.color, !yellow, !magenta, !green, !red]

  w = where( (struct.maskflags AND !FLAGS_MASKED_COMBINED) NE 0, nw)
  IF nw GT 0 THEN BEGIN 
      oplot, struct[w].clambda, struct[w].ceta, psym=3, color=colors[4]
  ENDIF 

  w = where( (struct.maskflags AND !FLAGS_MASKED_BOUND) NE 0, nw)
  IF nw GT 0 THEN BEGIN 
      oplot, struct[w].clambda, struct[w].ceta, psym=3, color=colors[3]
  ENDIF 

  w = where( (struct.maskflags AND !FLAGS_MASKED_SIMPLE) NE 0, nw)
  IF nw GT 0 THEN BEGIN 
      oplot, struct[w].clambda, struct[w].ceta, psym=3, color=colors[2]
  ENDIF 


  w = where( (struct.maskflags AND !FLAGS_MASKED_BASIC) NE 0, nw)
  IF nw GT 0 THEN BEGIN 
      oplot, struct[w].clambda, struct[w].ceta, psym=3, color=colors[1]
  ENDIF 




  message = ['all','basic','simple','bound','combined']
  legend, message, colors=colors, psym=8, box=0, /right

  pngfile = '~/plots/make_scat/masks_region'+ntostr(region)+'.png'
  print
  print,'Writing pngfile: ',pngfile
  write_png, pngfile, tvrd(), rmap, gmap, bmap
    

  setupplot,'X'


END 

PRO make_scat::plot_nrunave, struct, stripe82=stripe82, _extra=_extra

  simpctable, colorlist=colors

  self->plot_lameta, struct, stripe82=stripe82, _extra=_extra

  h = histogram( struct.nrunave, $
                 min=2, max=max(struct.nrunave), rev=rev )

  nh = n_elements(h)

  FOR i=0L, nh-1 DO BEGIN 

      IF rev[i] NE rev[i+1] THEN BEGIN 

          w = rev[ rev[i]:rev[i+1] -1 ]          
          nrunave = struct[w[0]].nrunave
          nw = n_elements(w)
          
          print,'nrunave = ',nrunave,nw


          IF nw EQ 1 THEN BEGIN 
              clambda = [struct[w].clambda]
              ceta = [struct[w].ceta]
          ENDIF ELSE BEGIN 
              clambda = struct[w].clambda
              ceta = struct[w].ceta
          ENDELSE 

          IF nw LT 200 THEN psym=8 ELSE psym=3

          oplot, clambda, ceta, $
            psym=psym, color=colors[nrunave-1]


      ENDIF 

  ENDFOR 

  maxnrunave = nh-1
  mess = 'nave = '+ntostr(lindgen(nh+1)+1)
  lcolors = colors[0:nh]

  legend, mess, psym=8, colors=lcolors, /right, box=0

END 





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Princeton stuff
; 
; Steps to create a catalog:
;   ->princeton_stuff
;   ->princeton_photoz_match
;   ->princeton_sub_sample
;   ->princeton_scat_write
;   ->princeton_make_htm_revind
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


FUNCTION make_scat::princeton_source_dir
  dir = concat_dir(esheldon_config("lensinput_dir"),'srcgal/princeton')
  dir = expand_tilde(dir)
  return,dir
END 
FUNCTION make_scat::princeton_bystripe_dir
  dir = self->princeton_source_dir()
  return,concat_dir(dir,'bystripe')
END 
FUNCTION make_scat::princeton_bystripe_infile, stripe, all=all

  dir = self->princeton_bystripe_dir()
  IF keyword_set(all) THEN BEGIN 
      fileregex = concat_dir(dir,'stripe-[0-9][0-9].extra.fits.gz')
      files = file_search(fileregex)
      return,files
  ENDIF 
  IF n_elements(stripe) EQ 0 THEN BEGIN 
      message,'-Syntax: file = ms->princeton_srcfile(stripe)'
  ENDIF 

  return,concat_dir(dir,'stripe-'+stripe2string(stripe)+'.extra.fits.gz')
END 
FUNCTION make_scat::princeton_bystripe_infile_get, stripe
  file = self->princeton_bystripe_infile(stripe)
  print,'Reading file: ',file
  return,mrdfits(file,1)
END 

FUNCTION make_scat::princeton_bystripe_file, stripe, all=all

  dir = self->princeton_bystripe_dir()
  IF keyword_set(all) THEN BEGIN 
      fileregex = concat_dir(dir,'stripe-[0-9][0-9].st')
      files = file_search(fileregex)
      return,files
  ENDIF 
  IF n_elements(stripe) EQ 0 THEN BEGIN 
      message,'-Syntax: file = ms->princeton_srcfile(stripe)'
  ENDIF 

  return,concat_dir(dir,'stripe'+stripe2string(stripe)+'.st')
END 
FUNCTION make_scat::princeton_bystripe_get, stripe,status=status
  file = self->princeton_bystripe_file(stripe)
  print,'Reading file: ',file
  return,read_idlstruct(file, status=status)
END 

FUNCTION make_scat::princeton_source_file, sample, photoz_match=photoz_match, equatorial=equatorial
    dir = self->princeton_source_dir()
    sample_string = self->source_sample_string(sample)
    if keyword_set(equatorial) then begin
        front = 'princeton-scat-eq'
        ext='.fit'
    endif else begin
        front = 'princeton-scat'
        ext='.st'
    endelse
    if keyword_set(photoz_match) then sample_string = sample_string+'-pzmatch'
    file = concat_dir(dir,front+'-'+sample_string+ext)
    return,file
END 
FUNCTION make_scat::princeton_source_get, sample, equatorial=equatorial
  file = self->princeton_source_file(sample, equatorial=equatorial)
  print
  print,'Reading source file: ',file
  struct = read_idlstruct(file)
  return,struct
END 


;; read in rachel's bystripe catalogs, make cuts, convert to mags, etc.
;; and re-output by stripe

pro make_scat::princeton_rotate2csurvey, $
  struct, e1csurvey, e2csurvey, $
  wgood=wgood, ngood=ngood, $
  e1pix=e1pix, e2pix=e2pix

  ;rotdir = '/net/cheops1/data7/imaging.local/princeton_surveyrot/'
  rotdir = '/global/early2/esheldon/princeton_surveyrot'
  nobj = n_elements(struct)
  e1csurvey = replicate(-9999.0, nobj)
  e2csurvey = e1csurvey

  ;ss = obj_new('sdss_survey_rot')

  ;; rotate the shapes to pixel coordinates
  rotate_e1e2, -struct.phi_offset, struct.e1, struct.e2, e1pix, e2pix
  
  ;; Make tree structure of run,rerun,camcol,field
  idstruct = sdss_histid(struct.run, struct.rerun, struct.camcol, struct.field)
  n_unique = idStruct.Nleaves
  ptrlist = ptrarr(n_unique)
  ptrIndex = 0L
      
  srot = obj_new('sdss_survey_rot')

  pruns = idStruct.runs
  for ri=0l, idstruct.nruns-1 do begin 
      run = (*pruns)[ri].run 

      preruns = (*pruns)[ri].reruns
      for rri=0l, (*pruns)[ri].nreruns-1 do begin 
          rerun = (*preruns)[rri].rerun

          ;query='select * from field_rotation where run='+ntostr(run)+' and rerun='+ntostr(rerun)
          ;print,query
          ;unrot = self->query(query, status=status)
          runrot = srot->read('survey', run, rerun=301, status=status)
          if status eq 0 then begin
          
              pcamcols = (*preruns)[rri].camcols
              for ci=0l, (*preruns)[rri].ncamcols-1 do begin 
                  camcol = (*pcamcols)[ci].camcol
                  fieldStructs = *(*pcamcols)[ci].fields
                  
                  nf = n_elements(fieldStructs)
                  for fi=0l, nf-1 do begin 
                      field = fieldStructs[fi].field
                      wfield = *fieldStructs[fi].indices
                      
                      w = where(runrot.camcol EQ camcol AND $
                                runrot.field EQ field, nw)
                      if nw ne 0 then begin 

                          ;; WARNING: just using r-band angle
                          rotate_e1e2, $
                            runrot[w].angle[2], $
                            e1pix[wfield], e2pix[wfield], $
                            e1surv, e2surv
                          
                          e1csurvey[wfield] = e1surv
                          e2csurvey[wfield] = e2surv
                          
                      endif else splog,form='("camcol/field not found: ",i0," ",i0)',camcol,field
                  endfor     ;; fields
              endfor         ;; camcols
          endif else splog,form='("rotation info not found for run: ",i0," rerun: ",i0)',run,rerun
      endfor                 ;; reruns
  endfor                     ;; runs

  wgood = where(e1csurvey NE -9999.0, ngood)
  ;obj_destroy,ss

end 

FUNCTION make_scat::princeton_bystripe_structdef
  st = self->bystripe_structdef()
  st = create_struct(st, 'OBJC_FLAGS', 0L, 'OBJC_FLAGS2', 0L)
  return,st
END 
FUNCTION make_scat::princeton_stuff_structdef, num

  st = {$
         photoid: 0LL, $
         run: 0, $
         rerun: 0, $
         camcol: 0, $
         field: 0, $
         id: 0, $
         stripe: 0, $
         objc_flags: 0L,$
         objc_flags2: 0L, $
         e1pix: 0.0, $
         e2pix: 0.0, $
         e1eq: 0.0, $
         e2eq: 0.0, $
         e1survey:0.0, $
         e2survey:0.0, $
         e1e1err:0.0, $
         e1e2err:0.0, $
         e2e2err:0.0, $
         seeing: 0.0, $
         r:0.0, $
         r_r: 0.0, $
         r_i: 0.0, $
         modelmag: fltarr(5), $  ; log mags
         cmodelmag: fltarr(5),$  ; log mags
         ra:0d, $
         dec: 0d, $
         clambda:0d, $
         ceta:0d, $
         maskflags: 0, $
         htm_index: 0L $
       }

  IF n_elements(num) NE 0 THEN BEGIN 
      st = replicate(st, num)
  ENDIF 

  return,st
END 

;; no galaxy probability cuts yet
pro make_scat::princeton_stuff


    dophotoz=0

    table = 'scat_princeton'

    nclrs = 2
    err_cut = 0.4*sqrt(nclrs)
    etot_cut = 4.0

    htm_depth = 10

    ;; define parstruct
    ps = self->parstruct_default(run=run, rerun=rerun, clrs=clrs)
    ;; set global parstruct
    self->set_parstruct, ps

    ps = *self.parstruct

    files = self->princeton_bystripe_infile(/all)
    tmpdir = self->princeton_source_dir()

    nstripe = n_elements(files)
    FOR i=0L, nstripe-1 DO BEGIN 

        ;; get stripe
        dirsep,files[i],tdir,tf
        stripe = fix( (strsplit((strsplit(tf,'-',/extract))[1],'.',/extract))[0])


        print,'-------------------------------------------------------------'
        print,'Reading file: ',files[i]
        t = mrdfits(files[i], 1)
        ncat = n_elements(t)
        nkeepold = ncat

        ;; Output structure
        outstruct = self->princeton_stuff_structdef(ncat)
        copy_struct, t, outstruct

        ;; unique identifier
        outstruct.photoid = sdss_photoid(t) 

        ;; copy in stripe
        outstruct.stripe = stripe

        ; errors: Rachel says there was an error processing, so we need
        ; to multiply by 2
        outstruct.e1e1err = t.sigmae*2
        outstruct.e1e2err = 0.0
        outstruct.e2e2err = t.sigmae*2

        ;; seeing in r
        outstruct.seeing = sqrt(t.m_rr_cc_psf[0]/2.0)*2.35*0.4

        ;; htm index
        outstruct.htm_index = htm_index(t.ra, t.dec, htm_depth)

        ;; corrected survey coords
        eq2csurvey, t.ra, t.dec, clambda, ceta
        outstruct.clambda = clambda 
        outstruct.ceta = ceta

        ;; correct e1/e2 for this coordinate system
        print,'Rotating with r-band angle'
        self->princeton_rotate2csurvey, t, e1csurvey, e2csurvey, $
            e1pix=e1pix, e2pix=e2pix, wgood=wgood


        outstruct.e1pix = e1pix
        outstruct.e2pix = e2pix
        outstruct.e1eq = t.e1
        outstruct.e2eq = t.e2
        outstruct.e1survey = e1csurvey
        outstruct.e2survey = e2csurvey

        outstruct.modelmag  = 22.5-2.5*alog10(t.modelflux > 0.001)
        outstruct.cmodelmag = 22.5-2.5*alog10(t.cmodelflux > 0.001)

        ;; average r
        outstruct.r = total(t.r2, 1)/2.0
        outstruct.r_r = t.r2[0]
        outstruct.r_i = t.r2[1]

        ;; pixel mask
        print,'Applying pixel mask'
        self->apply_pixel_mask, outstruct, /princeton

        keep = where((outstruct.maskflags and esheldon_config('flags_masked_basic')) eq 0, nkeep)

        pstr = ntostr( float(nkeepold-nkeep)/ncat*100., 5, /round)+'%'
        print,'    Threw out '+ntostr(nkeepold-nkeep)+' in pixel masks ('+pstr+')'

        IF nkeep NE 0 THEN BEGIN 

            nkeepold = nkeep
            t=t[keep]
            outstruct=outstruct[keep]

            ; Rachel already extinction corrected.  Note, there is actually already
            ; a modelmag[2] < 22 cut on the catalog, so this is redundant

            keep = where(outstruct.cmodelmag[2] LE ps.maxmag AND $
                          outstruct.cmodelmag[2] GE ps.minmag and $
                          outstruct.modelmag[2] le 22.5 and $
                          outstruct.modelmag[2] ge ps.minmag, nkeep)

            pstr = ntostr( float(nkeepold-nkeep)/ncat*100., 5, /round)+'%'
            print,'    Threw out '+ntostr(nkeepold-nkeep)+' in magnitude cut ('+pstr+')'

            IF nkeep GT 0 THEN BEGIN 

                nkeepold = nkeep
                t=t[keep]
                outstruct=outstruct[keep]

                ;; err_cut, etot_cut
                ;; r2 cut: they have already cut > 1/3 in at least one band,
                ;; but suggest in both
                etot = sqrt(t.e1^2 + t.e2^2)
                keep = where(t.sigmae*2 LT err_cut AND etot LT etot_cut AND $
                                t.r2[0] GT 0 AND t.r2[1] GT 0, nkeep)

                pstr = ntostr( float(nkeepold-nkeep)/ncat*100., 5, /round)+'%'
                print,'    Threw out '+ntostr(nkeepold-nkeep)+' in err/etot/r2 ('+pstr+')'
                IF nkeep GT 0 THEN BEGIN 
                    nkeepold = nkeep
                    t=t[keep]
                    outstruct=outstruct[keep]

                    print
                    print,'Finally kept '+ntostr(nkeep)

                    print,'Stuffing database'
                    self->struct2table, outstruct, table, $
                                    conn='user=postgres', $
                                    primary_key='photoid', status=stuff_status,$
                                    tmpdir=tmpdir
          
                    if stuff_status ne 0 then begin 
                        message,'could not stuff'
                    endif 

                    delvarx, t, outstruct
                endif else print,'-->none passed err/etot'
            endif else print,'-->none passed mag cuts'
        endif else print,'-->none passed pixel masks'

    endfor 

    self->princeton_create_indexes
    self->postgres::create_metatable, table, conn='user=postgres' 
    self->postgres::query, 'GRANT SELECT ON scat_princeton TO boss', conn='user=postgres'
    self->postgres::query, 'GRANT SELECT ON scat_princeton TO esheldon', conn='user=postgres'
    self->postgres::query, 'GRANT SELECT ON scat_princeton_meta TO boss', conn='user=postgres'
    self->postgres::query, 'GRANT SELECT ON scat_princeton_meta TO esheldon', conn='user=postgres'

end 


pro make_scat::princeton_create_indexes

    table='scat_princeton'
    columns = $
        ['stripe',$
         'run', $
         'cmodelmag[3]',$
         'e1eq','e2eq',$
         'e1pix','e2pix',$
         'e1survey','e2survey',$
         'r',$
         'r_r', $
         'r_i', $
         'seeing',$
         'maskflags']

    conn = 'user=postgres'
    pg=obj_new('postgres')
    self->postgres::create_index, table, columns, conn=conn

    self->postgres::query, 'analyze '+table, conn=conn, status=astatus
    if astatus ne self->postgres::status_val('no_result') then begin
        message,'analyze failed'
    endif

end 

FUNCTION make_scat::princeton_stripes_use, nstripe
  dir = self->princeton_bystripe_dir()
  freg = concat_dir(dir,'stripe[0-9][0-9].st')
  files = file_search(freg, count=nstripe)
  FOR i=0L, nstripe-1 DO BEGIN 
      hdr = read_idlheader(files[i])
      add_arrval, hdr.stripe, stripes
  ENDFOR 
  return,stripes
END 







function make_scat::source_sample_string, sample
  if n_elements(sample) eq 0 then $
    message,'You must send a sample number'
  return,'sample'+strn(sample, len=2, padchar='0')
end 

function make_scat::princeton_photoz_match_file, pztype
    if n_elements(pztype) eq 0 then begin
        message,'f=ms->princeton_photoz_match_file(pztype)'
    endif
    dir = self->princeton_source_dir()
    file = concat_dir(dir,'princeton-scat-pzmatch-'+pztype+'.st')
    return,file
end
function make_scat::princeton_photoz_match_read, pztype, hdr=hdr, status=status
    if n_elements(pztype) eq 0 then begin
        message,'st=ms->princeton_photoz_match_read(pztype)'
    endif
    file = self->princeton_photoz_match_file(pztype)
    print,'Reading file: ',file
    return, read_idlstruct(file, hdr=hdr, status=status)
end
; maybe should have picked zerr < 0.25
function make_scat::princeton_photoz_clause
    zmin = '0.01'
    zmax = '1.0'
    errmax = '0.6'
    clause='photoz_z between '+zmin+' and '+zmax+' and photoz_zerr < '+errmax
    return, clause
end

function make_scat::photoz_table, pztype
    if n_elements(pztype) eq 0 then begin
        message,'t=ms->photoz_table(pztype)'
    endif
    case pztype of
        'nn': return, 'zphot'
        'dr6cc2': return, 'zphotcc2'
        else: message,'dont yet support photoz type '+pztype
    endcase
end

; Use big memory machine
pro make_scat::princeton_photoz_match, pztype

    if n_elements(pztype) eq 0 then begin
        message,'ms->princeton_photoz_match, pztype'
    endif
    tm=systime(1)

    pg = obj_new('postgres')

    outfile = self->princeton_photoz_match_file(pztype)
    print
    print,'Will write to file: ',outfile

    ; Grab all the positions from the catalog
    squery = 'select photoid, clambda, ceta from scat_princeton'
    print
    print,squery
    scat=pg->query(squery)
    help, scat

    ; get positions and photozs from zphot table
    zclause=self->princeton_photoz_clause()
    ztable = self->photoz_table(pztype)
    zquery = $
        'select ra, dec, photoz_z, photoz_zerr from '+ztable+' where '+zclause
    print
    print,zquery
    zcat=pg->query(zquery)
    help, zcat

    obj_destroy, pg


    ;; convert positions from scat to ra/dec
    print,'converting to ra/dec'
    photoid=scat.photoid
    sclambda = scat.clambda
    sceta = scat.ceta
    nscat=n_elements(scat)
    scat=0
    csurvey2eq, sclambda, sceta, sra, sdec
    sclambda=0
    sceta=0

    ; this will allow us to micromanage the memory
    print,'Copying out data from zcat'
    zra = zcat.ra
    zdec = zcat.dec
    photoz = zcat.photoz_z
    photozerr = zcat.photoz_zerr
    zcat = 0

    ;; 1 arcsec match
    srad_arcsec = 1d
    print,'matching within ',srad_arcsec,' arcsec'
    srad = srad_arcsec/3600d*!dpi/180d
    htm_match, sra, sdec, zra, zdec, srad, smatch, zmatch, $
        maxmatch=1

    delvarx, zra, zdec, sra, sdec

    outstruct = { $
        photoid: 0LL, $
        photoz_z: 0.0, $
        photoz_zerr: 0.0 }

    nmatch = n_elements(smatch)
    print,'Matched '+ntostr(nmatch)+'/'+ntostr(nscat)
    outstruct = replicate(outstruct, nmatch)

    photoid = photoid[smatch]
    photoz=photoz[zmatch]
    photozerr=photozerr[zmatch]

    outstruct.photoid = photoid
    outstruct.photoz_z = photoz
    outstruct.photoz_zerr = photozerr

    hdr={zclause:zclause}

    delvarx, smatch, zmatch, photoid, photoz, photozerr
    print
    print,'Writing to file: ',outfile

    write_idlstruct, outstruct, outfile

    ptime,systime(1)-tm
end


function make_scat::princeton_sample_stuff_structdef, num
    st = {                      $
            stripe: 0,          $
            e1: 0.0,            $
            e2: 0.0,            $
            e1orig: 0.0,        $
            e2orig: 0.0,        $
            e1e1err: 0.0,       $
            e1e2err: 0.0,       $
            e2e2err: 0.0,       $
            phi_offset:0.0,     $
            ra: 0d,             $
            dec: 0d,            $
            clambda: 0d,        $
            ceta: 0d,           $
            photoz_z: 0.0,      $
            photoz_zerr: 0.0,   $
            modelmag: fltarr(5), $
            cmodelmag: fltarr(5), $
            htm_index: 0L       $
        }
        
    if n_elements(num) ne 0 then st=replicate(st, num)
    return, st
end

function make_scat::princeton_where_clause, sample

  if n_elements(sample) eq 0 then $
    message,'You must send a sample number'

  case sample of
      1: clause=''
      2: clause='( r >  0.585 )'  ; well resolved
      3: clause='( r <= 0.585 )'  ; less well resolved
      4: clause='( seeing < 1.5 )'
      ;; match flags on blended 8.  This should be very conservative
      5: clause='( (objc_flags & 8) = 0 )'
      6: clause='( r > 0.33333 )'
      7: clause='( r > 0.66666 )' ; this is actually the 2/3 point
      8: clause='( r > 0.4 )'
      9: clause='( r > 0.33333 )' ; here we use the new photozs
      else: message,'Unsupported sample '+ntostr(sample)
  endcase 
  return,clause

end 

; now requires large memory machine
; This reads from the main scat_princeton table and takes subsets.  Note,
; the photoz_match is done on the main scat_princeton table so it is correct
; to use it

function make_scat::pztype, sample
    if sample eq 9 then pztype = 'dr6cc2' else pztype='nn'
    return, pztype
end

pro make_scat::princeton_sub_sample, sample

    ; create a sub-sample
    ; see princeton_where_clause.  6 is r > 1/3

    tm=systime(1)
    if n_elements(sample) eq 0 then begin
        print,'-syntax: ms->princeton_stuff_sample, sample'
        on_error, 2
        message,'Halting'
    endif

    tmpdir = self->princeton_source_dir()

    sstr = ntostr(sample)
    table = 'scat_princeton'+sstr

    where_clause = self->princeton_where_clause(sample)

    query = $
        'SELECT '+$
            'photoid, '+$
            'stripe, '+$
            'e1survey, e2survey, e1pix, e2pix, e1eq, e2eq, '+$
            'e1e1err, e1e2err, e2e2err, '+$
            'clambda, ceta, '+$
            'htm_index '+$
        'FROM '+$
            'scat_princeton '+$
        'WHERE '+where_clause

    print,query
    tmpcat = self->query(query)
    nstart=n_elements(tmpcat)

    pztype=self->pztype(sample)
    zmatchcat = self->princeton_photoz_match_read(pztype)
    
    print,'Extracting matches'
    match, tmpcat.photoid, zmatchcat.photoid, mt, mz, /sort
    nmatch=n_elements(mt)
    print,'Keeping '+ntostr(nmatch)+'/'+ntostr(nstart)

    tmpcat = tmpcat[mt]
    zmatchcat = zmatchcat[mz]

    scat = self->princeton_sample_stuff_structdef(nmatch)
    struct_assign, tmpcat, scat, /verbose, /nozero
    scat.photoz_z = zmatchcat.photoz_z
    scat.photoz_zerr = zmatchcat.photoz_zerr

    delvarx, zmatchcat, tmpcat

    print,'Converting to ra/dec'
    csurvey2eq, scat.clambda, scat.ceta, ra, dec
    scat.ra = ra
    scat.dec = dec

    ; database operations
    nores = self->postgres::status_val('no_result')
    conn = 'user=postgres'
    self->postgres::struct2table, scat, table, conn=conn, status=stuff_status, tmpdir=tmpdir
    if stuff_status ne 0 then message,'struct2table on '+table+' failed'

    columns = ['stripe','photoz_z', 'photoz_zerr']
    self->postgres::create_index, table, columns, conn=conn

    query = 'analyze '+table
    print,query
    self->query, query, conn=conn, status=astatus
    if astatus ne nores then message,'analyze failed'

    ; Make sure the sdss user can select
    query = 'GRANT select ON '+table+' TO sdss'
    print
    print,query
    self->postgres::query,  query, conn=conn, status=gstatus
    if gstatus ne nores then message,'Could not grant select to sdss'

    ptime, systime(1)-tm

end 



; this version gets written to the file
function make_scat::princeton_source_structdef, num, equatorial=equatorial

    if not keyword_set(equatorial) then begin
        st = {                  $
            stripe: 0,          $
            e1: 0.0,            $
            e2: 0.0,            $
            e1e1err: 0.0,       $
            e1e2err: 0.0,       $
            e2e2err: 0.0,       $
            clambda: 0d,        $
            ceta: 0d,           $
            photoz_z: 0.0,      $
            photoz_zerr: 0.0,   $
            htm_index: 0L }
    endif else begin
        st = {                  $
            stripe: 0,          $
            e1: 0.0,            $
            e2: 0.0,            $
            e1e1err: 0.0,       $
            e1e2err: 0.0,       $
            e2e2err: 0.0,       $
            ra: 0d,             $
            dec: 0d,            $
            photoz_z: 0.0,      $
            photoz_zerr: 0.0,   $
            htm_index: 0L }
    endelse

    if n_elements(num) ne 0 then st=replicate(st, num)
    return, st
end


; now requires large memory machine for sure
pro make_scat::princeton_scat_write, sample, equatorial=equatorial, pix=pix

    if n_elements(sample) eq 0 then begin
        print,'-syntax: ms->princeton_scat_write, sample, /equatorial'
        on_error, 2
        message,'Halting'
    endif


    outfile = self->princeton_source_file(sample, equatorial=equatorial)
    sstr = ntostr(sample)
    table = 'scat_princeton'+sstr

    if keyword_set(equatorial) then begin
        query = $
            'select '+$
               'stripe, e1eq as e1, e2eq as e2, e1e1err, e1e2err, e2e2err, '+$
               'ra, dec, photoz_z, photoz_zerr, htm_index '+$
            'from '+$
               table
    endif else if keyword_set(pix) then begin
        query = $
            'select '+$
               'stripe, e1pix as e1, e2pix as e2, e1e1err, e1e2err, e2e2err, '+$
               'ra, dec, photoz_z, photoz_zerr, htm_index '+$
            'from '+$
               table
    endif else begin
        query = $
            'select '+$
               'stripe, e1survey as e1, e2survey as e2, e1e1err, e1e2err, e2e2err, '+$
               'clambda, ceta, photoz_z, photoz_zerr, htm_index '+$
            'from '+$
               table
    endelse 
           
    print,query
    scat = self->postgres::query(query)

    ; remove the -9999 default values
    w=where(scat.e1 gt -10 and scat.e2 gt -10, nw)
    if nw eq 0 then message,'No good e1 > -10'
    scat = scat[w]

    stripes = scat[ uniq(scat.stripe,sort(scat.stripe)) ].stripe
    where_clause = self->princeton_where_clause(sample)
    zclause=self->princeton_photoz_clause()
    hdr = { $
        stripes: stripes, $
        sample: sample, $
        clause: where_clause+':'+zclause }

    print
    print,'Writing to '+ntostr(n_elements(scat))+' objects to file: ',outfile
    if not keyword_set(equatorial) then begin
        write_idlstruct, scat, outfile, hdr=hdr
    endif else begin
        mwrfits2, scat, outfile, /create, /destroy
    endelse


end 



; create a sample for the weightin/histogram matching techniqe
pro make_scat::princeton_weight_input, number

    ; select randomly from the first 30 million
    ntot = '30000000'
    nrand = '1000000'

    dir = self->princeton_source_dir()
    dir = concat_dir(dir, 'weight_input')
    if not fexist(dir) then file_mkdir, dir

    nstr=ntostr(number,format='(I02)')
    file = 'princeton-scat-weight-input'+nstr+'.st'
    file = concat_dir(dir, file)
    print,'Will write to: ',file

    query = 'SELECT '+$
                'ra, dec, modelmag, photoz_z, photoz_zerr '+$
            'FROM scat_princeton6 '+$
            'ORDER BY random() LIMIT '+nrand
    print,query
    struct=self->postgres::query(query)
    print
    print,'Writing to file: ',file
    write_idlstruct, struct, file, /ascii

end









;; htm for princeton
FUNCTION make_scat::princeton_htmrev_file, sample
  file = self->princeton_source_file(sample)
  file = repstr(file,'.st', '-htmrev.bin')
  return,file
END 

PRO make_scat::princeton_make_htm_revind, sample

    rev_file = self->princeton_htmrev_file(sample)
  
    print
    print,'Will write to file: ',rev_file


    sfile = self->princeton_source_file(sample)
    print
    print,'Reading file: ',sfile
    st = read_idlstruct(sfile, columns='htm_index')

    htm_index = st.htm_index
    st = 0

    print
    print,'Creating reverse indices'

    minid = min(htm_index, max=maxid)
    h = histogram( htm_index-minid, min=0, rev=rev )
    nrev = n_elements(rev)

    rev_file = self->princeton_htmrev_file(sample)

    print
    print,'Writing to rev file: ',rev_file
    openw, lun, rev_file, /get_lun
    writeu, lun, nrev

    writeu, lun, minid
    writeu, lun, maxid

    writeu, lun, rev

    free_lun, lun

    rev = 0


END 






PRO make_scat::princeton_number_counts, sample, scat=scat, arcmin=arcmin, $
             stripelow=stripelow, stripehigh=stripehigh, $
             dops=dops, dopng=dopng

  if n_elements(sample) eq 0 then begin
      print,'-Syntax: ms->princeton_number_counts, sample, scat=, /arcmin, /stripelow, /stripehigh, /dops, /dopng'
      return
  endif

  plotdir = '~/plots/make_scat'
  sstr=ntostr(sample, f='(I0)')
  psfile = concat_dir(plotdir,'princeton-scat-number-counts-'+sstr+'.eps')
  pngfile = concat_dir(plotdir,'princeton-scat-number-counts-'+sstr+'.png')
  IF keyword_set(dops) THEN BEGIN 
      begplot, name=psfile, /encap, /color
      ecolor = !grey40
      ethick = !p.thick
  ENDIF ELSE BEGIN 
      ecolor = !yellow
      ethick = !p.thick
  ENDELSE 
  IF n_elements(scat) EQ 0 THEN BEGIN 
      table = 'scat_princeton'+ntostr(sample)
      query = $
        'SELECT '+$
        '   stripe, cmodelmag[2] as cmodelr, e1e1err '+$
        'FROM '+table

      print,query
	  pg=obj_new('postgres')
      scat = pg->query(query)
  ENDIF 

  area = 6325d                  ; square degrees

  IF keyword_set(arcmin) THEN BEGIN 
      area = area*3600d         ; square arcminutes
      ytitle = 'Density [#/arcmin!U2!N]'
  ENDIF ELSE BEGIN 
      ytitle = 'Density [#/degree!U2!N]'
  ENDELSE 
  num = n_elements(scat)

  IF keyword_set(stripelow) THEN BEGIN 
      w = where(scat.stripe LT 25, ngood)
  ENDIF ELSE IF keyword_set(stripehigh) THEN BEGIN 
      w = where(scat.stripe GE 25, ngood)
  ENDIF ELSE BEGIN 
      w=lindgen(n_elements(scat))
  ENDELSE 

  weight = 1.0/(0.32^2 + scat.e1e1err^2)

  binsize = 0.1
  bs=binner(scat[w].cmodelr,weight[w],binsize=binsize, min=18, max=22)

  wh = where(bs.hist NE 0 AND bs.xmean NE 0)
  xb = bs.xmean[wh]
  yb = bs.ymean[wh]
  yberr = bs.yerr[wh]
  hist = bs.hist[wh]

  colprint,xb,yb,hist

  density = float(hist)/area
  wdensity = yb/max(yb)*density

  wline = 2


  xtitle = 'r [mag]'
  ytitle='Number [Arbitrary Units]'

  aplot, !gratio, xb, density, psym=10, xtitle=xtitle, ytitle=ytitle
  oplot, xb, wdensity, psym=10, color=ecolor, thick=ethick
  oplot, xb, yb/max(yb)*max(density), line=wline

  legend,['Number', 'Effective Number', 'Relative Weight'], $
    line=[0,0,wline], color=[!p.color, ecolor, !p.color], $
    /right, box=0, charsize=1, thick=[!p.thick, ethick, !p.thick]

  IF keyword_set(dops) THEN BEGIN 
      endplot, /trim_bbox
  ENDIF ELSE IF keyword_set(dopng) THEN BEGIN 
      write_png, pngfile, tvrd(/true)
  ENDIF 

END 



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Intrinsic alignment catalogs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION make_scat::intrinsic_scat_extract, sample, query=query

  CASE sample OF
      1: BEGIN 
          ;; not selecting htm_index here because it is 64-bit in
          ;; the database
          query = $
            'SELECT '+$
            'stripe, '+$
            'm_e1_corr_h[2] as e1, m_e2_corr_h[2] as e2, m_r_h[2] as r, '+$
            'm_e1e1err[2] as e1e1err, m_e1e2err[2] as e1e2err, m_e2e2err[2] as e2e2err, '+$
            'm_e1_psf as e1psf, m_e2_psf as e2psf, '+$
            'clambda, ceta, '+$
            'z, z_err '+$
            'FROM specgal_lss'
          
          print,query
		  pg=obj_new('postgres')
          t = pg->query(query)
          
          outstruct = {stripe: 0, $
                       e1: 0.0, e2:0.0, $
                       e1e1err:0.0, e1e2err:0.0, e2e2err:0.0, $
                       clambda: 0d, ceta:0d, $
                       z:0.0, z_err:0.0, $
                       htm_index: 0L}
          
          print,'Making resolution cuts'
          keep = where(t.e1 GE -1 and $
                       t.r GT 0 AND $
                       t.r LT 0.7, nkeep)
          
          t = t[keep]

          ;; correct for slope
          fitstruct = egal_vs_epsf_read(2)
          
          correct_eslope, fitStruct, $
            t.r, t.e1, t.e2, $
            t.e1psf, t.e2psf, $
            e1c, e2c


          ;; Note hirata stuff already dilution corrected, just need
          ;; to scale errors
          corr = 1.0/(1.0-t.r)
          t.e1 = e1c
          t.e2 = e2c
          t.e1e1err = t.e1e1err*corr
          t.e1e2err = t.e1e2err*corr
          t.e2e2err = t.e2e2err*corr  
          
          outstruct = replicate(outstruct, nkeep)
          copy_struct, t, outstruct
          t=0          
          csurvey2eq, outstruct.clambda, outstruct.ceta, ra, dec
          
          depth = 10
          outstruct.htm_index = htm_index(ra,dec,depth)
          ra = 0
          dec = 0

          return,outstruct
      END 
      ELSE: message,'Unsupported intrinsic align sample '+strn(sample)
  ENDCASE 
END 

FUNCTION make_scat::intrinsic_source_dir
  return,concat_dir(esheldon_config('lensinput_dir'), 'srcgal/intrinsic_align')
END 
FUNCTION make_scat::intrinsic_source_file, sample
  IF n_elements(sample) EQ 0 THEN BEGIN 
      message,'You must enter sample number'
  ENDIF 
  CASE sample OF 
      1: BEGIN 
          ;; This is actually the intrinsic alignment
          ;; sources from sample 14
          dir = self->intrinsic_source_dir()
          file = concat_dir(dir,'sample14.st')
          return,file
      END 
      ELSE: message,'Unknown source sample: '+strn(sample)
  ENDCASE

END 

PRO make_scat::intrinsic_source_write, struct, sample, append=append, hdr=hdr
  file = self->intrinsic_source_file(sample)
  print
  print,'Writing source file: ',file
  write_idlstruct, struct, file, append=append, hdr=hdr
END 
FUNCTION make_scat::intrinsic_source_get, sample
  file = self->intrinsic_source_file(sample)
  print
  print,'Reading source file: ',file
  struct = read_idlstruct(file)
  return,struct
END 

PRO make_scat::intrinsic_scat_write, sample

  outfile = self->intrinsic_source_file(sample)

  print
  print,'Will write to file: ',outfile

  struct = self->intrinsic_scat_extract(sample, query=query)

  hdr = {stripes: -1, $
         sample: sample, $
         clause: query}

  self->intrinsic_source_write, struct, sample, hdr=hdr
  
END 

FUNCTION make_scat::intrinsic_htmrev_file, sample
  file = self->intrinsic_source_file(sample)
  file = repstr(file, '.st', '-htmrev.bin')
  return,file
END 

PRO make_scat::intrinsic_make_htm_revind, sample

  st = self->intrinsic_source_get(sample)
  htm_index = st.htm_index
  st = 0

  print
  print,'Creating reverse indices'
  
  minid = min(htm_index, max=maxid)
  h = histogram( htm_index-minid, min=0, rev=rev )
  nrev = n_elements(rev)
  
  ;; sample 1
  rev_file = self->intrinsic_htmrev_file(sample)
  
  print
  print,'Writing to rev file: ',rev_file
  openw, lun, rev_file, /get_lun
  writeu, lun, nrev
  
  writeu, lun, minid
  writeu, lun, maxid
  
  writeu, lun, rev
  
  free_lun, lun
  
  rev = 0


END 
























;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Pixel stuff
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


FUNCTION make_scat::pixel_struct
  struct = {htm_index: 0L, $
            nobj: 0, $
            e1: 0.0, $
            e2: 0.0, $
            e1e1err: 0.0, $
            e1e2err: 0.0, $
            e2e2err: 0.0, $
            aeta1: 0.0, $
            aeta2: 0.0, $
            clambda: 0d, $
            ceta: 0d}

  return,struct
END 

FUNCTION make_scat::pixel_file, region
  dir = self->byregion_dir()
  file = concat_dir(dir,'scat-pixel-region'+ntostr(long(region))+'.bin')
  return,file
END 
FUNCTION make_scat::pixel_htmrev_file, region
  file = self->pixel_file(region)
  file = repstr(file, '.bin','-htmrev.bin')
  return,file
END 
FUNCTION make_scat::read_pixel, region, npix

  file = self->pixel_file(region)
  openr, lun, file, /get_lun
  
  npix = 0L
  readu, lun, npix
  print
  print,'Reading '+ntostr(npix)+' from file: ',file

  struct = self->pixel_struct()
  struct = replicate(struct, npix)

  readu, lun, struct

  free_lun, lun

  return,struct

END 


PRO make_scat::pixelize_region, region, st=st

  IF n_elements(region) EQ 0 THEN BEGIN
      print,'-Syntax: obj->pixelize_region, region'
      return
  ENDIF 

  IF n_elements(st) EQ 0 THEN st = self->read_byregion_file(region)
  outfile = self->pixel_file(region)
  htmfile = self->pixel_htmrev_file(region)

  print
  print,'outfile = ',outfile
  print,'htmfile = ',htmfile

  print
  print,'Histogramming htm_index'
  hist = histogram(st.htm_index, min=0, rev=rev)

  wh = where(hist NE 0, nhist)
  help,wh

  struct = self->pixel_struct()

  struct = replicate(struct, nhist)

  omega_m = 0.27
  aeta_zero = aeta_lambda(0.0, omega_m)
  FOR i=0L, nhist-1 DO BEGIN 

      bin = wh[i]
      
      w = rev[ rev[bin]:rev[bin+1]-1 ]

      nobj = n_elements(w)

      struct[i].nobj = nobj

      ;; mean shapes
      shear_err = sqrt( 0.32^2 + st[w].e1e1err^2 )
      wmom, st[w].e1, shear_err, mean_e1, sig_e1, err_e1
      shear_err = sqrt( 0.32^2 + st[w].e2e2err^2 )
      wmom, st[w].e2, shear_err, mean_e2, sig_e2, err_e2

      ;; get aeta's and average.  Not using photoz errors yet

      aeta_s = aeta_lambda(st[w].photoz_z, omega_m)

      aeta1 = 1.0/(aeta_zero - aeta_s)
      aeta2 = aeta_s*aeta1

      err = sqrt( 0.32^2 + (st[w].e1e1err^2 + st[w].e2e2err^2)/2.0 )
      wmom, aeta1, err, mean_aeta1, sig, err
      wmom, aeta2, err, mean_aeta2, sig, err

      ;; mean positions
      wmom, st[w].clambda, err, mean_clambda, sig, err
      wmom, st[w].ceta, err, mean_ceta, sig, err

      struct[i].htm_index = st[w[0]].htm_index
      struct[i].e1 = mean_e1
      struct[i].e2 = mean_e2
      struct[i].e1e1err = err_e1
      struct[i].e1e2err = 0.0
      struct[i].e2e2err = err_e2
      struct[i].aeta1 = mean_aeta1
      struct[i].aeta2 = mean_aeta2
      struct[i].clambda = mean_clambda
      struct[i].ceta = mean_ceta

      IF ((i+1) MOD 10000) EQ 0 THEN BEGIN 
          print,ntostr(i+1)+'/'+ntostr(nhist)+'  nobj = '+ntostr(nobj)
      ENDIF 

  ENDFOR 

  ;; free memory
  delvarx, st

  ;; Will be sorted by htm_index
  min_htm_index = min(struct.htm_index, max=max_htm_index)

  print
  print,'Writing to file: ',outfile
  openw, lun, outfile, /get_lun
  writeu, lun, nhist
  writeu, lun, struct

  free_lun, lun




  ;; htm reverse indices for faster lookup
  print
  print,'Histogramming pixels'
  hist = histogram(struct.htm_index-min_htm_index, rev=rev)
  nrev = n_elements(rev)

  print,'Writing to htmrev file: ',htmfile

  openw, lun, htmfile, /get_lun

  writeu, lun, nrev
  writeu, lun, min_htm_index
  writeu, lun, max_htm_index
  writeu, lun, rev

  free_lun, lun





  print,'Done'

END 





PRO make_scat::free

  ;; Don't free parstruct
  ptr_free, $
    self.cat, $
    self.keep, $
    self.band_keep, $
    self.final_keep, $
    self.output_cat

END 






FUNCTION make_scat::cleanup

  self->free
  ptr_free, self.parstruct
  return,1
END 

PRO make_scat__define

  struct = {make_scat, $
            parstruct: ptr_new(), $
            $
            cat: ptr_new(), $
            $
            keep: ptr_new(), $
            nkeep: 0LL, $
            band_keep: ptrarr(5), $
            band_nkeep: lon64arr(5), $
            final_keep: ptr_new(), $
            final_nkeep: 0LL, $
            $
            output_cat: ptr_new(), $
            status: 0, $
            write_status: 0, $
            INHERITS postgres $
           }

END 
