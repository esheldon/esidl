FUNCTION postgres_adatc_extra::init, run, camcol, rerun=rerun

  nrun = n_elements(run)
  ncamcol = n_elements(camcol)

  IF nrun EQ 0 OR ncamcol EQ 0 THEN BEGIN 
      message,'You must initialize the run and camcol',/inf
      message,$
        "-Syntax: obj=obj_new('postgres_adatc_extra',run,camcol,rerun=')",/inf
      return,0
  ENDIF 

  tablename = 'adatc_extra'

  test = self->postgres_stuff_camcol::init(run,camcol,tablename, rerun=rerun)
  IF NOT test THEN BEGIN 
      message,'Could not initialize postgres_stuff_camcol',/inf
      return,0
  ENDIF 
  return,1

END 

PRO postgres_adatc_extra::reset, $
                  run=run, camcol=camcol, rerun=rerun

  self->postgres_stuff_camcol::reset, run=run, camcol=camcol, rerun=rerun

END 

FUNCTION postgres_adatc_extra::structdef, photo_struct=photo_struct

  val1 = -9999.
  val2 = -9999
  val3 = -9999L
  arrval1 = replicate(-9999., 5)
  arrval2 = replicate(-9999, 5)
  arrval3 = replicate(-9999L, 5)

  photo_struct = {$
                   parent: val2, $
                   nchild: val2, $
                   $
                   objc_flags: val3, $
                   objc_flags2: val3, $
                   objc_rowc: val1, $
                   objc_colc: val1, $
                   $
                   counts_model: arrval1, $
                   counts_modelerr: arrval1, $
                   reddening: arrval1, $
                   $
                   flags: arrval3, $
                   flags2: arrval3, $
                   $
                   status: val3, $
                   $
                   m_e1: arrval1, $
                   m_e2: arrval1, $
                   m_e1e1err: arrval1, $
                   m_e1e2err: arrval1, $
                   m_e2e2err: arrval1,  $
                   m_rr_cc: arrval1, $
                   m_rr_ccerr: arrval1, $
                   m_cr4: arrval1, $
                   m_e1_psf: arrval1, $
                   m_e2_psf: arrval1, $
                   m_rr_cc_psf: arrval1, $
                   m_cr4_psf: arrval1 $
                 }

  outst = create_struct('photoid', 0LL, photo_struct)

  photo_struct = $
    create_struct('run', val3, $
                  'rerun', val2, $
                  'camcol', val2, $
                  'field', val2, $
                  'id', val3, $
                  photo_struct)

  return,outst

END 

PRO postgres_adatc_extra::tabledef, sqlfile

  IF n_elements(sqlfile) EQ 0 THEN BEGIN 
      lun = -1
  ENDIF ELSE BEGIN 
      openw,lun,sqlfile,/get_lun
  ENDELSE 

  struct = self->structdef()
  tags = tag_names(struct)

  coldefs = self->struct2coldefs(struct)

  ;; Add the primary key
  coldefs = [coldefs, 'PRIMARY KEY (photoid)']

  ncoldefs = n_elements(coldefs)
  printf,lun,'CREATE TABLE adatc_extra'
  printf,lun,'('
  FOR i=0L, ncoldefs-2 DO BEGIN 
      printf,lun,coldefs[i]+', '
  ENDFOR 
  printf,lun,coldefs[i]
  printf,lun,');'

  IF n_elements(sqlfile) NE 0 THEN free_lun,lun


END 

PRO postgres_adatc_extra::write_input

  ;; status is 1 unless we reach end
  self.write_status = 1

  time=systime(1)
  
  sf = obj_new('sdss_files')

  run = self.run
  rerun = self.rerun
  camcol = self.camcol

  runstr = sf->run2string(run)
  cstr = ntostr(camcol)

  files = sf->filelist('tsObj',run,camcol,rerun=rerun,fields=fields)
  nfields = n_elements(fields)

  tmpfile = self->input_file()

  file_delete, tmpfile, /quiet

  print
  print,'Will write to temporary file: ',tmpfile

  ;; Get all photoid for this camcol.  
  query = $
    'SELECT photoid FROM adatc WHERE '+$
    'run = '+ntostr(run)+' AND camcol = '+ntostr(camcol)
  print,query

  astr = pgsql_query(query, status=status)

  self->query, query

  IF self.query_status NE self->status_val('success') THEN BEGIN 
      return
  ENDIF 

  nrows = self->nrows()
  astr = self->struct(/no_copy,/free)

  outst = self->structdef(photo_struct=photo_struct)
  taglist = tag_names(photo_struct)

  tsobj_struct = sdss_read('tsobj',run, camcol, rerun=rerun,$
                            /all, $
                            taglist=taglist, verbose=1, $
                            status=ts_status)

  IF n_elements(tag_names(tsobj_struct)) NE n_elements(taglist) THEN BEGIN 
      message,'Did not find all tags'
  ENDIF 

  IF ts_status EQ 0 THEN BEGIN 

      phid = photoid(tsobj_struct)
      match, phid, astr.photoid, mph, ma, /sort

      IF mph[0] EQ -1 THEN BEGIN 
          message,'No matches',/inf
          return
      ENDIF 

      nmatch = n_elements(mph)
      print,'Found '+ntostr(nmatch)+' matches'
      IF nmatch NE nrows THEN BEGIN 
          message,'Error: nmatch less than nrows',/inf
          return
      ENDIF ELSE BEGIN 
          print,'All matched'
      ENDELSE 

      outStruct = replicate(outst, nmatch)
      tsobj_struct = tsobj_struct[mph]
      phid = phid[mph]

      copy_struct, tsobj_struct, outStruct
      outStruct.photoid = phid

      tsobj_struct = 0

      print
      print,'Writing to file: ',tmpfile

      ascii_write, outStruct,  tmpfile, /bracket_arrays, status=status
;      write_idlstruct, outStruct, tmpfile, $
;        /ascii, delim='tab', /noheader, $
;        error=error

      outstruct = 0

      IF status NE 0 THEN BEGIN 
          message,'Error writing input file',/inf
          return
      ENDIF 

      self.write_status = 0
      
  ENDIF ELSE BEGIN 
      message,'Error reading from tsobj',/inf
      return
  ENDELSE 

END 




; What's the purpose of this?
PRO postgres_adatc_extra::write_newtable_input

  run = self.run
  rerun = self.rerun
  camcol = self.camcol

  minphotoid = photoid(run,rerun,camcol,0,0)
  maxphotoid = photoid(run,rerun,camcol,2000,5000)
  minstr = ntostr(minphotoid)
  maxstr = ntostr(maxphotoid)
  query = 'SELECT photoid FROM adatc WHERE photoid '+$
    'between '+minstr+' AND '+maxstr


  query = 'SELECT '+$
    'a.photoid,a.run,a.rerun,a.camcol,a.field,a.id,a.value_flags,'+$
    'a.corrselect_flags, a.objc_prob_gal, a.objc_prob_flags, '+$
    'a.cmodel_counts, a.cmodel_countserr, a.cmodel_counts_ext,'+$
    'a.m_e1_corr, a.m_e2_corr, a.m_r, '+$
    'a.m_e1_corr_h, a.m_e2_corr_h, a.m_r_h,'+$
    'ae.m_e1e1err, ae.m_e1e2err, ae.m_e2e2err,'+$
    'a.compea4flags, '+$
    'a.seeing, '+$
    'a.rotation, '+$ 
    'a.clambda, a.ceta, '+$
    'a.htm_index '+$
    'FROM adatc AS a, adatc_extra AS ae '+$
    'WHERE a.photoid BETWEEN '+minstr+' AND '+maxstr+' '+$
    'AND a.photoid = ae.photoid '


  input_file = self->input_file()
  print,'Writing to file: ', input_file
  print,query

  self->query, query, file=input_file
  IF self.query_status NE self->status_val('no_result') THEN BEGIN 
      return
  ENDIF 

  self.write_status = 0

END 

PRO postgres_adatc_extra__define

  struct = {$
             postgres_adatc_extra, $
             INHERITS postgres_stuff_camcol $
           }

END 
