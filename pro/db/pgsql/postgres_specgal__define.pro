;; Stuff Blanton's spec files into a postgres database

FUNCTION postgres_specgal::init, stripe=stripe

  self.tablename = 'vagc'
  self.stripe = -1

  ;; To avoid conflict with reset in inherited classes
  self->postgres_specgal::reset, stripe=stripe

  return,1
END 

PRO postgres_specgal::reset, stripe=stripe

  IF n_elements(stripe) NE 0 THEN BEGIN 
      self.stripe = stripe
  ENDIF 

  self.write_status = -1
  self.stuff_status = -1

END 

FUNCTION postgres_specgal::stripe
  return,self.stripe
END 
FUNCTION postgres_specgal::write_status
  return,self.write_status
END 
FUNCTION postgres_specgal::stuff_status
  return,self.stuff_status
END 

FUNCTION postgres_specgal::structdef

  ;; Read a struct
  name = spectra_name(9)
  struct=mrdfits(name, 1, rows=0)

  ;; add match info and photoid
  struct = create_struct('spectroid', 0LL,$
                         'photoid_p', 0LL, $
                         'stripe', 0, $
                         struct, $
                         'match_photoid', -1LL, $
                         'completeness', -1.0, $
                         'poly_id', -1L, $
                         'poly_area', -1d, $
                         'htm_index', 0LL)

  return,struct

END 

PRO postgres_specgal::tabledef, sqlfile

  IF n_elements(sqlfile) EQ 0 THEN BEGIN 
      lun = -1
  ENDIF ELSE BEGIN 
      openw,lun,sqlfile,/get_lun
  ENDELSE 

  struct = self->structdef()
  coldefs = self->struct2coldefs(struct)

  tags = tag_names(struct)
  w=where(tags EQ 'RERUN')
  coldefs[w[0]] = 'rerun SMALLINT NOT NULL'
  w=where(tags EQ 'CAMCOL')
  coldefs[w[0]] = 'camcol SMALLINT NOT NULL'
  w=where(tags EQ 'FIELD')
  coldefs[w[0]] = 'field SMALLINT NOT NULL'

  w=where(tags EQ 'OBJC_TYPE')
  coldefs[w[0]] = 'objc_type SMALLINT NOT NULL'
  w=where(tags EQ 'TYPE')
  coldefs[w[0]] = 'type SMALLINT[5] NOT NULL'
  

  ;; Add the primary key
  coldefs = [coldefs, 'PRIMARY KEY (spectroid)']

  ncoldefs = n_elements(coldefs)
  printf,lun,'CREATE TABLE specgal'
  printf,lun,'('
  FOR i=0L, ncoldefs-2 DO BEGIN 
      printf,lun,coldefs[i]+', '
  ENDFOR 
  printf,lun,coldefs[i]
  printf,lun,');'

  IF n_elements(sqlfile) NE 0 THEN free_lun,lun

END 

FUNCTION postgres_specgal::input_file

  dir = '/net/cheops1/data1/esheldon/tmp/'
  file = dir + 'spec_stripe'+stripe2string(self.stripe)+'_stuff.pgsql'
  return,file

END 

PRO postgres_specgal::write_input

  htm_depth = 17

  tm = systime(1)

  stripe = self.stripe
  IF stripe EQ -1 THEN BEGIN  
      message,'You must initialize stripe',/inf
      message,'obj->reset, stripe=stripe'
  ENDIF 

  tmpfile = self->input_file()
  print
  print,'Will write to temporary file: ',tmpfile

  lss_columns = ['run','rerun','camcol','field','id',$
                 'completeness', 'poly_id', 'poly_area']
  get_spectra_lcat, stripe, lss, /lss, columns=lss_columns
  get_spectra_lcat, stripe, vagc


  ;; Create the input struct
  struct =  self->structdef()
  nobj = n_elements(vagc)

  struct = replicate(struct, nobj)

  ;; Copy in from vagc
  copy_struct, vagc, struct

  ;; Copy in from lss.  
  sphoto_match, vagc, lss, mvagc, mlss
  struct[mvagc].completeness = lss[mlss].completeness
  struct[mvagc].poly_id = lss[mlss].poly_id
  struct[mvagc].poly_area = lss[mlss].poly_area

  ;; Set ids
  struct.spectroid = sdss_spectroid(struct)
  struct.photoid_p = sdss_photoid(struct)

  struct.htm_index = htm_index(struct.ra, struct.dec, htm_depth)

  w = where(struct.match_rerun GT 0, nw)

  IF nw NE 0 THEN BEGIN 
      struct[w].match_photoid = sdss_photoid(struct[w].run, $
                                             struct[w].match_rerun, $
                                             struct[w].camcol, $
                                             struct[w].field, $
                                             struct[w].match_id)
  ENDIF 

  struct.stripe = self.stripe

  ;; Write file
  print
  print,'Writing to file: ',tmpfile
  ascii_write, struct, tmpfile, /bracket_arrays, status=status

  ptime,systime(1)-tm

  IF status NE 0 THEN BEGIN 
      message,'Error writing input file',/inf
      return
  ENDIF 

  self.write_status = 0

END 

PRO postgres_specgal::stuff

  tm = systime(1)

  IF self.write_status EQ 1 THEN BEGIN 
      message,$
        'self.write_status indicates an error: '+ntostr(self.write_status),/inf
      return
  ENDIF 

  file = self->input_file()
  self.connect_info = 'user=postgres'
  query = "COPY specgal FROM '"+file+"'"
  print,query
  self->query, query

  IF self.query_status NE self->status_val('no_result') THEN BEGIN
      self.stuff_status = 1
  ENDIF ELSE BEGIN 
      self.stuff_status = 0
  ENDELSE 

  ptime,systime(1)-tm

END 

PRO postgres_specgal::delete_input
  print
  print,'Deleting file: ',self->input_file()
  file_delete, self->input_file(), /quiet
END 

FUNCTION postgres_specgal::stripes

  read_runlist, runstruct, /silent
  
  catname = vagc_catname(lss=lss, letter=letter, post=post, sample=sample)

  byrun_dir = sdssidl_config('spec_dir') + 'blanton/gal_collated/byrun_matched/'
  files = findfile(byrun_dir, count=nf)
  w=where( strmatch(files,'*'+catname+'*'), nmatch)
  IF nmatch EQ 0 THEN BEGIN 
      message,'No '+catname+' files found'
  ENDIF 

  nf = nmatch
  files = files[w]
  runs = lonarr(nf)

  FOR i=0L, nf-1 DO BEGIN 
      t=( strsplit(files[i],'run',/extract) )[0]
      runstr = ( strsplit(t,'-',/extract) )[0]
      runs[i] = long(runstr)
  ENDFOR 

  match, runs, runstruct.run, mruns, mstruct

  rmd = rem_dup(runstruct[mstruct].stripe)
  stripes = runstruct[mstruct[rmd]].stripe
  w=where(stripes LE 86 AND $
          stripes NE 61 AND $
          stripes NE 62, nstripe)
  stripes = stripes[w]

  return,stripes

END 

PRO postgres_specgal::stuff_all

  stripes = self->stripes()

  nstripe = n_elements(stripes)
  
  FOR i=0L, nstripe-1 DO BEGIN 

      self->reset, stripe=stripes[i]

      ;; Leave the inputs at first just for this test
      IF NOT fexist(self->input_file()) THEN BEGIN 
          self->write_input
          if self->write_status() ne 0 THEN message,'Write error'
          self->stuff
          IF self->stuff_status() NE 0 THEN message,'Stuff error'
      ENDIF ELSE BEGIN 
          print,'SKIPPING'
      ENDELSE 
  ENDFOR 

END 


PRO postgres_specgal__define

  struct = {$
             postgres_specgal, $
             stripe: 0, $
             tablename:'', $
             write_status: 0, $
             stuff_status: 0, $
             INHERITS postgres $
           }

END 
