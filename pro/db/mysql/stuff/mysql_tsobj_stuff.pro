PRO mysql_tsobj_stuff_files, run, camcol, files, fields, $
                             rerun=rerun, $
                             start=start, nfields=nfields

  files = sdss_filelist('tsObj', run, $
                        rerun=rerun, camcol=camcol, $
                        fields=fields)


  ;; Allow the user to choose a subset of the fields
  IF n_elements(start) NE 0 THEN BEGIN 
      w=where(fields GE start, nw)
      IF nw EQ 0 THEN BEGIN 
          message,'No fields <= start field = '+ntostr(start),/inf
          return
      ENDIF 
      fields = fields[w]
      files = files[w]
  ENDIF ELSE BEGIN 
      start = fields[0]
  ENDELSE 

  IF n_elements(nfields) NE 0 THEN BEGIN 
      w=where(fields LE (start + nfields - 1), nw)
      fields = fields[w]
      files = files[w]
  ENDIF 

  nf = n_elements(fields)
  ff = fields[0]
  lf = fields[nf-1]
  print,'Processing '
  print,'   Run: '+run2string(run)+$
    '  Rerun: '+ntostr(rerun)+'  Camcol:  '+ntostr(camcol)
  print,'   Fields ['+ntostr(ff)+', '+ntostr(lf)+']'

END 

PRO mysql_tsobj_stuff, run, camcol, rerun=rerun, $
                       start=start, nfields=nfields,$
                       status=status, outdir=outdir

  ;; status is 1 unless we reach end
  status =1

  IF n_params() LT 2 THEN BEGIN ;Help message
      print,'-syntax mysql_tsobj_stuff, run, camcol, rerun=rerun, $'
      print,'               start=start, nfields=nfields,$'
      print,'               status=status, outdir=outdir'
      return 
  ENDIF 

  time=systime(1)

  ;; HTM depth
  htm_depth = 17

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; temporary directory for mysql input files
  ;; This wasy we don't always try to write to the same disk, slows things
  ;; down
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; tmpdir = '/net/cheops2/home/esheldon/tmp/'
  tmpdir = '/net/'+sdssidl_config('hostname')+'/data0/esheldon/tmp/'



  ;; what tags to stuff
  taglist = mysql_tsobj_stuff_taglist()

  ;; The files and fields
  mysql_tsobj_stuff_files, run, camcol, infiles, fields, $
    rerun=rerun, $
    start=start, nfields=nfields

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Loop over fields and find psf moments and correct shapes
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nf = n_elements(fields)

  FOR fi = 0L, nf-1 DO BEGIN 

      field = fields[fi]
      fstr=ntostr(field)

      dirsep, infiles[fi], tdir, tname
      tname = repstr(tname, 'tsObj','tsObj-stuffdb')
      dbInputFile = tmpdir + repstr(tname, '.fit','.st')

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Read tsObj file
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      pstruct = sdss_read('tsobj',run, camcol, rerun=rerun,$
                                field=field, $
                                taglist=taglist, $
                                verbose=0, status=status)

      IF status EQ 0 THEN BEGIN 

          ;; We want photoid at the beginning, so don't use ex_struct
          IF n_elements(outst) EQ 0 THEN BEGIN 
              
              outst = create_struct('photoid',0ULL, $
                                    pstruct[0], $
                                    'htm_index', 0ULL)
          ENDIF 
          
          nout = n_elements(pstruct)
          
          IF nout NE 0 THEN BEGIN 
              
              outStruct = replicate(outst, nout)
              struct_assign, pstruct, outStruct, /nozero
              photoid = sdss_photoid(outStruct)
              outStruct.photoid = photoid
              
              outStruct.htm_index = $
                htm_index(outStruct.ra, outStruct.dec, htm_depth)
              
              ;; write file
              print,'Writing to db input file: ',dbInputFile
              write_idlstruct, outStruct, dbInputFile, /ascii, delimiter='tab'
              

              ;; Skip the header. This is necessary because the blank line
              ;; after the header gets read as a record.
              print,'Stuffing database'
              ignorelines = '124'
              query = $
                "LOAD DATA LOCAL INFILE '"+dbInputFile+"' "+$
                "INTO TABLE tsObj IGNORE "+ignorelines+" LINES"

              out = mysql_query(query, user='esheldon', pwd='cheopsDB')
              
              file_delete, dbInputFile
              
          ENDIF 
      ENDIF 

  ENDFOR 

  ptime,systime(1)-time

  status=0

END 
