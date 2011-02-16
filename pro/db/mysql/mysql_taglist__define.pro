FUNCTION mysql_taglist::init

  fivetags = ['ROWC', 'ROWCERR', 'COLC', 'COLCERR', $
              'SKY', 'SKYERR', $
              'PSFCOUNTS', 'PSFCOUNTSERR', $
              'FIBERCOUNTS', 'FIBERCOUNTSERR', $
              'PETROCOUNTS', 'PETROCOUNTSERR', $
              'PETRORAD', 'PETRORADERR', $
              'PETROR50', 'PETROR50ERR', 'PETROR90', 'PETROR90ERR', $
              'Q', 'QERR', 'U', 'UERR', $
              'M_E1', 'M_E2', 'M_E1E1ERR', 'M_E1E2ERR', 'M_E2E2ERR', $
              'M_RR_CC', 'M_RR_CCERR', 'M_CR4', $
              'M_E1_PSF', 'M_E2_PSF', $
              'M_RR_CC_PSF', 'M_CR4_PSF', $
              'ISO_ROWC', 'ISO_ROWCERR', 'ISO_ROWCGRAD', $
              'ISO_COLC', 'ISO_COLCERR', 'ISO_COLCGRAD', $
              'ISO_A', 'ISO_AERR', 'ISO_AGRAD', $
              'ISO_B', 'ISO_BERR', 'ISO_BGRAD', $
              'ISO_PHI', 'ISO_PHIERR', 'ISO_PHIGRAD', $
              'R_DEV', 'R_DEVERR', $
              'AB_DEV', 'AB_DEVERR', $
              'PHI_DEV', 'PHI_DEVERR', $
              'COUNTS_DEV', 'COUNTS_DEVERR', $
              'R_EXP', 'R_EXPERR', $
              'AB_EXP', 'AB_EXPERR', $
              'PHI_EXP', 'PHI_EXPERR', $
              'COUNTS_EXP', 'COUNTS_EXPERR', $
              'COUNTS_MODEL', 'COUNTS_MODELERR', $
              'TEXTURE', $
              'STAR_L', 'STAR_LNL', $
              'EXP_L', 'EXP_LNL', $
              'DEV_L', 'DEV_LNL', $
              'FRACPSF', $
              'FLAGS', 'FLAGS2', $
              'TYPE', $
              'PROB_PSF', 'NPROF', $
              'OFFSETRA', 'OFFSETDEC', $
              'REDDENING']

  self.fivetags_ptr = ptr_new(fivetags, /no_copy)

  self.isfive = 2^0
  self.isdec = 2^1

  return, 1

END 

PRO mysql_taglist::set, tags, type
  
  np = n_params()
  n_tags = n_elements(tags)
  n_type = n_elements(type)
  IF np LT 2 OR n_type EQ 0 OR n_tags EQ 0 THEN BEGIN 
      message,'-Syntax: a->set(tags, type)',/inf
      return
  ENDIF 

  IF n_tags NE 0 THEN BEGIN 

      ;; If tags already exists, then free it
      ptr_free, self.tags_ptr

      self.tags_ptr = ptr_new(strupcase(tags), /no_copy)
      self.n_tags = n_tags

      self.type = strupcase(type)
  ENDIF 

  IF self.type EQ 'SDSS' THEN BEGIN 
      self->matchfive
  ENDIF ELSE BEGIN 
      self.nvalues_ptr = ptr_new( replicate(1, self.n_tags), /no_copy)
  ENDELSE 

  return

END 

FUNCTION mysql_taglist::tags

  IF ptr_valid(self.tags_ptr) THEN BEGIN 
      IF n_elements( *(self.tags_ptr) ) NE 0 THEN BEGIN 
          return, *(self.tags_ptr)
      ENDIF ELSE BEGIN 
          message,'tags are uninitialized'
      ENDELSE 
  ENDIF ELSE BEGIN 
      message,'tags are uninitialized'
  ENDELSE 

END 

FUNCTION mysql_taglist::nvalues
  return, *(self.nvalues_ptr)
END 

FUNCTION mysql_taglist::type
  return, self.type
END 

PRO mysql_taglist::matchfive

  tags = *(self.tags_ptr)
  fivetags = *(self.fivetags_ptr)
  ntags = n_elements(tags)

  match, tags, fivetags, mtags, mfive
  
  nvalues = replicate(1, ntags)
  IF mtags[0] NE -1 THEN nvalues[mtags] = 5

  self.nvalues_ptr = ptr_new(nvalues)

END 

PRO mysql_taglist::convert, subscripts

  IF self.type EQ 'SDSS' THEN BEGIN 
      fivetags = *(self.fivetags_ptr)
      ISFIVE = self.isfive
      ISDEC = self.isdec

      ;; Find the tags that match the "fivetags"
      
      sdss_cols = *(self.tags_ptr)

      ncols = n_elements(sdss_cols)
      colflags = intarr(ncols)

      match, sdss_cols, fivetags, mcol, mfive

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; The ones that match are colflags = 1 and are replaced by _bandpass
      ;; tags in mysql call.  Otherwise, we will just return copy of 
      ;; the original
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      IF mfive[0] NE -1 THEN BEGIN 

          colflags[mcol] = ISFIVE
          nfive = n_elements(mfive)
          
          nmysql_cols = ncols + 4*nfive
          mysql_cols = strarr(nmysql_cols)
          nvalues = intarr(ncols)
          
          colors = ['U','G','R','I','Z']
          ii = 0L
          FOR i=0L, ncols-1 DO BEGIN 
              
              IF colflags[i] EQ ISFIVE THEN BEGIN 
                  mysql_cols[ii:ii+4] = sdss_cols[i] + '_'+colors
                  ii = ii+5
                  nvalues[i] = 5
              ENDIF ELSE BEGIN 
                  mysql_cols[ii] = sdss_cols[i]
                  ii = ii+1
                  nvalues[i] = 1
              ENDELSE 
              
          ENDFOR 

      ENDIF ELSE BEGIN 
          mysql_cols = sdss_cols
      ENDELSE 

      w=where(mysql_cols EQ 'DEC', nw)
      IF nw NE 0 THEN BEGIN 
          mysql_cols[w] = 'DECL'
          colflags[w] = colflags[w] + ISDEC
      ENDIF 

      self.n_tags = n_elements(mysql_cols)

      ptr_free, self.tags_ptr, self.nvalues_ptr
      self.tags_ptr = ptr_new(mysql_cols)
      self.nvalues_ptr = ptr_new( replicate(1, self.n_tags) )

      self.type = 'MYSQL'
      
  ENDIF ELSE IF self.type EQ 'MYSQL' THEN BEGIN 

      fivetags = *(self.fivetags_ptr)

      mysql_tags = *(self.tags_ptr)
      sdss_tags = strarr(self.n_tags)
      nvalues = intarr(self.n_tags)

      clrend = ['_U','_G','_R','_I','_Z']
      FOR i=0L, self.n_tags-1 DO BEGIN 

          send = strmid(mysql_tags[i], 1, /reverse_offset)
          match, send, clrend, msend, mclrend

          IF msend[0] NE -1 THEN BEGIN 
              slen = strlen(mysql_tags[i])
              sdss_tags[i] = strmid(mysql_tags[i], 0, slen-2)
              nvalues[i] = 5
          ENDIF ELSE BEGIN
              sdss_tags[i] = mysql_tags[i]
              nvalues[i] = 1
          ENDELSE 
      ENDFOR 

      ;; Get unique
      srt = sort(sdss_tags)
      unq = UNIQ(sdss_tags[srt])

      ;; Unsort
      uq_srt = srt[unq]
      uq_unsrt = uq_srt[ sort(uq_srt) ]
      sdss_tags = sdss_tags[uq_unsrt]
      nvalues = nvalues[uq_unsrt]

      ;; can output subscripts
      IF arg_present(subscripts) THEN subscripts = uq_unsrt

      ;; Convert decl to dec
      w = where(sdss_tags EQ 'DECL',nw)
      IF nw NE 0 THEN sdss_tags[w] = 'DEC'

      ptr_free, self.tags_ptr
      self.n_tags = n_elements(sdss_tags)
      self.tags_ptr = ptr_new(sdss_tags, /no_copy)
      self.nvalues_ptr = ptr_new(nvalues, /no_copy)

      self.type = 'SDSS'

  ENDIF ELSE BEGIN 
      message,'Cannot convert unknown type '+self.type
  ENDELSE 

END 


FUNCTION mysql_taglist::cleanup

;-- free memory allocated to pointer when destroying object

 ptr_free,self.tags_ptr, self.nvalues_ptr, self.fivetags_ptr

 return, 1

END 

PRO mysql_taglist__define

  struct = {$
             mysql_taglist, $
             tags_ptr: ptr_new(), $
             nvalues_ptr: ptr_new(), $
             n_tags:0L, $
             type: '', $
             fivetags_ptr: ptr_new(), $
             isfive: 2^0, $
             isdec: 2^1 $
           }

  return
END 
