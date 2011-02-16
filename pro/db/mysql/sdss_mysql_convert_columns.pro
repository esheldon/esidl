; Check for any of the 5-bandpass columns.  If they are entered,
; convert them to col_u, col_g, ... columns


PRO sdss_mysql_convert_columns, sdss_cols, colflags, nvalues, mysql_cols

  COMMON mysql_convert_block, fivetags, ISFIVE, ISDEC

  IF n_elements(fivetags) EQ 0 THEN BEGIN 

      ISFIVE = 2^0
      ISDEC = 2^1

      ;; This list is created by mysql_tsobj_stuff_fivetags
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
  ENDIF 

  ;; Find the tags that match the "fivetags"
  sdss_cols = strupcase(sdss_cols)
  
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

      colors = ['u','g','r','i','z']
      ii = 0L
      FOR i=0L, ncols-1 DO BEGIN 

          IF colflags[i] THEN BEGIN 
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

END 

