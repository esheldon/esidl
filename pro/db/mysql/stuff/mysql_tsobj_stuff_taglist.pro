PRO mysql_tsobj_stuff_taglist_readstruct, tags, struct

  ;; We will stuff everything into the database other than
  ;; The match tags (which will be stored separately) and the
  ;; profmean, proferr which are unwieldy
  
  tmp = sdss_read('tsobj',756,1,verbose=0)
  tmp = tmp[0]

  tags = tag_names(tmp)

  w=where(strmatch(tags,'USNO*') EQ 0 AND $
          strmatch(tags,'FIRST*') EQ 0 AND $
          strmatch(tags,'ROSAT*') EQ 0 AND $
          tags NE 'PROFMEAN' AND tags NE 'PROFERR' AND $
          tags NE 'MATCHID')

  IF arg_present(struct) THEN BEGIN 
      struct = sdss_read('tsobj',756,3,verbose=0, taglist=tags[w])
      struct = struct[0]
      zero_struct, struct
  ENDIF 

  tags = tags[w]

END 

FUNCTION mysql_tsobj_stuff_taglist, struct=struct

  IF arg_present(struct) THEN BEGIN 
      mysql_tsobj_stuff_taglist_readstruct, tags, struct
      return, tags
  ENDIF 

  tags = ['RUN', 'CAMCOL', 'RERUN', 'FIELD', 'PARENT', 'ID', $
          'NCHILD', $
          'OBJC_TYPE', $
          'OBJC_PROB_PSF', $
          'CATID', $
          'OBJC_FLAGS', 'OBJC_FLAGS2', $
          'OBJC_ROWC', 'OBJC_ROWCERR', $
          'OBJC_COLC', 'OBJC_COLCERR', $
          'ROWV', 'ROWVERR', 'COLV', 'COLVERR', $
          'ROWC', 'ROWCERR', 'COLC', 'COLCERR', $
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
          'PROB_PSF', $
          'NPROF', $
          'STATUS', $
          'RA', 'DEC', 'LAMBDA', 'ETA', 'L', 'B', $
          'OFFSETRA', 'OFFSETDEC', $
          'PRIMTARGET', 'SECTARGET', $
          'REDDENING', $
          'PROPERMOTIONMATCH', 'PROPERMOTIONDELTA', $
          'PROPERMOTION', 'PROPERMOTIONANGLE', 'PRIORITY']

  return, tags

END 
