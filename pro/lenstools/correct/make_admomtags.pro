pro make_admomtags,taglist, default=default

  if n_params() eq 0 then begin
      print,'-syntax make_admomtags, taglist, default=default'
      return
  endif
  
  taglist=[                    $
            'RUN',             $
            'RERUN',           $
            'CAMCOL',          $
            'FIELD',           $
            'ID',              $
            'PARENT',          $
            'NCHILD',          $
            'OBJC_TYPE',       $
            'TYPE',            $
            'FLAGS',           $
            'FLAGS2',          $
            'OBJC_FLAGS',      $
            'OBJC_FLAGS2',     $
            'OBJC_ROWC',       $
            'OBJC_COLC',       $
            'ROWC',            $
            'COLC',            $
            'COUNTS_MODEL',    $
            'COUNTS_MODELERR', $
            'COUNTS_EXP',      $
            'COUNTS_EXPERR',   $
            'COUNTS_DEV',      $
            'COUNTS_DEVERR',   $
            'FRACPSF',         $
            'PETROCOUNTS',     $
            'PETROCOUNTSERR',  $
            'PETRORAD',        $
            'PETRORADERR',     $
            'PETROR50',        $
            'PETROR50ERR',     $
            'PETROR90',        $
            'PETROR90ERR',     $
            'PSFCOUNTS',       $
            'PSFCOUNTSERR',    $
            'STATUS',          $
            'RA',              $
            'DEC',             $
            'PRIMTARGET',      $
            'SECTARGET',       $
            'REDDENING',       $
            'M_E1',            $
            'M_E2',            $
            'M_E1E1ERR',       $
            'M_E1E2ERR',       $
            'M_E2E2ERR',       $
            'M_RR_CC',         $
            'M_RR_CCERR',      $
            'M_CR4',           $
            'M_E1_PSF',        $
            'M_E2_PSF',        $
            'M_RR_CC_PSF',     $
            'M_CR4_PSF',       $
            'OBJC_PROB_PSF',   $
            'PROB_PSF',        $
            'SKY',             $
            'SKYERR',          $
            'FIRSTMATCH',      $
            'ROSATMATCH'       $
          ]

  IF NOT keyword_set(default) THEN BEGIN
      addtags = ['VALUE_FLAGS',      $
                 'CORRSELECT_FLAGS', $
                 'OBJC_PROB_FLAGS',  $
                 'M_E1_CORR',        $ ; lensing stuff
                 'M_E2_CORR',        $
                 'M_R',              $
                 'M_E1_CORR_H',      $
                 'M_E2_CORR_H',      $
                 'M_R_H',            $
                 'COMPEA4FLAGS',     $
                 'SEEING',           $
                 'ROTATION',         $ ; shape rotation
                 'CLAMBDA',          $ ; "corrected" survey coordinates
                 'CETA',             $
                 'PHOTOZ_Z',         $ ; photoz stuff
                 'PHOTOZ_ZERR',      $
                 'PHOTOZ_TYPE',      $
                 'PHOTOZ_TYPEERR',   $
                 'PHOTOZ_COVAR_TT',  $
                 'PHOTOZ_COVAR_TZ',  $
                 'PHOTOZ_COVAR_ZZ',  $
                 'PHOTOZ_CHISQ',     $
                 'PHOTOZ_QUALITY',   $
                 'PHOTOZ_DIST_MOD',  $
                 'PHOTOZ_KCORR',     $
                 'PHOTOZ_ABSCOUNTS', $
                 'PHOTO_QSO_Z',      $
                 'PHOTO_QSO_ZMIN',   $
                 'PHOTO_QSO_ZMAX',   $
                 'PHOTO_QSO_PROB',   $
                 'J_2MASS',          $
                 'J_2MASSERR',       $
                 'H_2MASS',          $
                 'H_2MASSERR',       $
                 'K_2MASS',          $
                 'K_2MASSERR'        $
                ]

      taglist = [taglist, addtags]

  ENDIF

  return
end
