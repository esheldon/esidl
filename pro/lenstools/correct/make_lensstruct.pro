PRO make_lensstruct, lensstr
  
  IF n_params() LT 1 THEN BEGIN
      print,'-Syntax: make_lensstruct, lensstr'
      return
  ENDIF 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make structure with new tags, including default values
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  val  = -9999.
  val2 =  9999.
  val3 = -9999
  val4 = -9999d
  val5 =  0b
  val6 = -1b
  arrval  = replicate(-9999., 5)
  arrval2 = replicate( 9999., 5)
  arrval3 = replicate( 2b^6,  5)
  lensstr = create_struct(                             $
                           'VALUE_FLAGS',        val5, $
                           'CORRSELECT_FLAGS',   val5, $
                           'OBJC_PROB_FLAGS',    val5, $
                           'CMODEL_COUNTS',    arrval, $ ; corrected model cnt 
                           'CMODEL_COUNTSERR', arrval, $
                           'M_E1_CORR',        arrval, $ ; lensing stuff
                           'M_E2_CORR',        arrval, $
                           'M_R',              arrval, $
                           'M_E1_CORR_H',      arrval, $ ; Hirata's stuff
                           'M_E2_CORR_H',      arrval, $
                           'M_R_H',            arrval, $
                           'COMPEA4FLAGS',    arrval3, $
                           'SEEING',          arrval2, $
                           'ROTATION',        arrval2, $ ; shape rotation
                           'CLAMBDA',            val4, $ ; corrected coords.
                           'CETA',               val4, $
                           'PHOTOZ_Z',            val, $ ; photoz stuff
                           'PHOTOZ_ZERR',         val, $
                           'PHOTOZ_TYPE',         val, $
                           'PHOTOZ_TYPEERR',      val, $
                           'PHOTOZ_COVAR_TT',     val, $
                           'PHOTOZ_COVAR_TZ',     val, $
                           'PHOTOZ_COVAR_ZZ',     val, $
                           'PHOTOZ_CHISQ',       val2, $
                           'PHOTOZ_QUALITY',     val3, $
                           'PHOTOZ_DIST_MOD',     val, $
                           'PHOTOZ_KCORR',     arrval, $
                           'PHOTOZ_ABSCOUNTS', arrval, $
                           'PHOTO_QSO_Z',         val, $
                           'PHOTO_QSO_ZMIN',      val, $
                           'PHOTO_QSO_ZMAX',      val, $
                           'PHOTO_QSO_PROB',      val, $
                           'J_2MASS',             val, $
                           'J_2MASSERR',          val, $
                           'H_2MASS',             val, $
                           'H_2MASSERR',          val, $
                           'K_2MASS',             val, $
                           'K_2MASSERR',          val  $
                         )

;yes/no flags: corrselect_flags (5 bits)
;              objc_prob_flags  (6 bits)
;              compea4flags     (7 bits)
;
;              photo_qso        (1 bits)
;              2mass            (2 bits)
;              valid photoz     (1 bits)
;              first match      (1 bits)
;              rosat match      (1 bits)


;2b^0: VALID BAYESIAN S/G
;2b^1: VALID PHOTOZ
;2b^2: PHOTOMETRIC QSO
;2b^3: 2MASS EXTENDED
;2b^4: 2MASS POINTSOURCE
;2b^5: FIRST MATCH
;2b^6: ROSAT MATCH


return
END 
