FUNCTION make_corrected_struct

  ;; added htm_index because we can't index tsObj right now

  val  = -9999.
  val2 =  9999.
  val3 = -9999
  val4 = -9999d
  val5 =  0b
  val6 = -1b
  val7 = 0ULL
  val8 = -9999L
  arrval  = replicate(-9999., 5)
  arrval2 = replicate( 9999., 5)
  arrval3 = replicate( 2b^6,  5)
  lensstr = create_struct(                             $
                           'photoid',            val7, $
                           'run',                val8, $
                           'rerun',              val3, $
                           'camcol',             val3, $
                           'field',              val3, $
                           'id',                 val8, $
                           'VALUE_FLAGS',        val5, $
                           'CORRSELECT_FLAGS',   val5, $
                           'OBJC_PROB_GAL',      val2, $
                           'OBJC_PROB_FLAGS',    val5, $
                           'CMODEL_COUNTS',    arrval, $ ; combined model cnt 
                           'CMODEL_COUNTSERR', arrval, $
                           'CMODEL_COUNTS_EXT',arrval, $ ; ext. corr. cmodel
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
                           'htm_index',          val7  $ 
                         )

  return,lensstr

END 
