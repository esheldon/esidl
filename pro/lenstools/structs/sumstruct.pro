FUNCTION sumstruct, arrval

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: sumstruct = sumstruct(arrval)'
      return, -1
  ENDIF 

  sumstruct = create_struct('nlenses', 0., $
                            'totpairs', 0., $
                            'binsize', 0., $
                            'rmin', 0., $
                            'rmax', 0., $
                            'rmax_act', arrval, $
                            'sshsum', 0., $
                            'lenswsum', 0., $
                            'rsum', arrval, $
                            'etansum', arrval, $
                            'etanerrsum', arrval, $
                            'eradsum', arrval, $
                            'eraderrsum', arrval)

  sumstruct = create_struct(sumstruct, $
                            'qrsum',arrval,$
                            'qrerrsum',arrval,$
                            'qisum',arrval,$
                            'qierrsum',arrval,$
                            'qradrsum',arrval,$
                            'qradrerrsum',arrval,$
                            'qradisum',arrval,$
                            'qradierrsum',arrval,$
                            $
                            'etan1sum',arrval,$
                            'etan1errsum',arrval,$
                            'etan2sum',arrval,$
                            'etan2errsum',arrval,$
                            'erad1sum',arrval,$
                            'erad1errsum',arrval,$
                            'erad2sum',arrval,$
                            'erad2errsum',arrval,$
                            'w1sum',arrval,$
                            'w2sum',arrval,$
                            $
                            'wsum', arrval, $
                            'npair', arrval)

  return, sumstruct

END 
