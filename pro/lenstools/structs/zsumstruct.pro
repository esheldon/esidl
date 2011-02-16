FUNCTION zsumstruct, arrval

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: sumstruct = zsumstruct(arrval)'
      return, -1
  ENDIF 

  tmparrval = arrval
  tmparrval[*] = 1.e7
  sumstruct = create_struct('nlenses', 0., $
                            'totpairs', 0., $
                            'binsize', 0., $
                            'compcut', 0.0, $
                            'rmin', 0., $
                            'rmax', 0., $
                            'rmax_act', arrval, $
                            'rmin_act', tmparrval,$
                            'sshsum', 0., $
                            'lenswsum', 0., $
                            'zsum', 0., $
                            'scritinvsum', 0., $
                            'h', 0., $
                            'angsum', arrval, $
                            'rsum', arrval, $
                            'etansum', arrval, $
                            'etanerrsum', arrval, $
                            'eradsum', arrval, $
                            'eraderrsum', arrval, $
                            'tansigsum', arrval, $
                            'tansigerrsum', arrval, $
                            'radsigsum', arrval, $
                            'radsigerrsum', arrval, $
                            'wsum', arrval, $
                            'wsum_ssh', 0., $
                            'npair', arrval)

  return, sumstruct

END 
