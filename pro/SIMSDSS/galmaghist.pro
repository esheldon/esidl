pro  galmaghist, yhist, xhist, plt=plt

  IF n_params() EQ 0 THEN BEGIN
    print,'-Syntax: galmaghist, yhist, xhist, plt=plt'
    return
  ENDIF
  
  file = '~/idl.lib/SIMSDSS/FIT/galmag.fit'
  galmag = mrdfits(file,1,hdr,/silent)

  xhist = galmag.mag
  yhist = galmag.maghist

  IF keyword_set(plt) THEN plot,xhist,yhist,psym=4,$
    xtitle='Mag_best', ytitle='Number',title='Galaxy Magnitude Distribution'


return
end










