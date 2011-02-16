pro  starmaghist, yhist, xhist, plt=plt


  if N_params() eq 0 then begin
    print,'Syntax: starmaghist,yhist,xhist,plt=plt'
    return
  endif

  file = '~/idl.lib/SIMSDSS/FIT/starmag.fit'
  starmag = mrdfits(file, 1,hdr,/silent)

  yhist = starmag.maghist
  xhist = starmag.mag

  IF keyword_set(plt) THEN plot,xhist,yhist,psym=4,$
    xtitle='Mag_best', ytitle='Number',title='Star Magnitude Distribution'


return
end










