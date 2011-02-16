PRO frac_overdense, lenspair, nlens, randpair, nrand, frac, fracerr, diff, differr

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: frac_overdense, lenspair, nlens, randpair, nrand, frac, fracerr, diff, differr'
      return
  ENDIF 

  lpair = lenspair/float(nlens)
  rpair = randpair/float(nrand)

  lerr = sqrt(lenspair)/nlens
  rerr = sqrt(randpair)/nrand

  frac = lpair/rpair - 1.

  fracerr = (lpair/rpair)*sqrt( (lerr/lpair)^2 + (rerr/rpair)^2 )

  diff = lpair-rpair
  differr = sqrt( (lerr)^2 + (rerr)^2 )

  return
END 
