PRO momname, run, camcol, cindex, iorder, fname

  IF n_params() LT 4 THEN BEGIN
      print,'-Syntax: momname, run, camcol, cindex, iorder, fname'
      return
  ENDIF 
  colors = ['u','g','r','i','z']
  rstr = ntostr(run)
  cstr = ntostr(camcol)
  iorderstr = ntostr(iorder)

  fname = 'mom'+rstr+'-col'+cstr+'-'+colors[cindex]+$
         '-order'+iorderstr+'.txt'

return
END 
