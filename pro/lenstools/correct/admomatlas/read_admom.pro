PRO read_admom, file,struct

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: read_admom, file,struct'
      return
  ENDIF 

  ;; read the output of the C-version of admomatlas

  t=systime(1)
  
  rdfloat,file,$
    run,rerun,camcol,field,id,$
    ixx1,ixx2,ixx3,$
    ixy1,ixy2,ixy3,$
    iyy1,iyy2,iyy3,$
    err1,err2,err3,$
    rho1,rho2,rho3,$
    why1,why2,why3

  ptime,systime(1)-t

  nline=n_elements(run)
  
  s=create_struct('run',0L,$
                  'rerun',0,$
                  'camcol',0,$
                  'field',0,$
                  'id',0L,$
                  'ixx',fltarr(5),$
                  'ixy',fltarr(5),$
                  'iyy',fltarr(5),$
                  'momerr',fltarr(5),$
                  'rho4',fltarr(5),$
                  'whyflag',lonarr(5) )
  
  struct = replicate(s,nline)

  struct.run = temporary(run)
  struct.rerun = temporary(rerun)
  struct.camcol = temporary(camcol)
  struct.field = temporary(field)
  struct.id = temporary(id)

  struct.ixx[1] = temporary(ixx1)
  struct.ixx[2] = temporary(ixx2)
  struct.ixx[3] = temporary(ixx3)
  struct.ixy[1] = temporary(ixy1)
  struct.ixy[2] = temporary(ixy2)
  struct.ixy[3] = temporary(ixy3)
  struct.iyy[1] = temporary(iyy1)
  struct.iyy[2] = temporary(iyy2)
  struct.iyy[3] = temporary(iyy3)

  struct.momerr[1] = temporary(err1)
  struct.momerr[2] = temporary(err2)
  struct.momerr[3] = temporary(err3)

  struct.rho4[1] = temporary(rho1)
  struct.rho4[2] = temporary(rho2)
  struct.rho4[3] = temporary(rho3)

  struct.whyflag[1] = temporary(why1)
  struct.whyflag[2] = temporary(why2)
  struct.whyflag[3] = temporary(why3)
  
  return



END 
