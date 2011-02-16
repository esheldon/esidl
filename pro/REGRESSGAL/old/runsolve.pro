PRO runsolve, matrixfile, evectfile, rsumfile, get=get

  IF n_params() LT 3 THEN BEGIN
      print,'-Syntax: runsolve, matrixfile, evectfile, rsumfile'
      return
  ENDIF 

  nbin = 16L

  IF NOT keyword_set(get) THEN get=0

  ;; Do filenames

  tmp=str_sep(matrixfile, 'matrix')

  atafile = tmp[0]+'ata'+tmp[1]
  btbfile = tmp[0]+'btb'+tmp[1]
  outfile = tmp[0]+'regshears'+tmp[1]

  IF (exist(atafile) OR exist(btbfile) OR exist(outfile) ) AND $
    (NOT get) THEN BEGIN 
      atafile = newname(atafile)
      btbfile = newname(btbfile)
      outfile = newname(outfile)
  ENDIF 

  print
  print,'atafile: ',atafile
  print,'btbfile: ',btbfile
  print,'outfile: ',outfile
  print

  IF NOT get THEN BEGIN 
      ;; create ata, btb
      print,'Running wm_trans_m'
      wm_trans_m,matrixfile,evectfile,ata,btb,nbin
      mwrfits,ata,atafile
      mwrfits,btb,btbfile
  ENDIF ELSE BEGIN 
      IF exist(atafile) THEN BEGIN 
          ata = mrdfits(atafile)
      ENDIF ELSE BEGIN 
          print,'No such atafile: ',atafile
          return 
      ENDELSE 
      IF exist(btbfile) THEN BEGIN 
          btb = mrdfits(btbfile)
      ENDIF ELSE BEGIN 
          print,'No such atafile: ',btbfile
          return 
      ENDELSE 
  ENDELSE 

  ;; get "uncertainties"  Only a crude approximation and only for 
  ;; normally distributed errors.

  uncert1 = fltarr(nbin)
  uncert2 = fltarr(nbin)

  a_inv=invert(ata)
  b_inv=invert(btb)

  FOR i=0, nbin-1 DO BEGIN 
      uncert1[i] = sqrt(a_inv[i,i])/2.
      uncert2[i] = sqrt(b_inv[i,i])/2.
  ENDFOR 

  print,'Running wsolve_reg'
  wsolve_reg,matrixfile,evectfile,nbin,ata,btb,$
    etan1,erad1,etan2,erad2

  shear1 = etan1/2.
  ortho1 = erad1/2.
  shear2 = etan2/2.
  ortho2 = erad2/2.

  readcol, rsumfile, meanr

  forprint, meanr, shear1, uncert1, shear2, uncert2, ortho1, ortho2

  openw, unit, outfile, /get_lun
  !textunit=unit
  fmt='(7(F0,:,1X))'
  forprint, meanr, shear1, uncert1, shear2, uncert2, ortho1, ortho2, $
    TEXT=5,F=fmt,/silent

  free_lun, unit


return
END 
