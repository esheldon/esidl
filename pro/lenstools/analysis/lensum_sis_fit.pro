PRO lensum_sis_fit, lensum, density, binsize, rminkpc, rmaxkpc, hval, outstruct

  IF n_params() LT 1 THEN BEGIN
      print,'-Syntax: lensum_sis_fit, lensum, density, binsize, rminkpc, rmaxkpc, hval, outstruct'
      return
  ENDIF 

  time=systime(1)

  IF NOT tag_exist(lensum,'z1d',index=wz) THEN BEGIN
      IF NOT tag_exist(lensum,'z',index=wz) THEN BEGIN 
          IF NOT tag_exist(lensum,'photoz',index=wz) THEN BEGIN
              print,'Lens structure must have "Z" or "PHOTOZ" or "Z1D" flag'
              return
          ENDIF 
      ENDIF 
  ENDIF 

  ;; fit mass to each lens
  nlens = n_elements(lensum)
  nrad = n_elements(lensum[0].rmax_act)

  IF NOT tag_exist(lensum, 'sismass') THEN BEGIN
      tstruct = create_struct('sismass',fltarr(nrad), $
                              'sismasserr', fltarr(nrad), $
                              'sissigma2', fltarr(nrad), $
                              'sissigmaerr2', fltarr(nrad), $
                              'corr', replicate(1., nrad))
      addstruct = replicate(tstruct, nlens)
      combine_structs, lensum, addstruct, outstruct
  ENDIF ELSE BEGIN 
      outstruct = lensum
  ENDELSE 

  nbin = n_elements(lensum[0].rsum)
  corr = fltarr(nbin)

  FOR i=0L, nlens-1 DO BEGIN 



      sumstruct = 0
      shstruct = 0
      combine_zlensum, outstruct[i], binsize, rminkpc, rmaxkpc, hval, $
        sumstruct, shstruct
      
      tsigma_sis_fit, shstruct, sismass, sismasserr, sissigma2, sissigmaerr2
      
      dlens = angdist_lambda(outstruct[i].(wz))*1000. ;kpc
      ;; area in arcminutes^2
      area = shstruct.tarea/dlens^2 ;radians^2
      area = area*(180./!pi*60.)^2 ;arcminutes^2
      randnum = area*density

      w=where(randnum GT 0.)
      corr[*] = 1.
      corr[w] = shstruct.tnpair[w]/randnum[w]
      ssh = shstruct.ssh
      outstruct[i].sismass = sismass*corr/ssh
      outstruct[i].sismasserr = sismasserr*corr/ssh
      outstruct[i].sissigma2 = sissigma2*corr/ssh
      outstruct[i].sissigmaerr2 = sissigmaerr2*corr/ssh
      outstruct[i].corr = corr

      IF (i MOD 20) EQ 0 THEN print,'.',format='($,a)'

  ENDFOR 
  print
  ptime,systime(1)-time

  return
END 
