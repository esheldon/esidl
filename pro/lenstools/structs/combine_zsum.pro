PRO combine_zsum, sumstruct, outstruct

  IF n_params() LT 1 THEN BEGIN
      print,'-Syntax: combine_zsum, sumstruct, outstruct'
      return
  ENDIF 

  nbin = n_elements(sumstruct[0].rsum)

  s=size(sumstruct.rsum)
  ;; assume 2 or fewer dimensions
  IF s[0] LT 2 THEN nang = 1 ELSE nang = s[2]
  nrad = s[1]

  arrval = sumstruct.rsum
  arrval[*] = 0.0
  outstruct = zshstruct(arrval)

  struct_assign, sumstruct, outstruct

  outstruct.totpairs = total(sumstruct.npair)
  outstruct.Ssh = sumstruct.Sshsum/sumstruct.wsum_ssh

  IF tag_exist(sumstruct, 'zsum') THEN $
    outstruct.zmean = sumstruct.zsum/sumstruct.lenswsum

  IF tag_exist(sumstruct, 'scritinvsum') THEN $
    outstruct.meanscritinv = sumstruct.scritinvsum/sumstruct.lenswsum

  IF tag_exist(sumstruct, 'rmin_act') THEN BEGIN
      IF total(sumstruct.rmin_act) NE 0. THEN DO_rminact=1 ELSE DO_rminact=0
  ENDIF ELSE DO_rminact=0

  npairold = 0
  FOR radbin=0L, nrad-1 DO BEGIN 
      FOR angbin = 0L, nang-1 DO BEGIN 
          ;; calculate area and density of background galaxies
          IF DO_rminact THEN BEGIN 
              R1 = sumstruct.rmin_act[radbin,angbin]
              R2 = sumstruct.rmax_act[radbin,angbin]
          ENDIF ELSE BEGIN 
              R1 = outstruct.rmin + radbin*outstruct.binsize
              R2 = outstruct.rmin + (radbin+1)*outstruct.binsize
          ENDELSE 
          outstruct.npair[radbin,angbin] = sumstruct.npair[radbin,angbin]
          outstruct.tnpair[radbin,angbin] = total(sumstruct.npair[0:radbin,angbin])
          outstruct.area[radbin,angbin] = !pi*(R2^2 - R1^2)
          outstruct.tarea[radbin,angbin] = !pi*(R2^2 - outstruct.rmin^2)
          
          IF outstruct.area[radbin,angbin] NE 0. THEN BEGIN 
              outstruct.density[radbin,angbin] = $
                outstruct.npair[radbin,angbin]/outstruct.area[radbin,angbin]/outstruct.nlenses
          ENDIF 
          outstruct.rmax_act[radbin,angbin] = sumstruct.rmax_act[radbin,angbin]
          IF DO_rminact THEN outstruct.rmin_act[radbin,angbin] = $
            sumstruct.rmin_act[radbin,angbin]
          outstruct.meanr[radbin,angbin]    = $
            sumstruct.rsum[radbin,angbin]/sumstruct.npair[radbin,angbin]
          IF tag_exist(sumstruct, 'angsum') THEN BEGIN 
              outstruct.meanang[radbin,angbin] = $
                sumstruct.angsum[radbin,angbin]/sumstruct.npair[radbin,angbin]
          ENDIF 
          
          outstruct.shear[radbin,angbin]    = $
            sumstruct.etansum[radbin,angbin]/sumstruct.wsum[radbin,angbin]/2.
          outstruct.ortho[radbin,angbin]    = $
            sumstruct.eradsum[radbin,angbin]/sumstruct.wsum[radbin,angbin]/2.
          outstruct.shearerr[radbin,angbin] = $
            sqrt( sumstruct.etanerrsum[radbin,angbin]/$
                  sumstruct.wsum[radbin,angbin]^2)/2.
          outstruct.orthoerr[radbin,angbin] = $
            sqrt( sumstruct.eraderrsum[radbin,angbin]/$
                  sumstruct.wsum[radbin,angbin]^2)/2.
          
          outstruct.sigma[radbin,angbin]    = $
            sumstruct.tansigsum[radbin,angbin]/sumstruct.wsum[radbin,angbin]/2.
          outstruct.sigmaerr[radbin,angbin] = $
            sqrt( sumstruct.tansigerrsum[radbin,angbin]/$
                  sumstruct.wsum[radbin,angbin]^2 )/2.
          outstruct.orthosig[radbin,angbin] = $
            sumstruct.radsigsum[radbin,angbin]/$
            sumstruct.wsum[radbin,angbin]/2.
          outstruct.orthosigerr[radbin,angbin] = $
            sqrt( sumstruct.radsigerrsum[radbin,angbin]/$
                  sumstruct.wsum[radbin,angbin]^2 )/2.
          

          ;; Now totals within each rmax_act[radbin,angbin]
          outstruct.tmeanr[radbin,angbin] = total(sumstruct.rsum[0:radbin,angbin])/total(sumstruct.npair[0:radbin,angbin])
          
          outstruct.tshear[radbin,angbin] = total(sumstruct.etansum[0:radbin,angbin])/total(sumstruct.wsum[0:radbin,angbin])/2.
          outstruct.tshearerr[radbin,angbin] = sqrt( total(sumstruct.etanerrsum[0:radbin,angbin])/total(sumstruct.wsum[0:radbin,angbin])^2 )/2.
          outstruct.tortho[radbin,angbin] = total(sumstruct.eradsum[0:radbin,angbin])/total(sumstruct.wsum[0:radbin,angbin])/2.
          outstruct.torthoerr[radbin,angbin] = sqrt( total(sumstruct.eraderrsum[0:radbin,angbin])/total(sumstruct.wsum[0:radbin,angbin])^2 )/2.
          
          outstruct.tsigma[radbin,angbin] = total(sumstruct.tansigsum[0:radbin,angbin])/total(sumstruct.wsum[0:radbin,angbin])/2.
          outstruct.tsigmaerr[radbin,angbin] = sqrt( total(sumstruct.tansigerrsum[0:radbin,angbin])/total(sumstruct.wsum[0:radbin,angbin])^2 )/2.
          outstruct.torthosig[radbin,angbin] = total(sumstruct.radsigsum[0:radbin,angbin])/total(sumstruct.wsum[0:radbin,angbin])/2.
          outstruct.torthosigerr[radbin,angbin] = sqrt( total(sumstruct.radsigerrsum[0:radbin,angbin])/total(sumstruct.wsum[0:radbin,angbin])^2 )/2.
      ENDFOR 

  ENDFOR 

  return
END 
