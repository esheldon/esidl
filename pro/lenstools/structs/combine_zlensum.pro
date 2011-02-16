PRO combine_zlensum, lensum, binsize, rminkpc, rmaxkpc, hval, shstruct, depth=depth, logbin=logbin, comoving=comoving, compcut=compcut

  IF n_params() LT 5 THEN BEGIN 
      print,'-Syntax: combine_zlensum, lensum, binsize, rminkpc, rmaxkpc, hval, shstruct, depth=, logbin=, comoving=, compcut='
      return
  ENDIF 

  wz = getztag(lensum)


  nlens = n_elements(lensum)

  s=size(lensum[0].rsum)
  ;; assume 2 or fewer dimensions
  IF s[0] LT 2 THEN nang = 1 ELSE nang = s[2]
  nrad = s[1]

  ;; Should we use the outer bin?
  IF total_int(lensum.npair[nrad-1,*]) LT total_int(lensum.npair[nrad-2,*]) THEN BEGIN 
      arrval = lensum[0].rsum[0:nrad-2,*]
      nrad = nrad-1
  ENDIF ELSE BEGIN
      arrval = lensum[0].rsum
  ENDELSE 

  arrval[*] = 0.0


  shstruct = zshstruct(arrval)

  IF tag_exist(lensum[0], 'rmin_act') THEN DO_rminact=1 ELSE DO_rminact=0

  ;; lens weight
  lensw = lensum.weight
  lenswsum = total( lensw ,/double)
  shstruct.lenswsum = lenswsum

  shstruct.nlenses = nlens
  shstruct.binsize = binsize
  shstruct.rmin    = rminkpc
  shstruct.rmax    = rmaxkpc
  shstruct.h       = hval

  ;; optional
  IF n_elements(depth)    NE 0 THEN shstruct.depth = depth
  IF n_elements(comoving) NE 0 THEN shstruct.comoving = comoving
  IF n_elements(logbin)   NE 0 THEN shstruct.logbin = logbin
  IF n_elements(compcut)  NE 0 THEN shstruct.compcut=compcut

  ;; scalar sums/means
  IF tag_exist(lensum[0], 'totpairs') THEN BEGIN 
      shstruct.totpairs = total_int( lensum.totpairs )
  ENDIF ELSE BEGIN 
      shstruct.totpairs = total_int( lensum.tot_pairs )
  ENDELSE 
  shstruct.wsum_ssh = total(lensum.wsum_ssh, /double)
  shstruct.sshsum = total(lensum.sshsum, /double)
  shstruct.ssh = shstruct.sshsum/shstruct.wsum_ssh


  IF wz NE -1 THEN BEGIN
      shstruct.zsum   = total(lensum.(wz[0])*lensw, /double)
      shstruct.zmean = shstruct.zsum/lenswsum
  ENDIF 

  IF tag_exist(lensum[0], 'scritinv') THEN BEGIN 
      shstruct.scritinvsum = total(lensum.scritinv*lensw, /double)
      shstruct.meanscritinv = shstruct.scritinvsum/lenswsum
  ENDIF 

  IF tag_exist(lensum[0], 'angsum') THEN doangsum=1 ELSE doangsum=0

  FOR radbin=0L, nrad-1 DO BEGIN 
      FOR angbin = 0L, nang-1 DO BEGIN 

          terr = lensum.sigmaerr[radbin,angbin]
          oterr = lensum.orthosigerr[radbin,angbin]
          w = where(terr GT 0.0 AND finite(terr) AND $
                    oterr GT 0.0 AND finite(oterr) AND $
                    abs(lensum.sigma[radbin,angbin]) LT 1000000, nw)
          shstruct.nlbin[radbin,angbin] = nw

          IF nw NE 0 THEN BEGIN 

              ;; sums
              shstruct.npair[radbin,angbin] = $
                total_int(lensum[w].npair[radbin,angbin])
              shstruct.tnpair[radbin,angbin] = $
                total_int(shstruct.npair[0:radbin,angbin])

              IF shstruct.npair[radbin,angbin] EQ 0 THEN BEGIN 
                  message,'npair = 0 when it should not'
              ENDIF 

              wsum = total(lensum[w].wsum[radbin,angbin], /double)
              wsum2 = total(lensum[w].wsum[radbin,angbin]^2, /double)

              shstruct.wsum[radbin,angbin] = wsum
              shstruct.wsum2[radbin,angbin] = wsum2

              wsum_mean =  wsum/nlens
              wsum_err = sqrt(wsum2/nlens - wsum_mean^2)/sqrt(nlens)

              shstruct.wsum_mean[radbin,angbin] = wsum_mean
              shstruct.wsum_err[radbin,angbin] = wsum_err

              shstruct.sigerrsum[radbin,angbin] = $
                total(lensum[w].sigerrsum[radbin,angbin], /double)
              shstruct.orthosigerrsum[radbin,angbin] = $
                total(lensum[w].orthosigerrsum[radbin,angbin], /double)

              shstruct.rsum[radbin,angbin] = total(lensum[w].rsum[radbin,angbin], /double)
              IF doangsum THEN BEGIN 
                  shstruct.angsum[radbin,angbin] = $
                    total(lensum[w].angsum[radbin,angbin], /double)
                  shstruct.meanang[radbin,angbin] = $
                    shstruct.angsum[radbin,angbin]/shstruct.npair[radbin,angbin]
              ENDIF 

              ;; means
              shstruct.meanr[radbin,angbin]    = $
                shstruct.rsum[radbin,angbin]/shstruct.npair[radbin,angbin]

              shstruct.rmax_act[radbin,angbin] = max(lensum[w].rmax_act[radbin,angbin])
              IF DO_rminact THEN BEGIN
                  wt=where(lensum[w].rmin_act[radbin,angbin] GT 0.)
                  shstruct.rmin_act[radbin,angbin] = $
                    min(lensum[w[wt]].rmin_act[radbin,angbin])

                  R1 = shstruct.rmin_act[radbin,angbin]
                  R2 = shstruct.rmax_act[radbin,angbin]

              ENDIF ELSE BEGIN 
                  R1 = shstruct.rmin + radbin*shstruct.binsize
                  R2 = shstruct.rmin + (radbin+1)*shstruct.binsize
              ENDELSE 

              shstruct.area[radbin,angbin] = !pi*(R2^2 - R1^2)
              shstruct.tarea[radbin,angbin] = !pi*(R2^2 - shstruct.rmin^2)
              IF shstruct.area[radbin,angbin] GT 0. THEN BEGIN 
                  shstruct.density[radbin,angbin] = $
                    shstruct.npair[radbin,angbin]/shstruct.area[radbin,angbin]/shstruct.nlenses
              ENDIF 

              IF nw GT 10 THEN calcerr=1 ELSE calcerr=0
calcerr=0
              wmom, $
                lensum[w].sigma[radbin,angbin], $
                lensum[w].sigmaerr[radbin,angbin], $
                wmean, wsig, werr, calcerr=calcerr

              shstruct.sigma[radbin,angbin] = wmean
              shstruct.sigmaerr[radbin,angbin] = werr

              wmom, $
                lensum[w].orthosig[radbin,angbin], $
                lensum[w].orthosigerr[radbin,angbin], $
                owmean, owsig, owerr, calcerr=calcerr

              shstruct.orthosig[radbin,angbin] = owmean
              shstruct.orthosigerr[radbin,angbin] = owerr

              ;; alternative error estimate

;              IF shstruct.wsum[radbin,angbin] EQ 0.0 THEN message,'wsum zero'
;              shstruct.sigmaerr2[radbin,angbin] = $
;                sqrt(shstruct.sigerrsum[radbin,angbin])/shstruct.wsum[radbin,angbin]
;              IF shstruct.owsum[radbin,angbin] EQ 0.0 THEN message,'owsum zero'
;              shstruct.orthosigerr2[radbin,angbin] = $
;                sqrt(shstruct.orthosigerrsum[radbin,angbin])/shstruct.owsum[radbin,angbin]

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; cumulative means
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              shstruct.tmeanr[radbin,angbin] = $
                total(shstruct.rsum[0:radbin,angbin],/double)/total_int(shstruct.npair[0:radbin,angbin])

              wmom, $
                shstruct.sigma[0:radbin,angbin],$
                shstruct.sigmaerr[0:radbin,angbin], $
                wmean, wsig, werr

              shstruct.tsigma[radbin,angbin] = wmean
              shstruct.tsigmaerr[radbin,angbin] = werr

              wmom, $
                shstruct.orthosig[0:radbin,angbin], $
                shstruct.orthosigerr[0:radbin,angbin], $
                owmean, owsig, owerr

              shstruct.torthosig[radbin,angbin] = owmean
              shstruct.torthosigerr[radbin,angbin] = owerr

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; alternative error estimate
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;              shstruct.tsigmaerr2[radbin,angbin] = $
;                sqrt(total(shstruct.sigerrsum[0:radbin,angbin],/double))/total(shstruct.wsum[0:radbin,angbin],/double)

;              shstruct.torthosigerr2[radbin,angbin] = $
;                sqrt(total(shstruct.orthosigerrsum[0:radbin,angbin],/double))/total(shstruct.owsum[0:radbin,angbin],/double)

          ENDIF 
                                          
      ENDFOR 
  ENDFOR 

  return
END 
