FUNCTION combine_shstructs, instruct

  outstruct = instruct[0]
  if not tag_exist(outstruct, 'meanscritinv') then begin
      outstruct = create_struct(outstruct, 'meanscritinv', 0d)
  endif

  IF tag_exist(outstruct,'clust_corr') THEN docorr=1 ELSE docorr=0

  arrval = outstruct.rsum
  arrval[*] = 0.0
  s=size(arrval)
  ;; assume 2 or fewer dimensions
  IF s[0] LT 2 THEN nang = 1 ELSE nang = s[2]
  nrad = s[1]

  ;; running scalar sums/means
  lenswsum = total(instruct.lenswsum, /double)
  outstruct.lenswsum = lenswsum

  nlens = total_int(instruct.nlenses)
  outstruct.nlenses = nlens

  outstruct.totpairs = total_int( instruct.totpairs )

  outstruct.wsum_ssh = total(instruct.wsum_ssh, /double)
  outstruct.sshsum = total(instruct.sshsum, /double)
  outstruct.ssh = outstruct.sshsum/outstruct.wsum_ssh

  outstruct.zsum   = total(instruct.zsum, /double)
  outstruct.zmean = outstruct.zsum/lenswsum

  outstruct.wscritinvsum = total(instruct.wscritinvsum, 2, /double)
  outstruct.meanscritinv = total(outstruct.wscritinvsum, /double)/lenswsum

  ;; now loop over the bins
  IF tag_exist(instruct[0], 'angsum') THEN doangsum=1 ELSE doangsum=0

  FOR radbin=0L, nrad-1 DO BEGIN 
      FOR angbin = 0L, nang-1 DO BEGIN 

          terr = instruct.sigmaerr[radbin,angbin]
          oterr = instruct.orthosigerr[radbin,angbin]
          w = where(terr NE 0.0 AND finite(terr) AND $
                    oterr NE 0.0 AND finite(oterr), nw)

          IF nw NE 0 THEN BEGIN 

              outstruct.nlbin[radbin,angbin] = total_int( instruct[w].nlbin[radbin,angbin] )

              ;; sums
              outstruct.npair[radbin,angbin] = total_int(instruct[w].npair[radbin,angbin])
              outstruct.tnpair[radbin,angbin] = total_int(outstruct.npair[0:radbin,angbin])

              wsum = total(instruct[w].wsum[radbin,angbin], /double)
              wsum2 = total(instruct[w].wsum2[radbin, angbin], /double)

              outstruct.wsum[radbin,angbin] = wsum
              outstruct.wsum2[radbin,angbin] = wsum2

              wsum_mean =  wsum/nlens
              wsum_err = sqrt(wsum2/nlens - wsum_mean^2)/sqrt(nlens)

              outstruct.wsum_mean[radbin,angbin] = wsum_mean
              outstruct.wsum_err[radbin,angbin] = wsum_err

              outstruct.sigerrsum[radbin,angbin] = $
                total(instruct[w].sigerrsum[radbin,angbin], /double)
              outstruct.orthosigerrsum[radbin,angbin] = $
                total(instruct[w].orthosigerrsum[radbin,angbin], /double)

              outstruct.rsum[radbin,angbin] = total(instruct[w].rsum[radbin,angbin], /double)
              IF doangsum THEN BEGIN 
                  outstruct.angsum[radbin,angbin] = $
                    total(instruct[w].angsum[radbin,angbin], /double)
                  outstruct.meanang[radbin,angbin] = $
                    outstruct.angsum[radbin,angbin]/outstruct.npair[radbin,angbin]
              ENDIF 

              ;; means
              outstruct.meanr[radbin,angbin]    = $
                outstruct.rsum[radbin,angbin]/outstruct.npair[radbin,angbin]

              outstruct.rmax_act[radbin,angbin] = max(instruct[w].rmax_act[radbin,angbin])


              wmom, instruct[w].sigma[radbin,angbin], instruct[w].sigmaerr[radbin,angbin], $
                    wmean, wsig, werr

              outstruct.sigma[radbin,angbin] = wmean
              outstruct.sigmaerr[radbin,angbin] = werr

              wmom, instruct[w].orthosig[radbin,angbin], instruct[w].orthosigerr[radbin,angbin], $
                    owmean, owsig, owerr

              outstruct.orthosig[radbin,angbin] = owmean
              outstruct.orthosigerr[radbin,angbin] = owerr

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; cumulative means
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              ;; meanr
              outstruct.tmeanr[radbin,angbin] = $
                total(outstruct.rsum[0:radbin,angbin], /double)/total_int(outstruct.npair[0:radbin,angbin])

              ;; sigma
              wmom, outstruct.sigma[0:radbin,angbin], outstruct.sigmaerr[0:radbin,angbin], $
                    wmean, wsig, werr

              outstruct.tsigma[radbin,angbin] = wmean
                  outstruct.tsigmaerr[radbin,angbin] = werr

              ;; orthosig
              wmom, outstruct.orthosig[0:radbin,angbin], outstruct.orthosigerr[0:radbin,angbin], $
                    owmean, owsig, owerr

              outstruct.torthosig[radbin,angbin] = owmean
              outstruct.torthosigerr[radbin,angbin] = owerr


              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; average the clustering correction
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              IF docorr THEN BEGIN 
                  wmom, instruct[w].clust_corr[radbin,angbin], $
                        instruct[w].clust_corr_err[radbin,angbin], $
                        wmean, wsig, werr
                  outstruct.clust_corr[radbin,angbin] = wmean
                  outstruct.clust_corr_err[radbin,angbin] = werr
              ENDIF 

          ENDIF ELSE BEGIN 
              message,'None passed err cuts for bin '+ntostr(radbin),/inf
              ;;wait,10
          ENDELSE 
                                          
      ENDFOR 
  ENDFOR 



  return,outstruct

END 
