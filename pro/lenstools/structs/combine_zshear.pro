PRO combine_zshear, files, outstruct, comoving=comoving, already_comoving=already_comoving

  IF n_params() LT 1 THEN BEGIN
      print,'-Syntax: combine_zshear, files, outstruct, comoving=comoving, already_comoving=already_comoving'
      return
  ENDIF 

  ;; /comoving: if its in physical or there is not hdr keyword for comoving,
  ;; then convert it

  nf=n_elements(files)
  ac = 0
  FOR i=0, nf-1 DO BEGIN

      print,'Reading file: ',files[i]
      struct1 = mrdfits(files[i], 1, hdr, /silent)

      ;; check if the comoving hdr keyword exists.  If not, it is 
      ;; assume it is physical and convert to co-moving

      IF keyword_set(comoving) THEN BEGIN 

          cm = sxpar(hdr, 'comoving', count=cmexist)

          docomv = 0
          IF cmexist THEN BEGIN 
              IF cm EQ 'no' THEN docomv = 1 ELSE ac = ac+1
          ENDIF ELSE docomv = 1

          IF docomv THEN BEGIN 
              print,'Converting to comoving'
              convert_lensout_comoving, struct1
          ENDIF 
      ENDIF 

      IF i EQ 0 THEN BEGIN
          instruct = temporary(struct1)     
      ENDIF ELSE BEGIN 
          concat_structs, temporary(struct1), temporary(instruct), tmpstruct
          instruct = temporary(tmpstruct)
      ENDELSE 
  ENDFOR 

  IF keyword_set(comoving) THEN BEGIN 
      IF ac EQ nf THEN BEGIN 
          already_comoving=1 
      ENDIF ELSE IF ac EQ 0 THEN BEGIN 
          already_comoving = 0
      ENDIF ELSE BEGIN 
          message,'Some files were co-moving but some were not!'
      ENDELSE 
  ENDIF 

  IF tag_exist(instruct[0], 'compcut') THEN BEGIN 
      IF n_elements(rem_dup(instruct.compcut)) NE 1 THEN message,'more than one compuct here!!'
  ENDIF 

  ;; now combine
  outstruct = instruct[0]

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

  outstruct.scritinvsum = total(instruct.scritinvsum, /double)
  outstruct.meanscritinv = outstruct.scritinvsum/lenswsum

  IF tag_exist(instruct[0], 'rmin_act') THEN DO_rminact=1 ELSE DO_rminact=0

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
;              wsum2 = total(instruct[w].wsum[radbin,angbin]^2, /double)
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
              IF DO_rminact THEN BEGIN
                  wt=where(instruct[w].rmin_act[radbin,angbin] GT 0.)
                  outstruct.rmin_act[radbin,angbin] = $
                    min(instruct[w[wt]].rmin_act[radbin,angbin])

                  R1 = outstruct.rmin_act[radbin,angbin]
                  R2 = outstruct.rmax_act[radbin,angbin]

              ENDIF ELSE BEGIN 
                  R1 = outstruct.rmin + radbin*outstruct.binsize
                  R2 = outstruct.rmin + (radbin+1)*outstruct.binsize
              ENDELSE 

              outstruct.area[radbin,angbin] = !pi*(R2^2 - R1^2)
              outstruct.tarea[radbin,angbin] = !pi*(R2^2 - outstruct.rmin^2)
              IF outstruct.area[radbin,angbin] NE 0. THEN BEGIN 
                  outstruct.density[radbin,angbin] = $
                    outstruct.npair[radbin,angbin]/outstruct.area[radbin,angbin]/outstruct.nlenses
              ENDIF 

              wmom, instruct[w].sigma[radbin,angbin], instruct[w].sigmaerr[radbin,angbin], $
                    wmean, wsig, werr
;              IF nf GE 2 THEN BEGIN 
;                  wmom, instruct[w].sigma[radbin,angbin], instruct[w].sigmaerr[radbin,angbin], $
;                        wmean2, wsig2, werr2, /calcerr
;                  werr = max([werr,werr2])
;              ENDIF
              outstruct.sigma[radbin,angbin] = wmean
              outstruct.sigmaerr[radbin,angbin] = werr

              wmom, instruct[w].orthosig[radbin,angbin], instruct[w].orthosigerr[radbin,angbin], $
                    owmean, owsig, owerr
;              IF nf GE 2 THEN BEGIN 
;                  wmom, instruct[w].orthosig[radbin,angbin], instruct[w].orthosigerr[radbin,angbin], $
;                        owmean2, owsig2, owerr2, /calcerr
;                  owerr = max([owerr,owerr2])
;              ENDIF
              outstruct.orthosig[radbin,angbin] = owmean
              outstruct.orthosigerr[radbin,angbin] = owerr

              ;; alternative error estimate
              outstruct.sigmaerr2[radbin,angbin] = $
                sqrt(total(instruct[w].sigerrsum[radbin,angbin], /double))/outstruct.wsum[radbin,angbin]

              outstruct.orthosigerr2[radbin,angbin] = $
                sqrt(total(instruct[w].orthosigerrsum[radbin,angbin], /double))/outstruct.owsum[radbin,angbin]

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

              ;; twsum, twsum_mean, etc.
;              outstruct.twsum[radbin,angbin] = total(outstruct.wsum[0:radbin,angbin])

;              wmom, outstruct.wsum_mean[0:radbin,angbin], outstruct.wsum_err[0:radbin,angbin], $
;                    wmean, wsig, werr

;              outstruct.twsum_mean[radbin,angbin] = wmean
;              outstruct.twsum_err[radbin,angbin] = werr

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; alternative error estimate
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              outstruct.tsigmaerr2[radbin,angbin] = $
                sqrt(total(outstruct.sigerrsum[0:radbin,angbin], /double))/total(outstruct.wsum[0:radbin,angbin], /double)

              outstruct.torthosigerr2[radbin,angbin] = $
                sqrt(total(outstruct.orthosigerrsum[0:radbin,angbin], /double))/total(outstruct.owsum[0:radbin,angbin], /double)

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



  return
END 
