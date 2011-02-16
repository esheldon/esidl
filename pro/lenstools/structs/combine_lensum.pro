
FUNCTION combine_lensum, lensum, hdrin=hdrin, silent=silent, index=index

    if n_params() lt 1 then begin 
        print,'-Syntax: sh = combine_lensum(lensum, hdr=, /silent, index=index)'
        print,'index is the ones you want to combine. Default is all'
        print,'hdr is the header you read from the input file'
        return,-1
    endif 

    nlensum = n_elements(lensum)
    IF n_elements(index) EQ 0 THEN index = lindgen(nlensum)

    nlens = n_elements(index)



    s=size(lensum[0].rsum)
    ;; assume 2 or fewer dimensions
    IF s[0] LT 2 THEN nang = 1 ELSE nang = s[2]
    nrad = s[1]

    ;; Should we use the outer bin?
    IF (total(lensum[index].npair[nrad-1,*],/int) LT $
        total(lensum[index].npair[nrad-2,*],/int) ) THEN BEGIN 
        arrval = lensum[0].rsum[0:nrad-2,*]
        nrad = nrad-1
    ENDIF ELSE BEGIN
        arrval = lensum[0].rsum
    ENDELSE 

    arrval = double(arrval)
    arrval[*] = 0d


    shst = shstruct(arrval)


    IF n_elements(hdrin) NE 0 THEN BEGIN 

        IF size(hdrin, /tname) EQ 'STRUCT' THEN BEGIN 
            hdr = idlstruct_hclean(hdrin)
            header_tags = tag_names(hdr)
            sh_tags = tag_names(shst)

            ;; either copy the header tag or add it
            FOR i=0L, n_elements(header_tags)-1 DO BEGIN 
                IF tag_exist(shst, header_tags[i], index=itag) THEN BEGIN 
                    shst.(itag) = hdr.(i)
                ENDIF ELSE BEGIN 
                    shst = create_struct(header_tags[i], hdr.(i), shst)
                ENDELSE 
            ENDFOR 
        ENDIF 

    ENDIF 



    shst.nlenses = nlens
    IF tag_exist(shst, 'nlens') THEN shst.nlens = nlens

    ;; lens weight
    lensw = lensum[index].weight
    lenswsum = total( lensw ,/double)
    shst.lenswsum = lenswsum

    ;; scalar sums/means
    IF tag_exist(lensum[0], 'totpairs') THEN BEGIN 
        shst.totpairs = total( lensum[index].totpairs,/int )
    ENDIF ELSE BEGIN 
        shst.totpairs = total( lensum[index].tot_pairs,/int )
    ENDELSE 
    shst.wsum_ssh = total(lensum[index].wsum_ssh, /double)
    shst.sshsum = total(lensum[index].sshsum, /double)
    shst.ssh = shst.sshsum/shst.wsum_ssh


    IF tag_exist(lensum[0], 'wscritinvsum') THEN BEGIN 
        DoScritinvMean = 1
    ENDIF ELSE BEGIN 
        DoScritinvMean = 0
    ENDELSE 

    IF tag_exist(lensum[0], 'angsum') THEN doangsum=1 ELSE doangsum=0

    FOR radbin=0L, nrad-1 DO BEGIN 
        FOR angbin = 0L, nang-1 DO BEGIN 

            terr = lensum[index].sigmaerr[radbin,angbin]
            oterr = lensum[index].orthosigerr[radbin,angbin]
            w = where(terr GT 0.0 AND finite(terr) AND $
                oterr GT 0.0 AND finite(oterr) AND $
                abs(lensum[index].sigma[radbin,angbin]) LT 1000000, nw)
            shst.nlbin[radbin,angbin] = nw

            IF nw NE 0 THEN BEGIN 

                w = index[w]

                ;; sums
                shst.npair[radbin,angbin] = $
                    total(lensum[w].npair[radbin,angbin],/int)
                shst.tnpair[radbin,angbin] = $
                    total(shst.npair[0:radbin,angbin],/int)

                IF shst.npair[radbin,angbin] EQ 0 THEN BEGIN 
                    message,'npair = 0 when it should not'
                ENDIF 

                wsum = total(lensum[w].wsum[radbin,angbin], /double)
                wsum2 = total(lensum[w].wsum[radbin,angbin]^2, /double)

                shst.wsum[radbin,angbin] = wsum
                shst.wsum2[radbin,angbin] = wsum2

                wsum_mean =  wsum/nlens
                wsum_err = sqrt(wsum2/nlens - wsum_mean^2)/sqrt(nlens)

                shst.wsum_mean[radbin,angbin] = wsum_mean
                shst.wsum_err[radbin,angbin] = wsum_err

                shst.sigerrsum[radbin,angbin] = $
                    total(lensum[w].sigerrsum[radbin,angbin], /double)
                shst.orthosigerrsum[radbin,angbin] = $
                    total(lensum[w].orthosigerrsum[radbin,angbin], /double)

                shst.rsum[radbin,angbin] = total(lensum[w].rsum[radbin,angbin], /double)
                IF doangsum THEN BEGIN 
                    shst.angsum[radbin,angbin] = $
                        total(lensum[w].angsum[radbin,angbin], /double)
                    shst.meanang[radbin,angbin] = $
                        shst.angsum[radbin,angbin]/shst.npair[radbin,angbin]
                ENDIF 

                ;; means
                shst.meanr[radbin,angbin]    = $
                    shst.rsum[radbin,angbin]/shst.npair[radbin,angbin]

                shst.rmax_act[radbin,angbin] = max(lensum[w].rmax_act[radbin,angbin])
                shst.rmin_act[radbin,angbin] = min(lensum[w].rmin_act[radbin,angbin])
                shst.area_act[radbin,angbin] = $
                    !dpi*( shst.rmax_act[radbin,angbin]^2 - shst.rmin_act[radbin,angbin]^2 )

                IF nw GT 10 THEN calcerr=1 ELSE calcerr=0

                wmom, $
                    lensum[w].sigma[radbin,angbin], $
                    lensum[w].sigmaerr[radbin,angbin], $
                    wmean, wsig, werr, calcerr=calcerr

                shst.sigma[radbin,angbin] = wmean
                shst.sigmaerr[radbin,angbin] = werr
                ; alternative error estimate just from ellipticities
                wmom, $
                    lensum[w].sigma[radbin,angbin], $
                    lensum[w].sigmaerr[radbin,angbin], $
                    wmean, wsig, werr, calcerr=0

                shst.sigmaerr2[radbin,angbin] = werr
 

                ; ortho
                wmom, $
                    lensum[w].orthosig[radbin,angbin], $
                    lensum[w].orthosigerr[radbin,angbin], $
                    owmean, owsig, owerr, calcerr=calcerr



                shst.orthosig[radbin,angbin] = owmean
                shst.orthosigerr[radbin,angbin] = owerr

                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                ;; cumulative means
                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                shst.tmeanr[radbin,angbin] = $
                    total(shst.rsum[0:radbin,angbin],/double)/total(shst.npair[0:radbin,angbin],/int)

                wmom, $
                    shst.sigma[0:radbin,angbin],$
                    shst.sigmaerr[0:radbin,angbin], $
                    wmean, wsig, werr

                shst.tsigma[radbin,angbin] = wmean
                shst.tsigmaerr[radbin,angbin] = werr

                wmom, $
                    shst.orthosig[0:radbin,angbin], $
                    shst.orthosigerr[0:radbin,angbin], $
                    owmean, owsig, owerr

                shst.torthosig[radbin,angbin] = owmean
                shst.torthosigerr[radbin,angbin] = owerr

                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                ;; Mean inverse critical density and sum
                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                IF DoScritinvMean THEN BEGIN 
                    shst.wscritinvsum[radbin,angbin] = $
                        total( lensum[w].wscritinvsum[radbin,angbin], /double )
                  shst.scritinv[radbin,angbin] = $
                    shst.wscritinvsum[radbin,angbin]/wsum
                ENDIF 


                IF NOT keyword_set(silent) THEN BEGIN 
                    mess = 'nuse/nlens = '+ntostr(nw)+'/'+ntostr(nlens)+$
                        '  radius = '+ntostr(shst.meanr[radbin,angbin])
                    print,mess
                ENDIF 

            ENDIF 
                                          
        ENDFOR 
    ENDFOR 

    return,shst

END 
