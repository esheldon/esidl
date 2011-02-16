PRO source_vs_rpetro,scat,nzstruct

  ;; make various plots binned versus rpetro

  ind = nzstruct.useind
  rpet = scat[ind].rpetro

  binsize = 0.1
  
  histo=histogram(rpet,bin=binsize,reverse_indices=rev_ind,min=min,max=max)
  nbins=n_elements(histo)

  rpetro = fltarr(nbins)
  rsmear = rpetro
  rsmear_err = rpetro
  phtz = rpetro
  phtz_err = rpetro
  err = rpetro
  err_err = rpetro

  FOR binnum=0L, nbins-1 DO BEGIN 
      IF rev_ind[binnum] NE rev_ind[binnum+1] THEN BEGIN 
          w=rev_ind[ rev_ind[binnum]:rev_ind[binnum+1]-1 ]
          
          nw=n_elements(w)
          IF nw GE 2 THEN BEGIN
              rpetro[binnum] = mean(rpet[w])

              ;; smear polarizability
              yi = scat[ind[w]].rsmear
              result=moment( yi , maxmoment=2, sdev=sdev)
              rsmear[binnum] = result[0]
              rsmear_err[binnum] = sdev/sqrt(nw)

              ;; photoz
              errsend = sqrt( scat[ind[w]].photoz_zerr^2 + 0.01^2 )
              wmom, scat[ind[w]].photoz_z, errsend, wmean, wsig, werr
              phtz[binnum] = wmean
              phtz_err[binnum] = werr

              ;; error
              yi = sqrt(scat[ind[w]].e1e1err^2 + scat[ind[w]].e2e2err^2)
              result=moment( yi , maxmoment=2, sdev=sdev)
              err[binnum] = result[0]
              err_err[binnum] = sdev/sqrt(nw)

          ENDIF 
      ENDIF 

  ENDFOR 

  xtitle = 'r'

  IF !d.name EQ 'PS' THEN !p.charsize = 1.2

  erase & multiplot, [1,3], /square
  plot, rpetro, rsmear, psym=8, ytitle='R!Dsmear!N', $
             ystyle=2, xstyle=2
  
  multiplot & plot, rpetro, phtz, psym=8, ytitle='Photoz', $
             ystyle=2, xstyle=1, yrange=[0,0.45]
  multiplot &  plot, rpetro, err, psym=8, ytitle=!csym.sigma+'(e)', xtitle=xtitle, $
             ystyle=2, xstyle=1, yrange=[0,0.45]
  oplot, rpetro, histo*0.1/max(histo), psym=10
  multiplot,/reset

END 
