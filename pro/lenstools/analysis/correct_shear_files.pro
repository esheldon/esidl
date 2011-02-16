PRO correct_shear_files, files, rfiles, corrfiles, wsubtract=wsubtract

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: correct_shear_files, files, rfiles [, corrfiles]'
      return
  ENDIF 

  nf = n_elements(files)
  nrf = n_elements(rfiles)

  IF nf NE nrf THEN message,'files and rfiles must be same size'

  corrfiles = strarr(nf)

  FOR i=0L, nf-1 DO BEGIN 

      fs = str_sep(files[i], '_N')

      IF n_elements(fs) NE 2 THEN message,'Filename format error: '+files[i]

      corrfiles[i] = fs[0]+'_corr_N'+fs[1]

      print,'Reading Shear file: ',files[i]
      t = mrdfits(files[i], 1, hdr, /silent)
      print,'Reading Rand file: ',rfiles[i]
      r = mrdfits(rfiles[i], 1, /silent)

      correct_shear, t, r, corr, corr_err

      ;; Might only want to subtract random for
      ;; radii with small errors in the randoms
      IF n_elements(wsubtract) NE 0 THEN BEGIN 
          t.sigma[wsubtract] = t.sigma[wsubtract] - r.sigma[wsubtract]
          t.sigmaerr[wsubtract] = sqrt(t.sigmaerr[wsubtract]^2 + r.sigmaerr[wsubtract]^2)
      ENDIF ELSE BEGIN 
          w=where(t.npair GT 0, nw)
          IF nw NE 0 THEN BEGIN 
              t.sigma[w] = t.sigma[w] - r.sigma[w]
              t.sigmaerr[w] = sqrt(t.sigmaerr[w]^2 + r.sigmaerr[w]^2)
          ENDIF 
      ENDELSE 


      nr = n_elements(corr)
      vals = replicate('fltarr('+ntostr(nr)+')', 2)
      add_tags, t, ['clust_corr', 'clust_corr_err'], vals, newt

      newt.clust_corr = corr
      newt.clust_corr_err = corr_err

      print,'Outputting corr file: ',corrfiles[i]
      fxhclean, hdr
      mwrfits, newt, corrfiles[i], hdr, /create

  ENDFOR 

END 
