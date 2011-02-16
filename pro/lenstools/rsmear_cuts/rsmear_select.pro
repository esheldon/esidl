PRO rsmear_select, run, rerun, purity, clr, rmag, Rsmear, keep, nkeep, $
                   camcol=camcol, maxrsmear=maxrsmear, rcuts=rcuts

  IF n_params() LT 7 THEN BEGIN  
      print,'-Syntax: rsmear_select, run, rerun, purity, clr, rmag, Rsmear, keep, nkeep, $'
      print,'         camcol=camcol, maxrsmear=maxrsmear, rcuts=rcuts'
      print,'Rsmear should be in the "clr" bandpass. rmag is reddening corrected petrosian magnitude.'
      return
  ENDIF 

  ;; hard cut
  IF n_elements(maxrsmear) EQ 0 THEN maxrsmear = 0.8

  maxrmag = 22.0
  minrmag = 18.0

  rsmear_cuts_files, run, rerun, purity, clr, fitfile, camcol=camcol

  print,'Reading file: ',fitfile
  cstr = mrdfits(fitfile, 1)
  
  ;; only use for rmag < 22.0
  nn = n_elements(rsmear)
  rcuts = fltarr(nn)
  keep = where(rmag LE maxrmag AND rmag GE minrmag, nkeep)
  IF nkeep NE 0 THEN BEGIN 
      ;; interpolate the Rsmear cut
      rcuts[keep] = interpol(cstr.rsmear_cuts, cstr.meanmag, rmag[keep]) > 0.0 < maxrsmear

      ;; make the cuts
      keep2 = where(Rsmear[keep] LT rcuts[keep] AND Rsmear [keep] GT 0.0, nkeep)

      ;; did any pass the cuts?
      IF nkeep NE 0 THEN keep = keep[keep2] ELSE keep=-1
  ENDIF
  
;  myusersym, 'fill_circle'

;  aplot, 1, [0], /nodata, xrange=[18,22], yrange=[-0.1,1]
;  oplot, cstr.meanmag, cstr.rsmear_cuts, psym=8
;  s = sort(rmag[keep])
;  oplot, rmag[keep[s]], rcuts[keep[s]], color=!green



END 
