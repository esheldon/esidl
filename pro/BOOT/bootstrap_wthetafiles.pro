PRO bootstrap_wthetafiles, file, rfile, nresample, radrange=radrange

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: bootstrap_wthetafiles, file, rfile, nresample, radrange=radrange'
      return
  ENDIF 

  ;; when more stripes, will combine files 

  outfile = repstr(file, 'lumlensum', 'bootstrap')
  ttfile = repstr(file, 'lumlensum_', '')
  print
  print,'Outputting file: ',outfile
  print

  tt=mrdfits(ttfile, 1, hdr)
  nrad=n_elements(tt.meanr)
  IF n_elements(radrange) NE 0 THEN BEGIN 
      wuse=where( (tt.meanr LT radrange[1]) AND (tt.meanr GT radrange[0]), nw)
      IF nw EQ 0 THEN message,'No bins passed radrange cut!'
  ENDIF ELSE wuse=lindgen(nrad)

  ;; read files and strip down for memory saving
  t=mrdfits(file, 1)
  s=create_struct('zindex', 0L, $
                  'wsum', fltarr(nrad), $
                  'lsum', fltarr(nrad) )

  nlens=n_elements(t)
  ts=replicate(s, nlens)
  ts.zindex = t.zindex
  ts.wsum = t.wsum
  ts.lsum = t.lsum
  delvarx, t

  ;; now rand
  tr=mrdfits(rfile, 1)

  nrand=n_elements(tr)
  trs=replicate(s, nrand)
  trs.zindex = tr.zindex
  trs.wsum = tr.wsum
  trs.lsum = tr.lsum
  delvarx, tr

  lumdensboot, ts, trs, nresample, meanlum, meanlumerr, covariance, wuse=wuse

  ss=create_struct('nresample', nresample,$
                   'meanr', tt.meanr[wuse],$
                   'meanlum', meanlum, $
                   'meanlumerr', meanlumerr, $
                   'covariance', covariance)
  print
  print,'Outputting file: ',outfile
  print
  mwrfits, ss, outfile,/create
  


END 
