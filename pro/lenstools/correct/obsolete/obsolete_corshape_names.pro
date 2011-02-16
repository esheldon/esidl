PRO corshape_names, run, rerun, camcol, clr, psfile, fitsfile, $
                    extno=extno, allcols=allcols, hirata=hirata

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: corshape_names, run, rerun, camcol, clr, psfile, fitsfile, '
      print,'          extno=extno, allcols=allcols'
      print,'If /allcols camcol is ignored'
      return
  ENDIF 

  IF n_elements(ext) EQ 0 THEN extno=1
  IF keyword_set(allcols) THEN BEGIN
      tcamcol = 1
      addstr='_allcols_' 
  ENDIF ELSE BEGIN
      addstr='_'+ntostr(camcol)+'_'
      tcamcol = camcol
  ENDELSE 

  ext='N'+ntostr(extno)+'.ps'
  IF keyword_set(hirata) THEN ext = 'h_'+ext
  fetch_dir, run, tcamcol, rerun, dir, atldir, corrdir=corrdir,$
    corratldir = fitdir
  psfile = fitdir+'corshape_'+run2string(run)+addstr+!colors[clr]+'_'+ext

  ext='N'+ntostr(extno)+'.fit'
  IF keyword_set(hirata) THEN ext = 'h_'+ext
  fitsfile = fitdir+'corshape_'+run2string(run)+addstr+!colors[clr]+'_'+ext

END 
