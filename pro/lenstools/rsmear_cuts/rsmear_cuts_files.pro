PRO rsmear_cuts_files, run, rerun, purity, clr, $
                       fitfile, psfile, comb_psfile, comb_fitfile, $
                       camcol=camcol, outdir=outdir, hirata=hirata

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: rsmear_cuts_files, run, rerun, purity, clr, fitfile, psfile, comb_psfile, comb_fitfile, camcol=camcol, outdir=outdir, hirata=hirata'
      return
  ENDIF 

  IF keyword_set(hirata) THEN rstr = '_h' ELSE rstr=''
  bstr = '_bayes'

  IF n_elements(camcol) NE 0 THEN BEGIN 

      IF n_elements(outdir) EQ 0 THEN BEGIN 
          fetch_dir,run,camcol,rerun,tdir,atldir,file,atlfile,$
                    corrdir=corrdir, corratldir=outdir
      ENDIF 
      psfile = outdir + 'run'+run2string(run)+$
        '_col'+ntostr(camcol)+$
        '_'+!colors[clr]+$
        bstr+$
        rstr+$
        '_purity'+ntostr(purity,4,/round)+'.ps'
      fitfile = repstr(psfile, '.ps','.fit')
  ENDIF ELSE BEGIN
      tcamcol = 1
      IF n_elements(outdir) EQ 0 THEN BEGIN 
          fetch_dir,run,tcamcol,rerun, tdir,atldir,file,atlfile,$
                    corrdir=corrdir, corratldir=outdir
      ENDIF 
      psfile = outdir + 'run'+run2string(run)+$
        '_allcols'+$
        '_'+!colors[clr]+$
        bstr+rstr+'_purity'+ntostr(purity,4,/round)+'.ps'
      fitfile = repstr(psfile, '.ps','.fit')
      comb_psfile = repstr(psfile, 'allcols', 'comb_allcols')
      comb_fitfile = repstr(fitfile, 'allcols', 'comb_allcols')
  ENDELSE 

END 
