PRO plot_envgmcf, color=color

  errfac=1.32

  indir='/sdss5/data0/lensout/stripe10/'
  outdir=indir

  IF keyword_set(color) THEN outfile = outdir+'environment_denscont_color.ps' $
  ELSE outfile = outdir+'environment_denscont.ps'
  begplot,name=outfile, color=color

  lfile=indir+'low_zgal_gal_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_comb_corr_N1.fit'
  hfile=indir+'high_zgal_gal_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_comb_corr_N1.fit'

  hcolor = !red
  lcolor = !blue

  low = mrdfits(lfile,1)
  high = mrdfits(hfile,1)

  low.sigmaerr = low.sigmaerr*errfac
  high.sigmaerr = high.sigmaerr*errfac

  yrange=prange(low.sigma,high.sigma,low.sigmaerr,high.sigmaerr)

  ytitle = '!S'+!tsym.sigma_cap+'!R!A'+!tsym.minus+'!N ('+!tsym.ltequal+'R) '+!tsym.minus+' !S'+$
    !tsym.sigma_cap+'!R!A'+!tsym.minus+$
    '!N (R) (h M'+sunsymbol()+' pc!U'+!tsym.minus+'2!N)!X'
  xtitle='Projected Radius (h!U'+!tsym.minus+'1!N kpc)'

  aploterror, !gratio, high.meanr, high.sigma, high.sigmaerr, $
    line=2, xtitle=xtitle, ytitle=ytitle, yrange=yrange
  oplot,[0,10000], [0,0]
  IF NOT keyword_set(color) THEN BEGIN 
      oploterror, low.meanr, low.sigma, low.sigmaerr, $
        line=0
      legend, ['High Density Regions','Low Density Regions'],$
        line=[2,0],/right,thick=[!p.thick,!p.thick]
  ENDIF ELSE BEGIN 
      oploterror, high.meanr, high.sigma, high.sigmaerr, $
        line=2, color=hcolor, errcolor=hcolor
      oploterror, low.meanr, low.sigma, low.sigmaerr, $
        line=0, color=lcolor, errcolor=lcolor
      legend, ['High Density Regions','Low Density Regions'],$
        line=[2,0],/right,thick=[!p.thick,!p.thick], $
        colors=[hcolor, lcolor]
  ENDELSE 
  endplot

END 
