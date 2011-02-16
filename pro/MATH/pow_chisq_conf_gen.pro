PRO pow_chisq_conf_gen, datax, data, covariance, powrange, normrange, $
                        npow, nnorm, $
                        chisq_surf, $
                        bestp, bestn, $
                        powlow, powhigh, $
                        normlow, normhigh, $
                        perrlow=perrlow,perrhigh=perrhigh,$
                        nerrlow=nerrlow,nerrhigh=nerrhigh, $
                        $
                        powallow = powallow, normallow=normallow, $
                        powvals=powvals, normvals=normvals,$
                        likelihood=likelihood, $
                        pow_like=pow_like, norm_like=norm_like, $
                        $
                        wuse=wuse,$
                        $
                        diagonals=diagonals,$
                        $
                        yfit=yfit,$
                        yallow_low=yallow_low, yallow_high=yallow_high,$
                        minchisq=minchisq,$
                        degfree=degfree,$
                        chisq_diff=chisq_diff, $
                        $
                        getpath=getpath, info=info, xy=xy, $
                        $
                        nodisplay=nodisplay, $
                        aspect=aspect, center=center, $
                        xtitle=xtitle, ytitle=ytitle, $
                        noplotmin=noplotmin,$
                        plot_both=plot_both, $
                        xtick_get=xtick_get, ytick_get=ytick_get, $
                        dolegend=dolegend, names=names,nkeep=nkeep, $
                        _extra=extra
  
    IF n_params() LT 7 THEN BEGIN
      print,'-Syntax: pow_chisq_conf_gen, datax, data, covariance, powrange, normrange, npow, nnorm, $'
      print,'       chisq_surf, $'
      print,'       bestp, bestn, $'
      print,'       powlow, powhigh, $'
      print,'       normlow, normhigh, $'
      print,'       perrlow=perrlow,perrhigh=perrhigh,$'
      print,'       nerrlow=nerrlow,nerrhigh=nerrhigh, $'
      print,'       $'
      print,'       likelihood=likelihood, $'
      print,'       pow_like=pow_like, norm_like=norm_like, $'
      print,'       $'
      print,'       wuse=wuse,$'
      print,'       $'
      print,'       diagonals=diagonals,$'
      print,'       $'
      print,'       yfit=yfit,$'
      print,'       yallow_low=yallow_low, yallow_high=yallow_high,$'
      print,'       minchisq=minchisq,$'
      print,'       degfree=degfree,$'
      print,'       $'
      print,'       chisq_diff=chisq_diff, $'
      print,'       $'
      print,'       powallow = powallow, normallow=normallow, $'
      print,'       powvals = powvals, normvals=normvals, $'
      print,'       $'
      print,'       getpath=getpath, info=info, xy=xy, $'
      print,'       $'
      print,'       nodisplay=nodisplay, $'
      print,'       aspect=aspect, center=center, $'
      print,'       noplotmin=noplotmin,$'
      print,'       xtick_get=xtick_get, ytick_get=ytick_get, $'
      print,'       dolegend=dolegend, names=names, nkeep=nkeep, $'
      print,'       project=project, $'
      print,'       _extra=extra'
      print
      return
    ENDIF 

  ;; generate the values
  powvals = arrscl( findgen(npow), powrange[0], powrange[1] )
  normvals = arrscl( findgen(nnorm), normrange[0], normrange[1] )

  ;; call pow_chisq_conf
  pow_chisq_conf, datax, data, covariance, powvals, normvals, $
                  chisq_surf, $
                  bestp, bestn, $
                  powlow, powhigh, $
                  normlow, normhigh, $
                  perrlow=perrlow,perrhigh=perrhigh,$
                  nerrlow=nerrlow,nerrhigh=nerrhigh, $
                  $
                  likelihood=likelihood, $
                  pow_like=pow_like, norm_like=norm_like, $
                  $
                  wuse=wuse,$
                  $
                  diagonals=diagonals,$
                  $
                  yfit=yfit,$
                  yallow_low=yallow_low, yallow_high=yallow_high,$
                  minchisq=minchisq,$
                  degfree=degfree,$
                  $
                  chisq_diff=chisq_diff, $
                  $
                  powallow = powallow, normallow=normallow, $
                  $
                  getpath=getpath, info=info, xy=xy, $
                  $
                  nodisplay=nodisplay, $
                  aspect=aspect, center=center, $
                  xtitle=xtitle, ytitle=ytitle, $
                  noplotmin=noplotmin,$
                  xtick_get=xtick_get, ytick_get=ytick_get, $
                  dolegend=dolegend, names=names, nkeep=nkeep, $
                  project=project, $
                  _extra=extra

  return
END 

