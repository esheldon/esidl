
PRO egal_vs_epsf_plotfits, fitstruct, rvals=rvals

  IF fitstruct.hirata THEN BEGIN 
      xtitle = 'R_h['+!colors[fitstruct.bandpass]+']'
  ENDIF ELSE BEGIN 
      xtitle = 'R['+!colors[fitstruct.bandpass]+']'
  ENDELSE 

  nsmooth = 5

  fitst = fitstruct.fitst

  erase & multiplot,[1,2]
  ploterror,$
    fitst.mean_r,fitst.slope_ge1pe1,fitst.slope_ge1pe1_error,psym=8, $
    ytitle = 'slope', title = 'gal e1 vs psf e1'
  oplot, fitst.mean_r, smooth(fitst.slope_ge1pe1, nsmooth, /edge_truncate)
  multiplot

  IF keyword_set(rvals) THEN BEGIN 
      rbin = 0.02
      cyrange = !y.crange
      plothist, rvals, xhist, yhist, bin=rbin, /noplot
      nyhist = yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
      oplot, xhist, nyhist, psym=10, color=!grey50
  ENDIF 


  ploterror,$
    fitst.mean_r,fitst.intercept_ge1pe1,fitst.intercept_ge1pe1_error,psym=8,$
    ytitle = 'intercept', xtitle = xtitle
  oplot, fitst.mean_r, smooth(fitst.intercept_ge1pe1, nsmooth, /edge_truncate)
  multiplot,/reset



  key=prompt_kbrd("Hit a key for the next plot")

  erase & multiplot,[1,2]
  ploterror,$
    fitst.mean_r,fitst.slope_ge1pe2,fitst.slope_ge1pe2_error,psym=8, $
    ytitle = 'slope', title = 'gal e1 vs psf e2'
  oplot, fitst.mean_r, smooth(fitst.slope_ge1pe2, nsmooth, /edge_truncate)

  multiplot
  ploterror,$
    fitst.mean_r,fitst.intercept_ge1pe2,fitst.intercept_ge1pe2_error,psym=8,$
    ytitle = 'intercept', xtitle = xtitle
  oplot, fitst.mean_r, smooth(fitst.intercept_ge1pe2, nsmooth, /edge_truncate)
  multiplot,/reset



  key=prompt_kbrd("Hit a key for the next plot")

  erase & multiplot,[1,2]
  ploterror,$
    fitst.mean_r,fitst.slope_ge2pe1,fitst.slope_ge2pe1_error,psym=8, $
    ytitle = 'slope', title = 'gal e2 vs psf e1'
  oplot, fitst.mean_r, smooth(fitst.slope_ge2pe1, nsmooth, /edge_truncate)
  multiplot
  ploterror,$
    fitst.mean_r,fitst.intercept_ge2pe1,fitst.intercept_ge2pe1_error,psym=8,$
    ytitle = 'intercept', xtitle = xtitle
  oplot, fitst.mean_r, smooth(fitst.intercept_ge2pe1, nsmooth, /edge_truncate)
  multiplot,/reset



  key=prompt_kbrd("Hit a key for the next plot")

  erase & multiplot,[1,2]
  ploterror,$
    fitst.mean_r,fitst.slope_ge2pe2,fitst.slope_ge2pe2_error,psym=8, $
    ytitle = 'slope', title = 'gal e2 vs psf e2'
  oplot, fitst.mean_r, smooth(fitst.slope_ge2pe2, nsmooth, /edge_truncate)

  multiplot
  ploterror,$
    fitst.mean_r,fitst.intercept_ge2pe2,fitst.intercept_ge2pe2_error,psym=8,$
    ytitle = 'intercept', xtitle = xtitle
  oplot, fitst.mean_r, smooth(fitst.intercept_ge2pe2, nsmooth, /edge_truncate)
  multiplot,/reset




END 

FUNCTION egal_vs_epsf_fitstruct

  struct = $
    { $
      min_r: 0.0, $
      max_r: 0.0, $
      mean_r: 0.0, $
      label: '', $
      intercept_ge1pe1:0d, $
      intercept_ge1pe1_error:0d, $
      slope_ge1pe1:0d, $
      slope_ge1pe1_error:0d,$
      $
      intercept_ge2pe1:0d, $
      intercept_ge2pe1_error:0d, $
      slope_ge2pe1:0d, $
      slope_ge2pe1_error:0d,$
      $
      intercept_ge1pe2:0d, $
      intercept_ge1pe2_error:0d, $
      slope_ge1pe2:0d, $
      slope_ge1pe2_error:0d,$
      $
      intercept_ge2pe2:0d, $
      intercept_ge2pe2_error:0d, $
      slope_ge2pe2:0d, $
      slope_ge2pe2_error:0d $
    }

  return,struct
END 

PRO egal_vs_epsf_plot_fitlin, struct, w, nperEbin, title, fitst, $
                      hirata=hirata, e1c=e1c, e2c=e2c

  yrange = [-0.10, 0.10]
  xrange = [-0.23, 0.23]
  IF n_elements(e1c) EQ 0 THEN BEGIN 
      IF keyword_set(hirata) THEN BEGIN 
          e1c = struct[w].m_e1_corr_h
          e2c = struct[w].m_e2_corr_h
      ENDIF ELSE BEGIN 
          e1c = struct[w].m_e1_corr/(1. - struct[w].m_r)
          e2c = struct[w].m_e2_corr/(1. - struct[w].m_r)
      ENDELSE 
  ENDIF 

  legend_charsize = 1
  charsize = 1
  tl = 0.04

  erase & multiplot, [2,2], /square

  e1err = sqrt( struct[w].m_e1e1err^2 + 0.32^2 )
  e2err = sqrt( struct[w].m_e2e2err^2 + 0.32^2 )

  plot_fitlin, $
    struct[w].m_e1_psf, e1c, nperEbin, $
    intercept_ge1pe1, intercept_ge1pe1_error, $
    slope_ge1pe1, slope_ge1pe1_error,$
    ytitle="e!D1!N Gal",$
    xrange=xrange, yrange=yrange,ystyle=1, xstyle=1+2, $
    title=title, $
    charsize=charsize, $
    xticklen=tl, yticklen=tl, $
    legend_charsize=legend_charsize, yin_err=e1err
  multiplot

  ebin = 0.015
  cyrange = !y.crange
  plothist, struct[w].m_e1_psf, xhist, yhist, bin=ebin, /noplot
  nyhist = yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  plot_fitlin, $
    struct[w].m_e2_psf, e1c, nperEbin, $
    intercept_ge2pe1, intercept_ge2pe1_error, $
    slope_ge2pe1, slope_ge2pe1_error,$
    xrange=xrange, yrange=yrange,ystyle=1, xstyle=1+2, $
    charsize=charsize, $
    xticklen=tl, yticklen=tl, $
    legend_charsize=legend_charsize, yin_err=e1err
  axis, yaxis=1, ystyle=1, ytitle="e!D1!N Gal", charsize=charsize
  multiplot

  cyrange = !y.crange
  plothist, struct[w].m_e2_psf, xhist, yhist, bin=ebin, /noplot
  nyhist = yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50


  plot_fitlin, $
    struct[w].m_e1_psf, e2c, nperEbin, $
    intercept_ge1pe2, intercept_ge1pe2_error, $
    slope_ge1pe2, slope_ge1pe2_error,$
    xrange=xrange, yrange=yrange,ystyle=1, xstyle=1+2, $
    charsize=charsize, $
    ytitle="e!D2!N Gal",xtitle="e!D1!N PSF",$
    xticklen=tl, yticklen=tl, $
    legend_charsize=legend_charsize, yin_err=e2err

  multiplot

  plot_fitlin, $
    struct[w].m_e2_psf, e2c, nperEbin, $
    intercept_ge2pe2, intercept_ge2pe2_error, $
    slope_ge2pe2, slope_ge2pe2_error,$
    xtitle="e!D2!N PSF",$
    xrange=xrange, yrange=yrange, ystyle=1, xstyle=1+2, $
    charsize=charsize, $
    xticklen=tl, yticklen=tl, $
    legend_charsize=legend_charsize, yin_err=e2err
  axis, yaxis=1, ystyle=1, ytitle="e!D2!N Gal", charsize=charsize
  multiplot,/reset

  fitst = egal_vs_epsf_fitstruct()

  fitst.intercept_ge1pe1=intercept_ge1pe1
  fitst.intercept_ge1pe1_error=intercept_ge1pe1_error
  fitst.slope_ge1pe1=slope_ge1pe1
  fitst.slope_ge1pe1_error=slope_ge1pe1_error
  
  fitst.intercept_ge2pe1=intercept_ge2pe1
  fitst.intercept_ge2pe1_error=intercept_ge2pe1_error
  fitst.slope_ge2pe1=slope_ge2pe1
  fitst.slope_ge2pe1_error=slope_ge2pe1_error
  
  fitst.intercept_ge1pe2=intercept_ge1pe2
  fitst.intercept_ge1pe2_error=intercept_ge1pe2_error
  fitst.slope_ge1pe2=slope_ge1pe2
  fitst.slope_ge1pe2_error=slope_ge1pe2_error
  
  fitst.intercept_ge2pe2=intercept_ge2pe2
  fitst.intercept_ge2pe2_error=intercept_ge2pe2_error
  fitst.slope_ge2pe2=slope_ge2pe2
  fitst.slope_ge2pe2_error=slope_ge2pe2_error
  

END 

FUNCTION egal_vs_epsf_run_photoid_limits, run, rerun=rerun

  IF n_elements(rerun) EQ 0 THEN rerun = sdss_rerun(run)
  minid = ntostr( photoid(run,rerun,1,0,0) )
  maxid = ntostr( photoid(run,rerun,6,5000,5000) )

  clause = 'BETWEEN '+minid+' AND '+maxid

  return,clause
END 

FUNCTION egal_vs_epsf_run_photoid_limits2, run1, run2

  rerun1 = sdss_rerun(run1)
  rerun2 = sdss_rerun(run2)

  minid = ntostr( photoid(run1,rerun1,1,0,0) )
  maxid = ntostr( photoid(run2,rerun2,6,5000,5000) )

  clause = 'BETWEEN '+minid+' AND '+maxid

  return,clause
END 

PRO egal_vs_epsf_getcat, clr, struct, requery=requery

  cstr = '['+ntostr(clr)+']'

  out_dir = '/net/cheops1/data0/esheldon/tmp/'
  outfile = out_dir + 'test_corr_'+!colors[clr]+'.st'

  IF NOT fexist(outfile) OR keyword_set(requery) THEN BEGIN

      print
      print,'Will write to file: ',outfile

      run1 = 2855
      run2 = 2987
      
      w=where(!run_status.tsObj_photo_v GE 5.4 AND $
              !run_status.run GE run1 AND !run_status.run LE run2, nw)
      w2 = where(!run_status.run NE 2964)
      w = w[w2]
      

      FOR i=0L, nw-1 DO BEGIN 

          run = !run_status[w[i]].run

          query  = 'SELECT '+$
            'a.photoid,'+$
            'a.objc_prob_gal, a.cmodel_counts_ext[2], '+$
            'a.seeing'+cstr+', '+$
            $
            'a.m_r'+cstr+', '+$
            'a.m_e1_corr'+cstr+', a.m_e2_corr'+cstr+', '+$
            $
            'a.m_r_h'+cstr+', '+$
            'a.m_e1_corr_h'+cstr+', a.m_e2_corr_h'+cstr+', '+$
            $
            'ae.m_e1e1err'+cstr+', '+$
            'ae.m_e2e2err'+cstr+', '+ $
            $
            'ae.m_e1_psf'+cstr+', ae.m_e2_psf'+cstr+', ae.m_cr4_psf'+cstr+' '+$
            'FROM adatc AS a, adatc_extra as ae '+$
            'WHERE a.photoid '+egal_vs_epsf_run_photoid_limits(run)+' AND '+$
            'a.photoid = ae.photoid'

          print,query
          
          print,'sending query'
          tm=systime(1)
          struct = pgsql_query(query)
          ptime,systime(1)-tm
          
          print,'Writing to file: ',outfile
          write_idlstruct, struct, outfile, /append

          struct = 0
      ENDFOR 

  ENDIF 

  print
  print,'Reading file: ',outfile
  struct = read_idlstruct(outfile)
  return

END 

PRO egal_vs_epsf_process, clr, struct, fitstruct_h, fitstruct, $
  requery=requery

  ;; Try plotting mean galaxy e versus various things

  IF n_elements(struct) EQ 0 OR keyword_set(requery) THEN BEGIN 
      egal_vs_epsf_getcat, clr, struct, requery=requery
  ENDIF 

  plot_dir = '~/plots/corr_egal_vs_epsf/'

  maxmag = 22.0
  minmag = 18.0

  legend_charsize=0.7

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Plot rsmear versus mag
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  wmag = where(struct.cmodel_counts_ext LT maxmag AND $
               struct.cmodel_counts_ext GT minmag, nmag)
  help,wmag,struct

  begplot,name=plot_dir + 'rsmear_'+!colors[clr]+'_vs_rmag.ps'
  !p.multi = [0,0,2]
  ploth, struct[wmag].cmodel_counts_ext, struct[wmag].m_r, /asinh, $
    xrange=[18, 22], yrange = [0,1.2], $
    xtitle = 'cmodel_mags['+!colors[2]+']', ytitle = 'R['+!colors[clr]+']'
  ploth, struct[wmag].cmodel_counts_ext, struct[wmag].m_r_h, /asinh,$
    xrange=[18, 22], yrange = [0,1.2], $
    xtitle = 'cmodel_mags['+!colors[2]+']', ytitle = 'R_h['+!colors[clr]+']'
  !p.multi=0
  key = prompt_kbrd("hit a key")
  endplot

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;
  ;; Hirata r
  ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  begplot,name=plot_dir + 'egal_vs_epsf_h_'+!colors[clr]+'.ps', /color

  wm=where(struct.cmodel_counts_ext LT maxmag AND $
           struct.cmodel_counts_ext GT minmag AND $
           struct.m_r_h GT 0.0 AND struct.m_r_h LT 0.8 AND $
           struct.objc_prob_gal GT 0.8, nobj_used)  



  help,wm,struct
  
  s = sort( struct[wm].m_r_h )
  wm = wm[s]

  ind = lindgen(nobj_used)
  nperRbin = 100000
  nperEbin = 10000

  hist = histogram(ind, bin=nperRbin, reverse_indices=mrev_ind)

  wrh = where(hist EQ nperRbin, nbin)
  help,wrh


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Now plot things in each of the r_h bins
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  fitstruct_h = $
    {bandpass: clr, $
     hirata: 1, $
     nobj_used: nobj_used, $
     fitst: replicate(egal_vs_epsf_fitstruct(), nbin) $
    }
  FOR mi=0L, nbin-1 DO BEGIN 

      ;; Objects in this R bin
      wmi = mrev_ind[ mrev_ind[mi]:mrev_ind[mi+1]-1 ]
      wmi = wm[wmi]

      min_r = min(struct[wmi].m_r_h, max=max_r)
      mean_r = mean(struct[wmi].m_r_h)
      label = ntostr(min_r, 4, /round)+' < R_h < '+ntostr(max_r, 4, /round)

      egal_vs_epsf_plot_fitlin, struct, wmi, nperEbin, label, fitst, /hirata

      fitst.min_r = min_r
      fitst.max_r = max_r
      fitst.mean_r = mean_r
      fitst.label = label

      fitstruct_h.fitst[mi] = fitst

      key = prompt_kbrd("hit a key for the next R bin")
      IF key EQ 'q' THEN return

  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Output the fit structures 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  file = $
    sdssidl_config('shapecorr_dir') + 'corr_egal_vs_epsf/'+$
    'egal_vs_epsf_h_'+!colors[clr]+'.st'
  print
  print,'Writing fit struct: ',file
  write_idlstruct, $
    fitstruct_h.fitst, file, $
    hdrStruct = {bandpass: fitstruct_h.bandpass, $
                 hirata: fitstruct_h.hirata, $
                 nobj_used: fitstruct_h.nobj_used}

  ;; Correct the shapes and re-plot
  nsmooth = 5
  correct_eslope, $
    fitStruct_h.fitst, $
    struct[wm].m_r_h, $
    struct[wm].m_e1_corr_h, $
    struct[wm].m_e2_corr_h, $
    struct[wm].m_e1_psf, $
    struct[wm].m_e2_psf, $
    e1recorr_h, e2recorr_h

  egal_vs_epsf_plot_fitlin, struct, wm, 100000, 'recorr h', $
    e1c=e1recorr_h, e2c=e2recorr_h
  
  endplot  

  ;; Plot the fits
  begplot,name=plot_dir + 'egal_vs_epsf_h_'+!colors[clr]+'_fits.ps', /color
  egal_vs_epsf_plotfits, fitstruct_h, rvals=struct[wm].m_r_h
  endplot






  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;
  ;; Original r
  ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  begplot,name=plot_dir + 'egal_vs_epsf_'+!colors[clr]+'.ps', /color

  wm=where(struct.cmodel_counts_ext LT maxmag AND $
           struct.cmodel_counts_ext GT minmag AND $
           struct.m_r GT 0.0 AND struct.m_r LT 0.8 AND $
           struct.objc_prob_gal GT 0.8, nobj_used)  

  s = sort( struct[wm].m_r )
  wm = wm[s]

  ind = lindgen(nobj_used)

  hist = histogram(ind, bin=nperRbin, reverse_indices=mrev_ind)

  wrh = where(hist EQ nperRbin, nbin)
  help,wrh


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Now plot things in each of the r bins
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  fitstruct = $
    {bandpass: clr, $
     hirata: 0, $
     nobj_used: nobj_used, $
     fitst: replicate(egal_vs_epsf_fitstruct(), nbin) $
    }

  FOR mi=0L, nbin-1 DO BEGIN 

      ;; Objects in this R bin
      wmi = mrev_ind[ mrev_ind[mi]:mrev_ind[mi+1]-1 ]
      wmi = wm[wmi]

      min_r = min(struct[wmi].m_r, max=max_r)
      mean_r = mean(struct[wmi].m_r)
      label = ntostr(min_r, 4, /round)+' < R < '+ntostr(max_r, 4, /round)

      egal_vs_epsf_plot_fitlin, struct, wmi, nperEbin, 'orig '+label, fitst

      fitst.min_r = min_r
      fitst.max_r = max_r
      fitst.mean_r = mean_r
      fitst.label = label

      fitstruct.fitst[mi] = fitst

      key = prompt_kbrd("hit a key for the next R bin")
      IF key EQ 'q' THEN return

  ENDFOR 

  file = $
    sdssidl_config('shapecorr_dir') + 'corr_egal_vs_epsf/'+$
    'egal_vs_epsf_'+!colors[clr]+'.st'
  print,'Writing fit struct: ',file
  write_idlstruct, $
    fitstruct.fitst, file, $
    hdrStruct = {bandpass: fitstruct.bandpass, $
                 hirata: fitstruct.hirata, $
                 nobj_used: fitstruct.nobj_used}

  ;; Correct the shapes and re-plot
  corr = 1.0/(1.0 - struct[wm].m_r)
  nsmooth = 5
  correct_eslope, $
    fitstruct.fitst, $
    struct[wm].m_r, $
    struct[wm].m_e1_corr*corr, $
    struct[wm].m_e2_corr*corr, $
    struct[wm].m_e1_psf, $
    struct[wm].m_e2_psf, $
    e1recorr, e2recorr, /original

  egal_vs_epsf_plot_fitlin, struct, wm, 100000, 'recorr h',$
    e1c=e1recorr, e2c=e2recorr

  endplot

  begplot,name=plot_dir + 'egal_vs_epsf_'+!colors[clr]+'_fits.ps', /color
  egal_vs_epsf_plotfits, fitstruct, rvals=struct[wm].m_r
  endplot




  return

END 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Main program
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO egal_vs_epsf, clr, struct, fitstruct_h, fitstruct, requery=requery

  IF n_elements(clr) EQ 0 THEN BEGIN 
      print,'-Syntax: egal_vs_epsf, clr [, struct, fitstruct_h, fitstruct, /requery]'
      return
  ENDIF 

  IF keyword_set(princeton) THEN BEGIN 

  ENDIF ELSE BEGIN 

      egal_vs_epsf_process, clr, struct, fitstruct_h, fitstruct, $
        requery=requery

  ENDELSE 

END 
