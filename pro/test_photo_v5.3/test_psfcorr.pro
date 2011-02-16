PRO test_psfcorr, struct, run, camcol, clr, gind, e1, e2, diff_e1, diff_e2, cutzero=cutzero

  ;; only need the "new" to see how many total are thrown out
  ;; via various flags

  IF n_elements(struct) EQ 0 THEN BEGIN 

      tags = ['run','rerun','camcol','field','id', 'ra','dec',$
              'ixx','iyy','ixy','rho4','momerr',$
              'psfixx','psfiyy','psfixy', 'psfrho4',$
              'm_rr_cc','m_rr_ccerr','m_cr4',$
              'm_rr_cc_psf','m_cr4_psf',$
              'm_e1','m_e2','m_e1e1err', 'm_e1e2err', 'm_e2e2err',$
              'm_e1_psf', 'm_e2_psf',$
              'e1','e2','r','petrocounts','counts_model','flags','flags2','objc_flags']

      rerun=9
      fetch_dir, run, camcol, rerun, dir, corrdir=corrdir
      read_tsobj, corrdir, struct, /all, taglist=tags, front='adatc'

  ENDIF
 
  IF n_elements(gind) EQ 0 THEN BEGIN 

      ;; get bright ones that pass flag cuts
      nstruct = n_elements(struct)
      print,'Making cuts on struct'
      wbright=where(struct.petrocounts[2] LE 22.,nbright)
      ;;help,new,wbright
      print,'Removed '+ntostr(nstruct - nbright)+'/'+ntostr(nstruct)+' = '+$
            ntostr((nstruct - nbright)/float(nstruct))+' faint objects'
      
      make_flag_struct,fs
      fs.amoment_faint = 'N'
      flag_select, struct[wbright], fs, clr, si_faint
      nkeep = n_elements(si_faint)
      print,'Removed '+ntostr(nbright - nkeep)+'/'+ntostr(nbright)+' = '+$
            ntostr((nbright-nkeep)/float(nbright))+' amoment_faint'
      
      make_flag_struct,fs
      fs.amoment_shift = 'N'
      flag_select, struct[wbright], fs, clr, si_shift
      nkeep = n_elements(si_shift)
      print,'Removed '+ntostr(nbright - nkeep)+'/'+ntostr(nbright)+' = '+$
            ntostr((nbright-nkeep)/float(nbright))+' amoment_shift'

      delvarx,fs
      make_flag_struct,fs
      fs.amoment_maxiter = 'N'
      flag_select, struct[wbright], fs, clr, si_maxiter
      ;;help,new[wbright],si_maxiter
      nkeep = n_elements(si_maxiter)
      print,'Removed '+ntostr(nbright - nkeep)+'/'+ntostr(nbright)+' = '+$
            ntostr((nbright-nkeep)/float(nbright))+' amoment_maxiter'
      
      
      fs.amoment_faint = 'N'
      fs.amoment_shift = 'N'
      fs.amoment_maxiter = 'N'
      flag_select, struct[wbright], fs, clr, si_amomentcheck
      ;;help,new[wbright],si_amomentcheck
      nkeep = n_elements(si_amomentcheck)
      print,'Removed '+ntostr(nbright - nkeep)+'/'+ntostr(nbright)+' = '+$
            ntostr((nbright-nkeep)/float(nbright))+' total bad amoment'

      gind = wbright[si_amomentcheck]

      gind2 = where(struct[gind].e1[clr] NE 1.e10, nkeep2)
      gind = gind[gind2]
      
      print,'Removed '+ntostr(nkeep-nkeep2)+'/'+ntostr(nkeep)+' = '+$
            ntostr((nkeep-nkeep2)/float(nkeep))+' total bad amoment'
      nkeep = nkeep2

      ;; get good galaxies
      corr = struct[gind].m_rr_cc_psf[clr]/struct[gind].m_rr_cc[clr]*(4./struct[gind].m_cr4_psf[clr] -1.)/(4./struct[gind].m_cr4[clr]-1.)
      gind_2 = where(corr LT 0.8 AND corr GT 0.0 AND $
                     struct[gind].r[clr] LT 0.8 AND struct[gind].r[clr] GT 0.0 AND $
                     struct[gind].m_e1e1err[clr] LT 0.6 AND $
                     struct[gind].m_e2e2err[clr] LT 0.6 AND $
                     struct[gind].m_e1e1err[clr] GT 0.0 AND $
                     struct[gind].m_e2e2err[clr] GT 0.0, ngal)
      gind = gind[gind_2]
      print,'Removed '+ntostr(nkeep-ngal)+'/'+ntostr(nkeep)+' = '+$
            ntostr((nkeep-ngal)/float(nkeep))+' "stars" and bad measurements'
      print,'Total Gals: ',ngal
      print

  ENDIF 
  ;; plot various things new vs old as a function of magnitude


  ;; correct for PSF anisotropy
  corr = struct[gind].m_rr_cc_psf[clr]/struct[gind].m_rr_cc[clr]*(4./struct[gind].m_cr4_psf[clr] -1.)/(4./struct[gind].m_cr4[clr]-1.)

  e1 = struct[gind].m_e1[clr]
  e2 = struct[gind].m_e2[clr]
  e1_psf = struct[gind].m_e1_psf[clr]
  e2_psf = struct[gind].m_e2_psf[clr]
  e1corr = struct[gind].m_e1[clr] - corr*struct[gind].m_e1_psf[clr]
  e2corr = struct[gind].m_e2[clr] - corr*struct[gind].m_e2_psf[clr]

  frac = 1./15.
  snperbin = long( round(frac*n_elements(e1_psf)) )

  ;; correct for PSF anisotropy
  size = struct[gind].ixx[clr] + struct[gind].iyy[clr]
  oe1 = (struct[gind].ixx[clr] - struct[gind].iyy[clr])/size
  oe2 = 2.*struct[gind].ixy[clr]/size
  size_psf = struct[gind].psfixx[clr] + struct[gind].psfiyy[clr]
  oe1_psf = (struct[gind].psfixx[clr] - struct[gind].psfiyy[clr])/size_psf
  oe2_psf = 2.*struct[gind].psfixy[clr]/size_psf
  oe1corr = struct[gind].e1[clr]
  oe2corr = struct[gind].e2[clr]

  ;;;;;;;;;;;;;
  ;; now plots
  ;;;;;;;;;;;;;

  !p.multi=[0,0,2]
  ;; need to compile corrtest_plots to get this func
  corrtest_plots_getrange, oe1_psf, oe1, oe2_psf, oe2, snperbin,$
                           xrange, yrange

  plot_fitlin, e1_psf, e1, snperbin, $
               xtitle="e!D1!N Star",ytitle='e!D1!N Gal', $
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1,$
               title=!colors[clr]+"-band  Uncorrected  PHOTO", /iso

  plot_fitlin, oe1_psf, oe1, snperbin, $
               xtitle="e!D1!N Star",ytitle='e!D1!N Gal', $
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1,$
               title=!colors[clr]+"-band  Uncorrected  Michigan", /iso

  IF !d.name EQ 'X' THEN key=get_kbrd(1)

  plot_fitlin, e2_psf, e2, snperbin, $
               xtitle="e!D2!N Star",ytitle='e!D2!N Gal', $
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1,$
               title=!colors[clr]+"-band  Uncorrected  PHOTO", /iso

  plot_fitlin, oe2_psf, oe2, snperbin, $
               xtitle="e!D2!N Star",ytitle='e!D2!N Gal', $
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1,$
               title=!colors[clr]+"-band  Uncorrected  Michigan", /iso

  ;; plot corrected on same scale first

  IF !d.name EQ 'X' THEN key=get_kbrd(1)

  corrtest_plots_getrange, oe1_psf, oe1corr, oe2_psf, oe2corr, snperbin,$
                           cxrange, cyrange

  plot_fitlin, e1_psf, e1corr, snperbin, $
               xtitle="e!D1!N Star",ytitle='e!D1!N Gal', $
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1,$
               title=!colors[clr]+"-band  Corrected  PHOTO", /iso

  plot_fitlin, oe1_psf, oe1corr, snperbin, $
               xtitle="e!D1!N Star",ytitle='e!D1!N Gal', $
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1,$
               title=!colors[clr]+"-band  Corrected  Michigan", /iso

  IF !d.name EQ 'X' THEN key=get_kbrd(1)

  plot_fitlin, e2_psf, e2corr, snperbin, $
               xtitle="e!D2!N Star",ytitle='e!D2!N Gal', $
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1,$
               title=!colors[clr]+"-band  Ccorrected  PHOTO", /iso

  plot_fitlin, oe2_psf, oe2corr, snperbin, $
               xtitle="e!D2!N Star",ytitle='e!D2!N Gal', $
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1,$
               title=!colors[clr]+"-band  Ccorrected  Michigan", /iso

  ;; Now on smaller range

  IF !d.name EQ 'X' THEN key=get_kbrd(1)

  corrtest_plots_getrange, oe1_psf, oe1corr, oe2_psf, oe2corr, snperbin,$
                           cxrange, cyrange

  plot_fitlin, e1_psf, e1corr, snperbin, $
               xtitle="e!D1!N Star",ytitle='e!D1!N Gal', $
               xrange=cxrange,yrange=cyrange,ystyle=1,xstyle=1,$
               title=!colors[clr]+"-band  Corrected  PHOTO"

  plot_fitlin, oe1_psf, oe1corr, snperbin, $
               xtitle="e!D1!N Star",ytitle='e!D1!N Gal', $
               xrange=cxrange,yrange=cyrange,ystyle=1,xstyle=1,$
               title=!colors[clr]+"-band  Corrected  Michigan"

  IF !d.name EQ 'X' THEN key=get_kbrd(1)

  plot_fitlin, e2_psf, e2corr, snperbin, $
               xtitle="e!D2!N Star",ytitle='e!D2!N Gal', $
               xrange=cxrange,yrange=cyrange,ystyle=1,xstyle=1,$
               title=!colors[clr]+"-band  Ccorrected  PHOTO"

  plot_fitlin, oe2_psf, oe2corr, snperbin, $
               xtitle="e!D2!N Star",ytitle='e!D2!N Gal', $
               xrange=cxrange,yrange=cyrange,ystyle=1,xstyle=1,$
               title=!colors[clr]+"-band  Ccorrected  Michigan"



return










END 
