PRO test_error, prun1, prun2, camcol, clr, mprun1, mprun2, wprun1, wprun2, e1, e2, diff_e1, diff_e2, cutzero=cutzero

  run1 = 756
  rerun1 = 9
  run2 = 745
  rerun2 = 9
  IF n_elements(prun1) EQ 0 THEN BEGIN 

      tags = ['run','rerun','camcol','field','id', 'ra','dec',$
              'm_rr_cc','m_rr_ccerr','m_cr4',$
              'm_rr_cc_psf','m_cr4_psf',$
              'm_e1','m_e2','m_e1e1err', 'm_e1e2err', 'm_e2e2err',$
              'm_e1_psf', 'm_e2_psf',$
              'petrocounts','counts_model','flags','flags2','objc_flags']

      fetch_dir, run1, camcol, rerun1, dir, corrdir=corrdir
      read_tsobj, corrdir, prun1, /all, taglist=tags, front='adatc'

  ENDIF 

  IF n_elements(prun2) EQ 0 THEN BEGIN

      tags = ['id', 'ra','dec',$
              'm_rr_cc','m_rr_ccerr','m_cr4',$
              'm_rr_cc_psf','m_cr4_psf',$
              'm_e1','m_e2','m_e1e1err', 'm_e1e2err', 'm_e2e2err',$
              'm_e1_psf', 'm_e2_psf',$
              'petrocounts','counts_model','flags','flags2','objc_flags']


      fetch_dir, run2, camcol, rerun2, dir, corrdir=corrdir
      read_tsobj, corrdir, prun2, /all, taglist=tags, front='adatc'

  ENDIF 

  IF n_elements(mprun1) EQ 0 OR n_elements(mprun2) EQ 0 THEN BEGIN 

      ;; get bright ones that pass flag cuts
      ;; Prun2
      nprun2 = n_elements(prun2)
      print,'Making cuts on prun2 struct'
      wbright=where(prun2.petrocounts[2] LE 22,nbright)
      ;;help,prun2,wbright
      print,'Removed '+ntostr(nprun2 - nbright)+'/'+ntostr(nprun2)+' = '+$
            ntostr((nprun2 - nbright)/float(nprun2))+' faint objects'
      
      make_flag_struct,fs
      fs.amoment_faint = 'N'
      flag_select, prun2[wbright], fs, clr, si_faint
      nkeep = n_elements(si_faint)
      print,'Removed '+ntostr(nbright - nkeep)+'/'+ntostr(nbright)+' = '+$
            ntostr((nbright-nkeep)/float(nbright))+' amoment_faint'
      
      make_flag_struct,fs
      fs.amoment_shift = 'N'
      flag_select, prun2[wbright], fs, clr, si_shift
      nkeep = n_elements(si_shift)
      print,'Removed '+ntostr(nbright - nkeep)+'/'+ntostr(nbright)+' = '+$
            ntostr((nbright-nkeep)/float(nbright))+' amoment_shift'

      delvarx,fs
      make_flag_struct,fs
      fs.amoment_maxiter = 'N'
      flag_select, prun2[wbright], fs, clr, si_maxiter
      ;;help,prun2[wbright],si_maxiter
      nkeep = n_elements(si_maxiter)
      print,'Removed '+ntostr(nbright - nkeep)+'/'+ntostr(nbright)+' = '+$
            ntostr((nbright-nkeep)/float(nbright))+' amoment_maxiter'
      
      
      fs.amoment_faint = 'N'
      fs.amoment_shift = 'N'
      fs.amoment_maxiter = 'N'
      flag_select, prun2[wbright], fs, clr, si_amomentcheck
      ;;help,prun2[wbright],si_amomentcheck
      nkeep = n_elements(si_amomentcheck)
      print,'Removed '+ntostr(nbright - nkeep)+'/'+ntostr(nbright)+' = '+$
            ntostr((nbright-nkeep)/float(nbright))+' total bad amoment'

      gprun2 = wbright[si_amomentcheck]
      gprun2_2 = where(prun2[gprun2].m_e1_psf[clr] NE 0.0, npsf)
      gprun2 = gprun2[gprun2_2]
      print,'Removed '+ntostr(nkeep - npsf)+'/'+ntostr(nkeep)+' = '+$
            ntostr((nkeep-npsf)/float(nkeep))+' total bad PSF amoment'

      ;; get galaxies
      corr2 = prun2[gprun2].m_rr_cc_psf[clr]/prun2[gprun2].m_rr_cc[clr]*(4./prun2[gprun2].m_cr4_psf[clr] -1.)/(4./prun2[gprun2].m_cr4[clr]-1.)
      gprun2_2 = where(corr2 LT 0.8 AND corr2 GT 0.0, ngal)
      gprun2 = gprun2[gprun2_2]
      print,'Removed '+ntostr(npsf - ngal)+'/'+ntostr(npsf)+' = '+$
            ntostr((npsf-ngal)/float(npsf))+' stars'
      print,'Total Gals: ',ngal
      print

      ;; get bright ones that pass flag cuts
      ;; Prun1
      nprun1 = n_elements(prun1)
      print,'Making cuts on prun1 struct'
      wbright=where(prun1.petrocounts[2] LE 22,nbright)
      ;;help,prun1,wbright
      print,'Removed '+ntostr(nprun1 - nbright)+'/'+ntostr(nprun1)+' = '+$
            ntostr((nprun1 - nbright)/float(nprun1))+' faint objects'
      
      make_flag_struct,fs
      fs.amoment_faint = 'N'
      flag_select, prun1[wbright], fs, clr, si_faint
      nkeep = n_elements(si_faint)
      print,'Removed '+ntostr(nbright - nkeep)+'/'+ntostr(nbright)+' = '+$
            ntostr((nbright-nkeep)/float(nbright))+' amoment_faint'
      
      make_flag_struct,fs
      fs.amoment_shift = 'N'
      flag_select, prun1[wbright], fs, clr, si_shift
      nkeep = n_elements(si_shift)
      print,'Removed '+ntostr(nbright - nkeep)+'/'+ntostr(nbright)+' = '+$
            ntostr((nbright-nkeep)/float(nbright))+' amoment_shift'

      delvarx,fs
      make_flag_struct,fs
      fs.amoment_maxiter = 'N'
      flag_select, prun1[wbright], fs, clr, si_maxiter
      ;;help,prun1[wbright],si_maxiter
      nkeep = n_elements(si_maxiter)
      print,'Removed '+ntostr(nbright - nkeep)+'/'+ntostr(nbright)+' = '+$
            ntostr((nbright-nkeep)/float(nbright))+' amoment_maxiter'
      
      
      fs.amoment_faint = 'N'
      fs.amoment_shift = 'N'
      fs.amoment_maxiter = 'N'
      flag_select, prun1[wbright], fs, clr, si_amomentcheck
      ;;help,prun1[wbright],si_amomentcheck
      nkeep = n_elements(si_amomentcheck)
      print,'Removed '+ntostr(nbright - nkeep)+'/'+ntostr(nbright)+' = '+$
            ntostr((nbright-nkeep)/float(nbright))+' total bad amoment'

      gprun1 = wbright[si_amomentcheck]

      gprun1 = wbright[si_amomentcheck]
      gprun1_2 = where(prun1[gprun1].m_e1_psf[clr] NE 0.0, npsf)
      gprun1 = gprun1[gprun1_2]
      print,'Removed '+ntostr(nkeep - npsf)+'/'+ntostr(nkeep)+' = '+$
            ntostr((nkeep-npsf)/float(nkeep))+' total bad PSF amoment'

      ;; get galaxies
      corr1 = prun1[gprun1].m_rr_cc_psf[clr]/prun1[gprun1].m_rr_cc[clr]*(4./prun1[gprun1].m_cr4_psf[clr] -1.)/(4./prun1[gprun1].m_cr4[clr]-1.)
      gprun1_2 = where(corr1 LT 0.8 AND corr1 GT 0.0, ngal)
      gprun1 = gprun1[gprun1_2]
      print,'Removed '+ntostr(npsf - ngal)+'/'+ntostr(npsf)+' = '+$
            ntostr((npsf-ngal)/float(npsf))+' stars'
      print,'Total Gals: ',ngal
      print

      print
      print,'Matching by (ra,dec)'
      rad = 1d/3600d
      allow = 1
      close_match_radec, prun1[gprun1].ra, prun1[gprun1].dec, $
                         prun2[gprun2].ra, prun2[gprun2].dec, mprun1,mprun2,$
                         rad, allow
      
      mprun1 = gprun1[mprun1]
      mprun2 = gprun2[mprun2]
      print
  ENDIF 
  ;; plot various things prun2 vs prun1 as a function of magnitude

  ;; bin by 0.5 in mag
  binsize = 1
  rmin = 15.0
  rmax = 22.0
  nbin = 7

  magptr = ptrarr(nbin)

  hist = histogram(prun2[mprun2].petrocounts[2], min=rmin, max=rmax,$
                       binsize = binsize, reverse=rev_ind)

  ;; 
  FOR binnum=0L,nbin-1 DO BEGIN 

      IF rev_ind[binnum] NE rev_ind[binnum+1] THEN BEGIN 

          w = rev_ind[ rev_ind[binnum]:rev_ind[binnum+1]-1 ]
          magptr[binnum] = ptr_new(w, /no_copy)

      ENDIF 

  ENDFOR 

  xrange=[-1,1]
  yrange=[-1,1]
  hxrange = [-0.5,0.5]
  dbinsize = 0.01

  ;; e1 vs. e1
  print
  print,'e1 vs e1'
  !p.multi = [0,2,4]
  xtitle = 'e!D1!N run'+ntostr(run1)
  ytitle = 'e!D1!N run'+ntostr(run2)
  
  FOR binnum=0L, nbin-1 DO BEGIN 

      mmin = rmin + binnum*binsize
      mmax = rmin + (binnum+1)*binsize

      title = 'r petro ['+ntostr(mmin,4)+', '+ntostr(mmax,4)+']'

      wprun1 = mprun1[*magptr[binnum]]
      wprun2 = mprun2[*magptr[binnum]]

      e1_run2 = prun2[wprun2].m_e1[clr]
      e1_run1 = prun1[wprun1].m_e1[clr]

      corr1 = prun1[wprun1].m_rr_cc_psf[clr]/prun1[wprun1].m_rr_cc[clr]*(4./prun1[wprun1].m_cr4_psf[clr] -1.)/(4./prun1[wprun1].m_cr4[clr]-1.)
      corr2 = prun2[wprun2].m_rr_cc_psf[clr]/prun2[wprun2].m_rr_cc[clr]*(4./prun2[wprun2].m_cr4_psf[clr] -1.)/(4./prun2[wprun2].m_cr4[clr]-1.)

      e1_run1 = (e1_run1 - corr1*prun1[wprun1].m_e1_psf[clr])/(1. - corr1)
      e1_run2 = (e1_run2 - corr2*prun2[wprun2].m_e1_psf[clr])/(1. - corr2)

      diff_e1 = e1_run2 - e1_run1

      print,title
      print,'  Median diff_e1 = ',median(diff_e1)
      ;;print,'Mean diff_e1 = ',mean(diff_e1)
      print,'  Sdev diff_e1 = ',sdev(diff_e1)
      print,'  Skewness diff_e1 = ',skewness(diff_e1)

      ;; chi-squared
      ;;print,min(prun1[wprun1].m_e1e1err[clr]),min(prun2[wprun2].m_e1e1err[clr])
      error2 = (prun1[wprun1].m_e1e1err[clr]/(1.-corr1)) + (prun2[wprun2].m_e1e1err[clr]/(1.-corr2))
      chi_squared_per = total( diff_e1^2/error2 )/n_elements(wprun1)
      print,'---Chi-squared: ',chi_squared_per

      error2 = (prun1[wprun1].m_e1e1err[clr]/(1.-corr1))^2 + (prun2[wprun2].m_e1e1err[clr]/(1.-corr2))^2
      chi_squared_per = total( diff_e1^2/error2 )/n_elements(wprun1)
      print,'---Chi-squared: ',chi_squared_per
      print

      aplot, 1, e1_run1, e1_run2, psym=3,$
            xtitle=xtitle, ytitle=ytitle, title=title, $
            xrange=xrange, yrange=yrange
      oplot, [-1000, 1000],[-1000,1000], color=!blue

      plothist, diff_e1, bin=dbinsize, xrange=hxrange,$
                xtitle = ytitle + ' - '+ xtitle,title=!csym.chi+'!U2!N/'+!csym.nu+' = '+ntostr(chi_squared_per)


  ENDFOR 


  ;; e2 vs. e2
  print
  print,'e2 vs e2'
  !p.multi = [0,2,4]
 
  xtitle = 'e!D2!N run'+ntostr(run1)
  ytitle = 'e!D2!N run'+ntostr(run2)
  
  FOR binnum=0L, nbin-1 DO BEGIN 

      mmin = rmin + binnum*binsize
      mmax = rmin + (binnum+1)*binsize

      title = 'r petro ['+ntostr(mmin,4)+', '+ntostr(mmax,4)+']'

      wprun1 = mprun1[*magptr[binnum]]
      wprun2 = mprun2[*magptr[binnum]]

      e2_run2 = prun2[wprun2].m_e2[clr]
      e2_run1 = prun1[wprun1].m_e2[clr]

      corr1 = prun1[wprun1].m_rr_cc_psf[clr]/prun1[wprun1].m_rr_cc[clr]*(4./prun1[wprun1].m_cr4_psf[clr] -1.)/(4./prun1[wprun1].m_cr4[clr]-1.)
      corr2 = prun2[wprun2].m_rr_cc_psf[clr]/prun2[wprun2].m_rr_cc[clr]*(4./prun2[wprun2].m_cr4_psf[clr] -1.)/(4./prun2[wprun2].m_cr4[clr]-1.)

      e2_run1 = (e2_run1 - corr1*prun1[wprun1].m_e2_psf[clr])/(1. - corr1)
      e2_run2 = (e2_run2 - corr2*prun2[wprun2].m_e2_psf[clr])/(1. - corr2)

      diff_e2 = e2_run2 - e2_run1

      print,title
      print,'  Median diff_e2 = ',median(diff_e2)
      ;;print,'Mean diff_e2 = ',mean(diff_e2)
      print,'  Sdev diff_e2 = ',sdev(diff_e2)
      print,'  Skewness diff_e2 = ',skewness(diff_e2)

      ;; chi-squared
      error2 = (prun1[wprun1].m_e2e2err[clr]/(1.-corr1)) + (prun2[wprun2].m_e2e2err[clr]/(1.-corr2))
      chi_squared_per = total( diff_e2^2/error2 )/n_elements(wprun1)
      print,'---Chi-squared: ',chi_squared_per

      error2 = (prun1[wprun1].m_e2e2err[clr]/(1.-corr1))^2 + (prun2[wprun2].m_e2e2err[clr]/(1.-corr2))^2
      chi_squared_per = total( diff_e2^2/error2 )/n_elements(wprun1)
      print,'---Chi-squared: ',chi_squared_per
      print

      aplot, 1, e2_run1, e2_run2, psym=3,$
            xtitle=xtitle, ytitle=ytitle, title=title, $
            xrange=xrange, yrange=yrange
      oplot, [-1000, 1000],[-1000,1000], color=!blue

      plothist, diff_e2, bin=dbinsize, xrange=hxrange,$
                xtitle = ytitle + ' - '+ xtitle,title=!csym.chi+'!U2!N/'+!csym.nu+' = '+ntostr(chi_squared_per)


  ENDFOR 




  ptrarr_free, magptr

END 
