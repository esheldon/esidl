PRO test_vsold_samefile, struct, run, camcol, clr, gind, e1, e2, diff_e1, diff_e2, cutzero=cutzero, gal=gal

  ;; only need the "new" to see how many total are thrown out
  ;; via various flags

  IF n_elements(struct) EQ 0 THEN BEGIN 

      tags = ['id', 'ra','dec',$
              'ixx','iyy','ixy','rho4','momerr',$
              'psfixx','psfiyy','psfixy', 'psfrho4',$
              'm_rr_cc','m_rr_ccerr','m_cr4',$
              'm_rr_cc_psf','m_cr4_psf',$
              'm_e1','m_e2','m_e1e1err', 'm_e1e2err', 'm_e2e2err',$
              'm_e1_psf', 'm_e2_psf',$
              'e1','e2','r', 'petrocounts','counts_model','flags','flags2','objc_flags']

      rerun=9
      fetch_dir, run, camcol, rerun, dir, corrdir=corrdir
      read_tsobj, corrdir, struct, /all, taglist=tags, front='adatc'

  ENDIF
 
  IF n_elements(gind) EQ 0 THEN BEGIN 

      ;; get bright ones that pass flag cuts
      nstruct = n_elements(struct)
      print,'Making cuts on struct'
      wbright=where(struct.petrocounts[2] LE 22,nbright)
      ;;help,new,wbright
      print,'Removed '+ntostr(nstruct - nbright)+'/'+ntostr(nstruct)+' = '+$
            ntostr((nstruct - nbright)/float(nstruct))+' faint objects'
      
;      make_flag_struct,fs
;      fs.amoment_faint = 'N'
;      flag_select, struct[wbright], fs, clr, si_faint
;      nkeep = n_elements(si_faint)
;      print,'Removed '+ntostr(nbright - nkeep)+'/'+ntostr(nbright)+' = '+$
;            ntostr((nbright-nkeep)/float(nbright))+' amoment_faint'
      
;      make_flag_struct,fs
;      fs.amoment_shift = 'N'
;      flag_select, struct[wbright], fs, clr, si_shift
;      nkeep = n_elements(si_shift)
;      print,'Removed '+ntostr(nbright - nkeep)+'/'+ntostr(nbright)+' = '+$
;            ntostr((nbright-nkeep)/float(nbright))+' amoment_shift'

;      delvarx,fs
;      make_flag_struct,fs
;      fs.amoment_maxiter = 'N'
;      flag_select, struct[wbright], fs, clr, si_maxiter
;      nkeep = n_elements(si_maxiter)
;      print,'Removed '+ntostr(nbright - nkeep)+'/'+ntostr(nbright)+' = '+$
;            ntostr((nbright-nkeep)/float(nbright))+' amoment_maxiter'
      
      make_flag_struct, fs
      fs.amoment_faint = 'N'
      fs.amoment_shift = 'N'
      fs.amoment_maxiter = 'N'
      flag_select, struct[wbright], fs, clr, si_amomentcheck
      ;;help,new[wbright],si_amomentcheck
      nkeep = n_elements(si_amomentcheck)
      print,'Removed '+ntostr(nbright - nkeep)+'/'+ntostr(nbright)+' = '+$
            ntostr((nbright-nkeep)/float(nbright))+' total bad amoment'

      gind = wbright[si_amomentcheck]

      gind2 = where(struct[gind].e1[clr] NE 1.e10)
      gind = gind[gind2]
      nkeep = n_elements(gind)

      w=where(struct[gind].m_e1[clr] NE 0.0 AND struct[gind].m_e2[clr] NE 0.0,ngood) 
      gind = gind[w]
      ;;IF nw NE 0 THEN message,'STILL ZEROS!!!  '+ntostr(nw)+'/'+ntostr(nkeep)
      
      print,'Removed '+ntostr(nkeep - ngood)+'/'+ntostr(nkeep)+' = '+$
            ntostr((nkeep-ngood)/float(nkeep))+' STILL ZERO!!'

      ;; get good galaxies
      IF keyword_set(gal) THEN BEGIN 
          corr = struct[gind].m_rr_cc_psf[clr]/struct[gind].m_rr_cc[clr]*(4./struct[gind].m_cr4_psf[clr] -1.)/(4./struct[gind].m_cr4[clr]-1.)
          gind_2 = where(corr LT 0.8 AND corr GT 0.0 AND $
                         struct[gind].r[clr] LT 0.8 AND struct[gind].r[clr] GT 0.0 AND $
                         struct[gind].m_e1e1err[clr] LT 0.6 AND $
                         struct[gind].m_e2e2err[clr] LT 0.6 AND $
                         struct[gind].m_e1e1err[clr] GT 0.0 AND $
                         struct[gind].m_e2e2err[clr] GT 0.0, ngal)
          gind = gind[gind_2]
          print,'Removed '+ntostr(ngood-ngal)+'/'+ntostr(ngood)+' = '+$
                ntostr((ngood-ngal)/float(ngood))+' "stars" and bad measurements'
          print,'Total Gals: ',ngal
          print
      ENDIF 
  ENDIF 
  ;; plot various things new vs old as a function of magnitude

  ;; bin by 0.5 in mag
  binsize = 1
  rmin = 15.0
  rmax = 22.0
  nbin = 7

  magptr = ptrarr(nbin)

  hist = histogram(struct[gind].petrocounts[2], min=rmin, max=rmax,$
                       binsize = binsize, reverse=rev_ind)

  ;; 
  FOR binnum=0L,nbin-1 DO BEGIN 

      IF rev_ind[binnum] NE rev_ind[binnum+1] THEN BEGIN 

          w = rev_ind[ rev_ind[binnum]:rev_ind[binnum+1]-1 ]
          magptr[binnum] = ptr_new(w, /no_copy)

      ENDIF 

  ENDFOR 

  dbinsize = 0.001

  ;; e1 vs. e1
  print
  print,'e1 vs e1'
  !p.multi = [0,2,4]
  xrange=[-1,1]
  yrange=[-1,1]
  xtitle = 'e!D1!N Michigan'
  ytitle = 'e!D1!N PHOTO v5.3'
  
  FOR binnum=0L, nbin-1 DO BEGIN 

      mmin = rmin + binnum*binsize
      mmax = rmin + (binnum+1)*binsize

      title = 'r petro ['+ntostr(mmin,4)+', '+ntostr(mmax,4)+']'

      wuse = gind[*magptr[binnum]]

      e1old = (struct[wuse].ixx[clr] - struct[wuse].iyy[clr])/(struct[wuse].ixx[clr] + struct[wuse].iyy[clr])
      ;;e1old = struct[wuse].e1[clr]
      e1new = struct[wuse].m_e1[clr]
      ;;corr = struct[wuse].m_rr_cc_psf[clr]/struct[wuse].m_rr_cc[clr]*(4./struct[wuse].m_cr4_psf[clr] -1.)/(4./struct[wuse].m_cr4[clr]-1.)
      ;;e1new = e1new - corr*struct[wuse].m_e1_psf[clr]
      aplot, 1, e1old, e1new, psym=3,$
            xtitle=xtitle, ytitle=ytitle, title=title, $
            xrange=xrange
      oplot, [-1000, 1000],[-1000,1000], color=!blue

      diff_e1 = e1new - e1old
      plothist, diff_e1, bin=dbinsize, xrange=[-.05,.05],$
                xtitle = 'e!D1!N PHOTO - e!D1!N Michigan'

      print,title
      print,'  Median diff_e1 = ',median(diff_e1)
      ;;print,'Mean diff_e1 = ',mean(diff_e1)
      print,'  Sdev diff_e1 = ',sdev(diff_e1)
      print,'  Skewness diff_e1 = ',skewness(diff_e1)

  ENDFOR 


  ;; e2 vs. e2
  !p.multi = [0,2,4]
  xrange=[-1,1]
  yrange=[-1,1]
  xtitle = 'e!D2!N Michigan'
  ytitle = 'e!D2!N PHOTO v5.3'
  print
  print,'e2 vs e2'

  FOR binnum=0L, nbin-1 DO BEGIN 

      mmin = rmin + binnum*binsize
      mmax = rmin + (binnum+1)*binsize

      title = 'r petro ['+ntostr(mmin,4)+', '+ntostr(mmax,4)+']'

      wuse = gind[*magptr[binnum]]

      e2 = 2.*struct[wuse].ixy[clr]/(struct[wuse].ixx[clr] + struct[wuse].iyy[clr])
      aplot, 1, e2, struct[wuse].m_e2[clr], psym=3,$
            xtitle=xtitle, ytitle=ytitle, title=title, $
            xrange=xrange
      oplot, [-1000, 1000],[-1000,1000], color=!blue

      diff_e2 = struct[wuse].m_e2[clr] - e2
      plothist, diff_e2, bin=dbinsize, xrange=[-.05,.05],$
                xtitle = 'e!D2!N PHOTO - e!D2!N Michigan'

      print,title
      print,'  Median diff_e2 = ',median(diff_e2)
      ;;print,'Mean diff_e2 = ',mean(diff_e2)
      print,'  Sdev diff_e2 = ',sdev(diff_e2)
      print,'  Skewness diff_e2 = ',skewness(diff_e2)
           
  ENDFOR 

  ptrarr_free, magptr

  !p.multi=[0,2,4]
  ;; now PSF info
  xtitle = 'e!D1!N Michigan'
  ytitle = 'e!D1!N PHOTO v5.3'

  xrange=[-0.5,0.5]
  yrange=[-0.5,0.5]
  psf_e1old = (struct[gind].psfixx[clr] - struct[gind].psfiyy[clr])/(struct[gind].psfixx[clr] + struct[gind].psfiyy[clr])
  psf_e1new = struct[gind].m_e1_psf[clr]
  aplot, 1, psf_e1old, psf_e1new, psym=3, title='PSF', xrange=xrange,yrange=yrange,xtit=xtitle,ytit=ytitle
  oplot, [-1000, 1000],[-1000,1000], color=!blue

  diff_e1_psf = psf_e1new-psf_e1old
  plothist, diff_e1_psf, bin=dbinsize, xrange=[-.05,.05],$
            xtitle = 'e!D2!N PHOTO - e!D2!N Michigan'
  oplot,[0,0],[0,1.e5],color=!blue

  xtitle = 'e!D2!N Michigan'
  ytitle = 'e!D2!N PHOTO v5.3'
  psf_e2old = 2.*struct[gind].psfixy[clr]/(struct[gind].psfixx[clr] + struct[gind].psfiyy[clr])
  psf_e2new = struct[gind].m_e2_psf[clr]
  aplot, 1, psf_e2old, psf_e2new, psym=3, title='PSF', xrange=xrange,yrange=yrange,xtit=xtitle,ytit=ytitle
  oplot, [-1000, 1000],[-1000,1000], color=!blue

  diff_e2_psf = psf_e2new-psf_e2old
  plothist, diff_e2_psf, bin=dbinsize, xrange=[-.05,.05],$
            xtitle = 'e!D2!N PHOTO - e!D2!N Michigan'
  oplot,[0,0],[0,1.e5],color=!blue
  print
  print,'PSF'
  print,'  Median diff_e2_psf = ',median(diff_e2_psf)
  ;;print,'Mean diff_e2_psf = ',mean(diff_e2_psf)
  print,'  Sdev diff_e2_psf = ',sdev(diff_e2_psf)
  print,'  Skewness diff_e2_psf = ',skewness(diff_e2_psf)
  !p.multi=0

END 
