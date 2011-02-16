PRO test_vsold, old, new, camcol, clr, mold, mnew, wold, wnew, e1, e2, diff_e1, diff_e2, cutzero=cutzero

  ;; only need the "new" to see how many total are thrown out
  ;; via various flags

  IF n_elements(old) EQ 0 THEN BEGIN 

      tags = ['ra','dec',$
              'ixx','iyy','ixy','rho4','momerr',$
              'psfixx','psfiyy','psfixy', 'psfrho4',$
              'e1','e2','petrocounts','counts_model','flags','flags2','objc_flags']

      rerun=9
      fetch_dir, 756, camcol, rerun, dir, corrdir=corrdir
      read_tsobj, corrdir, old, /all, taglist=tags, front='adatc'

  ENDIF 

  IF n_elements(new) EQ 0 THEN BEGIN

      tags = ['ra','dec',$
              'm_rr_cc','m_rr_ccerr','m_cr4',$
              'm_rr_cc_psf','m_cr4_psf',$
              'm_e1','m_e2','m_e1e1err', 'm_e1e2err', 'm_e2e2err',$
              'petrocounts','counts_model','flags','flags2','objc_flags']

      indir = '/net/cheops3/home/esheldon/test_photo_v5.3/756/9/calibChunks/'+$
        ntostr(camcol)+'/'

      read_tsobj, indir, new, /all, taglist=tags

  ENDIF 

  IF n_elements(mold) EQ 0 OR n_elements(mnew) EQ 0 THEN BEGIN 

      ;; get bright ones that pass flag cuts
      ;; NEW
      nnew = n_elements(new)
      print,'Making cuts on new struct'
      wbright=where(new.petrocounts[2] LE 22,nbright)
      ;;help,new,wbright
      print,'Removed '+ntostr(nnew - nbright)+'/'+ntostr(nnew)+' = '+$
            ntostr((nnew - nbright)/float(nnew))+' faint objects'
      
      make_flag_struct,fs
      fs.amoment_faint = 'N'
      flag_select, new[wbright], fs, clr, si_faint
      nkeep = n_elements(si_faint)
      print,'Removed '+ntostr(nbright - nkeep)+'/'+ntostr(nbright)+' = '+$
            ntostr((nbright-nkeep)/float(nbright))+' amoment_faint'
      
      make_flag_struct,fs
      fs.amoment_shift = 'N'
      flag_select, new[wbright], fs, clr, si_shift
      nkeep = n_elements(si_shift)
      print,'Removed '+ntostr(nbright - nkeep)+'/'+ntostr(nbright)+' = '+$
            ntostr((nbright-nkeep)/float(nbright))+' amoment_shift'

      delvarx,fs
      make_flag_struct,fs
      fs.amoment_maxiter = 'N'
      flag_select, new[wbright], fs, clr, si_maxiter
      ;;help,new[wbright],si_maxiter
      nkeep = n_elements(si_maxiter)
      print,'Removed '+ntostr(nbright - nkeep)+'/'+ntostr(nbright)+' = '+$
            ntostr((nbright-nkeep)/float(nbright))+' amoment_maxiter'
      
      
      fs.amoment_faint = 'N'
      fs.amoment_shift = 'N'
      fs.amoment_maxiter = 'N'
      flag_select, new[wbright], fs, clr, si_amomentcheck
      ;;help,new[wbright],si_amomentcheck
      nkeep = n_elements(si_amomentcheck)
      print,'Removed '+ntostr(nbright - nkeep)+'/'+ntostr(nbright)+' = '+$
            ntostr((nbright-nkeep)/float(nbright))+' total bad amoment'

      gnew = wbright[si_amomentcheck]

      IF keyword_set(cutzero) THEN BEGIN 
          keep = where(new[gnew].m_e1[clr] NE 0.0,nkeep2)
          print,'Removed '+ntostr(nkeep - nkeep2)+'/'+ntostr(nkeep)+' = '+$
                ntostr((nkeep-nkeep2)/float(nkeep))+' because e1 zero'
          gnew = gnew[keep]
      ENDIF 
      nkeep = n_elements(gnew)

      ;; now match to old way by photoid

      ;; OLD
;      nold = n_elements(new)
;      print
;      print,'Making cuts on old struct'
;      wbright = where(old.petrocounts[2] LE 22,nbright)
      ;;help,old,wbright
;      print,'Removed '+ntostr(nold - nbright)+'/'+ntostr(nold)+' = '+$
;            ntostr((nold - nbright)/float(nold))+' faint objects'
;      
;      gold = where(old[wbright].e1[clr] NE 1.e10)
;      gold = wbright[gold]
;      nkeep = n_elements(gold)
;      print,'Removed '+ntostr(nbright - nkeep)+'/'+ntostr(nbright)+' = '+$
;            ntostr((nbright-nkeep)/float(nbright))+' total bad amoment'
      
;      print
;      print,'Matching by (ra,dec)'
;      rad = 1d/3600d
;      allow = 1
;      close_match_radec, old[gold].ra, old[gold].dec, $
;                         new[gnew].ra, new[gnew].dec, mold,mnew,$
;                         rad, allow
      
;      mold = gold[mold]
;      mnew = gnew[mnew]
;      print


      sphoto_match, new[gnew], old, mnew,mold
      mnew = gnew[mnew]
      nmatch = n_elements(mnew)
      print,'  Matched '+ntostr(nmatch)+'/'+ntostr(nkeep)+' with "old" catalog'

  ENDIF 
  ;; plot various things new vs old as a function of magnitude

  ;; bin by 0.5 in mag
  binsize = 1
  rmin = 15.0
  rmax = 22.0
  nbin = 7

  magptr = ptrarr(nbin)

  hist = histogram(new[mnew].petrocounts[2], min=rmin, max=rmax,$
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

      wold = mold[*magptr[binnum]]
      wnew = mnew[*magptr[binnum]]

      e1 = (old[wold].ixx[clr] - old[wold].iyy[clr])/(old[wold].ixx[clr] + old[wold].iyy[clr])
      aplot, 1, e1, new[wnew].m_e1[clr], psym=3,$
            xtitle=xtitle, ytitle=ytitle, title=title, $
            xrange=xrange
      oplot, [-1000, 1000],[-1000,1000], color=!blue

      diff_e1 = new[wnew].m_e1[clr] - e1
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

      wold = mold[*magptr[binnum]]
      wnew = mnew[*magptr[binnum]]

      e2 = 2.*old[wold].ixy[clr]/(old[wold].ixx[clr] + old[wold].iyy[clr])
      aplot, 1, e2, new[wnew].m_e2[clr], psym=3,$
            xtitle=xtitle, ytitle=ytitle, title=title, $
            xrange=xrange
      oplot, [-1000, 1000],[-1000,1000], color=!blue

      diff_e2 = new[wnew].m_e2[clr] - e2
      plothist, diff_e2, bin=dbinsize, xrange=[-.05,.05],$
                xtitle = 'e!D2!N PHOTO - e!D2!N Michigan'

      print,title
      print,'  Median diff_e2 = ',median(diff_e2)
      ;;print,'Mean diff_e2 = ',mean(diff_e2)
      print,'  Sdev diff_e2 = ',sdev(diff_e2)
      print,'  Skewness diff_e2 = ',skewness(diff_e2)

  ENDFOR 

  ptrarr_free, magptr

END 
