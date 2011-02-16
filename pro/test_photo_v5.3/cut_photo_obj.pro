PRO cut_photo_obj, struct, clr, gind, gal=gal, cutzero=cutzero

  ;; !!!!!!! ONLY CUTS ON PHOTO STUFF  !!!!!!!

  nstruct = n_elements(struct)
  print,'Making cuts on struct'
  gind=where(struct.petrocounts[2] LE 22,nbright)
  ;;help,new,wbright
  print,'Removed '+ntostr(nstruct - nbright)+'/'+ntostr(nstruct)+' = '+$
        ntostr((nstruct - nbright)/float(nstruct))+' faint objects'
  
  make_flag_struct, fs
  fs.amoment_faint = 'N'
  fs.amoment_shift = 'N'
  fs.amoment_maxiter = 'N'
  flag_select, struct[gind], fs, clr, si_amomentcheck
  ;;help,new[wbright],si_amomentcheck
  nkeep = n_elements(si_amomentcheck)
  print,'Removed '+ntostr(nbright - nkeep)+'/'+ntostr(nbright)+' = '+$
        ntostr((nbright-nkeep)/float(nbright))+' total bad amoment'
  
  gind = gind[si_amomentcheck]
  
  IF keyword_set(cutzero) THEN BEGIN 
      w=where(struct[gind].m_e1[clr] NE 0.0 AND struct[gind].m_e2[clr] NE 0.0,ngood) 
      gind = gind[w]
  
      print,'Removed '+ntostr(nkeep - ngood)+'/'+ntostr(nkeep)+' = '+$
            ntostr((nkeep-ngood)/float(nkeep))+' STILL ZERO!!'
  ENDIF 
  ;; get good galaxies
  IF keyword_set(gal) THEN BEGIN 
      corr = struct[gind].m_rr_cc_psf[clr]/struct[gind].m_rr_cc[clr]*(4./struct[gind].m_cr4_psf[clr] -1.)/(4./struct[gind].m_cr4[clr]-1.)
      gind_2 = where(corr LT 0.8 AND corr GT 0.0 AND $
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

  ;; now check our stuff
  ;;gind2 = where(struct[gind].e1[clr] NE 1.e10)
  ;;gind = gind[gind2]
  ;;nkeep = n_elements(gind)

END 
