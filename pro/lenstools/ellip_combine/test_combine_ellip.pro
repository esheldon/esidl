PRO test_combine_ellip, struct, $
                              e1, e2, e1e1err, e1e2err, e2e2err, $
                              new_e1, new_e2, $
                              new_e1e1err, new_e1e2err, new_e2e2err,$
                              new_smear, $
                              combine_flag, good, redo=redo

  IF n_elements(struct) EQ 0 THEN BEGIN 
      
      
      
      tags = ['m_e1','m_e2','m_e1e1err', 'm_e1e2err','m_e2e2err',$
              'm_e1_corr', 'm_e2_corr',$
              'm_rr_cc', 'm_rr_ccerr', $
              'm_cr4','m_r','petrocounts','reddening']
      
      run = 1458
      rerun = 20
      FOR camcol=1,6 DO BEGIN 
          read_tsobj, [run,rerun,camcol], tstruct, /all, /corr, $
                      taglist=tags,/add_sky

          IF n_elements(struct) EQ 0 THEN BEGIN
              struct=temporary(tstruct) 
          ENDIF ELSE BEGIN
              concat_structs, temporary(struct),temporary(tstruct),tmp
              struct=temporary(tmp)
          ENDELSE 

      ENDFOR 
      return
  ENDIF 

  begplot,name='/net/cheops2/home/esheldon/plots/compare_combine_ellip/combine3.ps',$
          /color

  rcut = 1.2
  rminmag = 18.0
  rmaxmag = 21.5

  rmag = struct.petrocounts[2] - struct.reddening[2]

  w=where(rmag GE rminmag AND rmag LE rmaxmag AND $
         struct.m_r[2] GT 0.0 AND struct.m_r[2] LT rcut)


  e1 = struct[w].m_e1_corr[1:3]
  e2 = struct[w].m_e2_corr[1:3]
 
  e1e1err = struct[w].m_e1e1err[1:3]
  e1e2err = struct[w].m_e1e2err[1:3]
  e2e2err = struct[w].m_e2e2err[1:3]

  smear = struct[w].m_r[1:3]

  IF n_elements(new_e1) EQ 0 OR keyword_set(redo) THEN BEGIN 
      combine_ellip_cove1e2, e1, e2, e1e1err, e1e2err, e2e2err, smear, $
                             new_e1, new_e2, $
                             new_e1e1err, new_e1e2err, new_e2e2err,$
                             new_smear, $
                             combine_flag, good
  ENDIF 

  ;; objects for which there were 3 good measurements
  
  w2=where( (combine_flag AND 2b^0) NE 0 AND $
            (combine_flag AND 2b^1) NE 0 AND $
            (combine_flag AND 2b^2) NE 0, nw2)

  !p.multi=[0,0,2]

  loadct,0
  ploth,smear[1,w2],new_smear[w2],xrange=[0,1.2],yrange=[0,1.2],/tvim,$
        xtitle='R[2]', ytitle='R combined'
  simpctable
  oplot,[-10,10],[-10,10],color=!blue

  loadct,0
  ploth,smear[1,w2],new_smear[w2]-smear[1,w2],xrange=[0,1.2],yrange=[-0.2,0.2],/tvim,$
        xtitle='R[2]', ytitle='R combined '+!csym.minus+' R[2]'
  simpctable
  oplot,[-10,10],[0,0],color=!blue

  !p.multi=0

  endplot


END 
