
PRO test_combine_ellip_error, struct, $
                              e1, e2, e1e1err, e1e2err, e2e2err, $
                              new_e1, new_e2, $
                              new_e1e1err, new_e1e2err, new_e2e2err,$
                              new_smear, $
                              combine_flag, good

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

  begplot,name='/net/cheops2/home/esheldon/plots/compare_combine_ellip/errors.ps',$
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

  combine_ellip_cove1e2, e1, e2, e1e1err, e1e2err, e2e2err, smear, $
                         new_e1, new_e2, $
                         new_e1e1err, new_e1e2err, new_e2e2err,$
                         new_smear, $
                         combine_flag, good


  ;; objects for which there were 3 good measurements
  
  w2=where( (combine_flag AND 2b^0) NE 0 AND $
            (combine_flag AND 2b^1) NE 0 AND $
            (combine_flag AND 2b^2) NE 0, nw2)

  !p.multi=[0,0,2]

  ratio = new_e1e1err[w2]/e1e1err[1,w2]
  ratmed = median(ratio)
  plothist, ratio, bin=0.05, min=0, max=2, $
            title='3 bandpasses averaged', $
            xtitle='new_e1e1err / old_e1e1err'
  oplot, [1./sqrt(3.),1./sqrt(3.)], [0, 1.e6], line=0
  oplot, [ratmed,ratmed], [0, 1.e6], line=2
  legend, [ntostr(ratmed),'1/'+!csym.sqrt+'3 = '+ntostr(1./sqrt(3.))], $
          /right, box=0,charsize=1,line=[2,0]

  ratio = new_e2e2err[w2]/e2e2err[1,w2]
  ratmed = median(ratio)
  plothist, ratio, bin=0.05, min=0, max=2, $
            xtitle='new_e2e2err / old_e2e2err'
  oplot, [1./sqrt(3.),1./sqrt(3.)], [0, 1.e6], line=0
  oplot, [ratmed,ratmed], [0, 1.e6], line=2
  legend, [ntostr(ratmed),'1/'+!csym.sqrt+'3 = '+ntostr(1./sqrt(3.))], $
          /right, box=0,charsize=1,line=[2,0]

  ;; multiply 
  plothist, new_e1e1err[w2]/e1e1err[1,w2]*sqrt(3.), bin=0.05, min=0, max=2, $
            title='3 bandpasses averaged', $
            xtitle='new_e1e1err / old_e1e1err * '+!csym.sqrt + '3'
  oplot, [1.,1.], [0, 1.e6], line=0
  plothist, new_e2e2err[w2]/e2e2err[1,w2]*sqrt(3.), bin=0.05, min=0, max=2, $
            xtitle='new_e2e2err / old_e2e2err * '+!csym.sqrt + '3'
  oplot, [1.,1.], [0, 1.e6], line=0

  ;; objects for which there were only 2 good measurements: r band
  ;; and another

  w2=where( ((combine_flag AND 2b^0) NE 0 AND $
             (combine_flag AND 2b^1) NE 0 AND $
             (combine_flag AND 2b^2) EQ 0) OR $
            $
            ((combine_flag AND 2b^0) NE 0 AND $
             (combine_flag AND 2b^1) EQ 0 AND $
             (combine_flag AND 2b^2) NE 0), nw2)

  ratio = new_e1e1err[w2]/e1e1err[1,w2]
  ratmed = median(ratio)
  plothist, ratio, bin=0.05, min=0, max=2, $
            title='2 bandpasses averaged', $
            xtitle='new_e1e1err / old_e1e1err'
  oplot, [1./sqrt(2.),1./sqrt(2.)], [0, 1.e6], line=0
  oplot, [ratmed,ratmed], [0, 1.e6], line=2
  legend, [ntostr(ratmed),'1/'+!csym.sqrt+'2 = '+ntostr(1./sqrt(2.))], $
          /right, box=0,charsize=1,line=[2,0]

  ratio = new_e2e2err[w2]/e2e2err[1,w2]
  ratmed = median(ratio)
  plothist, ratio, bin=0.05, min=0, max=2, $
            xtitle='new_e2e2err / old_e2e2err'
  oplot, [1./sqrt(2.),1./sqrt(2.)], [0, 1.e6], line=0
  oplot, [ratmed,ratmed], [0, 1.e6], line=2
  legend, [ntostr(ratmed),'1/'+!csym.sqrt+'2 = '+ntostr(1./sqrt(2.))], $
          /right, box=0,charsize=1,line=[2,0]

  ;; multiply 
  plothist, new_e1e1err[w2]/e1e1err[1,w2]*sqrt(2.), bin=0.05, min=0, max=2, $
            title='2 bandpasses averaged', $
            xtitle='new_e1e1err / old_e1e1err * '+!csym.sqrt + '2'
  oplot, [1.,1.], [0, 1.e6], line=0
  plothist, new_e2e2err[w2]/e2e2err[1,w2]*sqrt(2.), bin=0.05, min=0, max=2, $
            xtitle='new_e2e2err / old_e2e2err * '+!csym.sqrt + '2'
  oplot, [1.,1.], [0, 1.e6], line=0

  !p.multi=0

  endplot



END 
