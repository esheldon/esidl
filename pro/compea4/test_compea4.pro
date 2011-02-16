PRO test_compea4, struct, w, R_out

  IF n_elements(struct) EQ 0 THEN BEGIN 

      read_tsobj, [2662, 20, 3], struct, /corr, /all

  ENDIF 

  clr = 2

;  w=where(struct.m_R[clr] GT 0.0 AND struct.m_R[clr] LT 0.8 AND $
;          struct.petrocounts[clr] LT 21. AND struct.petrocounts[clr] GT 18., nw)
  w=where(struct.m_R[clr] GT 0.0 AND $
          struct.petrocounts[clr] LT 21. AND struct.petrocounts[clr] GT 18., nw)

  compea4_struct, struct[w], clr, e1_out, e2_out, R_out, flags
  help,struct,where(flags NE 0)

  e1_corr = double(struct[w].m_e1_corr[clr])/(1. - struct[w].m_R[clr])
  e2_corr = double(struct[w].m_e2_corr[clr])/(1. - struct[w].m_R[clr])

  !p.multi=[0,0,2]

  color = !green

  xrange = [-1, 1]
  yrange = [-1, 1]
  aplot, 1, e1_corr, e1_out, psym=3, xrange=xrange, yrange=yrange,$
         xtit='e1 old', ytit='e1 new'
  oplot, [-10, 10], [-10, 10], color=color
  aplot, 1, e2_corr, e2_out, psym=3, xrange=xrange, yrange=yrange,$
         xtit='e2 old', ytit='e2 new'
  oplot, [-10, 10], [-10, 10], color=color

  key=get_kbrd(1)


  R = struct[w].m_R[clr]

  dil_corr = 1./(1.-R)
  dil_corr_out = 1./(1.-R_out)

  percinc = (dil_corr_out - dil_corr)/dil_corr

  aplot,1,R,R_out,psym=3,xrange=[0,1],yrange=[0,1],$
        xtit='R old',ytit='R new'
  oplot,[0,10],[0,10],color=color
  aplot,1,R,percinc,psym=3,xrange=[0,1],yrange=[0,1],$
        xtit='R old', ytit='(corr_new-corr_old)/corr_old'

  !p.multi=0


END 
