PRO compea4_struct, struct, clr, e1_out, e2_out, R_out, flags

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: compea4_struct, struct, clr, e1_out, e2_out, R_out, flags'
      return
  ENDIF 

  ntot = n_elements(struct)

  defval = -9999d

  e1_out = replicate(defval, ntot)
  e2_out = e1_out
  R_out = e1_out

  flags = lonarr(ntot)

  defval2 = -1000.
  w=where(struct.m_e1[clr] NE defval2, Ngood, comp=out, ncomp=nbad)

  IF ngood EQ 0 THEN BEGIN 
      print,'No good measurements'
      return
  ENDIF 

  e1 = struct[w].m_e1[clr]
  e2 = struct[w].m_e2[clr]

  e1_psf = struct[w].m_e1_psf[clr]
  e2_psf = struct[w].m_e2_psf[clr]

  ;; factor makes a differenc on rho4 but not 
  ;; on m_rr_cc
  ;;fac = 2.0
  m_rr_cc = struct[w].m_rr_cc[clr]
  m_rr_cc_psf = struct[w].m_rr_cc_psf[clr]

  rho4 = struct[w].m_cr4[clr]
  rho4_psf = struct[w].m_cr4_psf[clr]

  ;compea4old, e1, e2, e1_psf, e2_psf, m_rr_cc, m_rr_cc_psf, $
  ;  rho4, rho4_psf, te1_out, te2_out, tR_out, tflags
	
	compea4, $
		e1, e2, m_rr_cc, rho4, $
		e1_psf, e2_psf, m_rr_cc_psf, rho4_psf, $
		te1_out, te2_out, tR_out, tflags

  e1_out[w] = te1_out
  e2_out[w] = te2_out
  R_out[w] = tR_out
  flags[w] = tflags
  
  ;; last flag set in compea4 is 2L^5
  IF nbad NE 0 THEN flags[out] = flags[out] + 2L^6

END 
