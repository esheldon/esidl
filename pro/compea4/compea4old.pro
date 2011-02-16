PRO compea4old, e1, e2, e1_psf, e2_psf, m_rr_cc, m_rr_cc_psf, $
             rho4, rho4_psf, e1_out, e2_out, R_out, flags

  IF n_params() LT 8 THEN BEGIN 
      print,'-Syntax:  compea4, e1, e2, e1_psf, e2_psf, m_rr_cc, m_rr_cc_psf,$'
      print,'             rho4, rho4_psf [, e1_out, e2_out, R_out, flags]'
  ENDIF 

  Ndata = n_elements(e1)

  e1 = double(e1)
  e2 = double(e2)
  e1_psf = double(e1_psf)
  e2_psf = double(e2_psf)
  m_rr_cc = double(m_rr_cc)
  m_rr_cc_psf = double(m_rr_cc_psf)
  rho4 = double(rho4)
  rho4_psf = double(rho4_psf)


  e1_out = replicate(-9999d, Ndata)
  e2_out = replicate(-9999d, Ndata)
  R_out  = replicate(-9999d, Ndata)

  flags  = lonarr(Ndata)

  compea4_sofile = esheldon_config("compea4_sofile")
  compea4_entry = esheldon_config("compea4_entry")

  retval = call_external(value=[0b, 0b, 0b, 0b, 0b, 0b, $
                                0b, 0b, 0b, 0b, 0b, 0b, 0b], $
                         compea4_sofile, compea4_entry, $
                         m_rr_cc, m_rr_cc_psf, $
                         e1, e2, e1_psf, e2_psf, $
                         rho4, rho4_psf, Ndata, $
                         e1_out, e2_out, R_out, flags)

END 
