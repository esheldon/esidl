;; Budavari photozs

PRO make_photoz_struct, struct

  val = -9999.
  val2 = -9999
  struct = create_struct('unique_id', ulong64(val2), $
                         'chisq', val, $
                         'z', val, $
                         'zErr', val, $
                         't', val, $
                         'tErr', val, $
                         'covar_tt', val, $
                         'covar_tz', val, $
                         'covar_zz', val, $
                         'diff_idxZ', val2, $
                         'thres', val, $
                         'quality', val2, $
                         'distmod', val, $
                         'abs_umg', val, $
                         'abs_gmr', val, $
                         'abs_rmi', val, $
                         'abs_imz', val, $
                         'kcorr_u', val, $
                         'kcorr_g', val, $
                         'kcorr_r', val, $
                         'kcorr_i', val, $
                         'kcorr_z', val, $
                         'abscounts_u', val, $
                         'abscounts_g', val, $
                         'abscounts_r', val, $
                         'abscounts_i', val, $
                         'abscounts_z', val, $
                         'rmag', val, $
                         'ra', val, $
                         'dec', val) ; They used floats for ra/dec !! why!!??!!

END 
