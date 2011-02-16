FUNCTION clrarr2string, clr

  nclr = n_elements(clr)
  clr_string = ''
  FOR iclr=0L, nclr-1 DO clr_string = clr_string + !colors[clr[iclr]]

  return, clr_string

END 
