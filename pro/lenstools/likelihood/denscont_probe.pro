FUNCTION denscont_probe, etan, ewidth,siginv,denscont

  ;; uses h=1.0, pc^2/Msolar
  model = 2.*denscont*siginv
  pofe = gaussprob(model, etan, ewidth)

  return, pofe

END 
