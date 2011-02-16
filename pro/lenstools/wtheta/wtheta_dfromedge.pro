PRO wtheta_dfromedge, llambda, leta, DL, etarange, dfrommax, dfrommin

  d2r = !dpi/180d

  maxeta = interpol(etarange.maxeta, etarange.lambda, llambda)
  mineta = interpol(etarange.mineta, etarange.lambda, llambda)

  gcirc, 0, leta*d2r, llambda*d2r, maxeta*d2r, llambda*d2r, dfrommax
  dfrommax=dfrommax*DL
  gcirc, 0, leta*d2r, llambda*d2r, mineta*d2r, llambda*d2r, dfrommin
  dfrommin=dfrommin*DL

END 
