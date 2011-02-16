PRO shecter, nbins, func, minL, maxL, Lum=Lum, aeta=aeta, alpha=alpha, h=h, $
             plot=plot

  IF n_params() EQ 0 THEN BEGIN
    print,'-Syntax: shecter, nbins, func, [minL, maxL, Lum=Lum, aeta=aeta, alpha=alpha, h=h, plot=plot]'
    print,' Returns shecter function (#/cubic Mpc) in nbins'
    print,' Optionally returns L in units of L*'
    return
  ENDIF


  ; parameters for shecter function.
  IF n_elements(minL) EQ 0 THEN minL = .02
  IF n_elements(maxL) EQ 0 THEN maxL = 2.
  IF NOT keyword_set(h) THEN h = .7
  IF NOT keyword_set(aeta) THEN aeta  = 1.2e-2*h^3 ; per cubic Mpc
  IF NOT keyword_set(alpha) THEN alpha = 1.15


  ;scaled step and luminosity
  step = (maxL - minL)/(nbins-1)
  Lum =  step*dindgen(nbins) + minL
 
  func = aeta*Lum^(-alpha)*exp(-Lum)*step
  

  IF keyword_set(plot) THEN BEGIN
    xtitle='L/L*'
    ytitle='#/ (cubic Mpc)
    title='Shecter Function'
    plot,Lum, func,xtitle=xtitle,ytitle=ytitle,title=title
  ENDIF 
  return
END 
