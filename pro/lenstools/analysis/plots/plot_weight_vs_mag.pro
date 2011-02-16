PRO plot_weight_vs_mag, struct, weight, xo, yo, sig, hist, binsize=binsize, $
                        dops=dops, encap=encap, lrg=lrg

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: plot_weight_vs_z, struct, weight, xo, yo, sig, his, binsize=binsize, /dops, /encap, /lrg'
      return
  ENDIF 

  IF keyword_set(lrg) THEN addstr='_lrg' ELSE addstr=''

  asciifile = '~/lensout/weights/weight_vs_mag'+addstr+'.dat'
  fitfile = '~/lensout/weights/weight_vs_mag'+addstr+'.fit'
  psfile = '~/plots/weights/weight_vs_mag'+addstr
  IF keyword_set(encap) THEN psfile = psfile+'.eps' ELSE psfile=psfile+'.ps'

  IF keyword_set(dops) THEN BEGIN 
      IF keyword_set(encap) THEN BEGIN 
          begplot,name=psfile,xsize=7, ysize=7,/encap
      ENDIF ELSE BEGIN 
          begplot,name=psfile
      ENDELSE 
  ENDIF
  nlens = n_elements(struct)
  
  weight = struct.weight

  

  mmax_act = max(struct.absmag[2], min=mmin_act)

  IF keyword_set(lrg) THEN BEGIN 
      mmin = (-24) > mmin_act
      mmax = (-21) < mmax_act
      yrange = [0, 0.005]
      IF n_elements(binsize) EQ 0 THEN binsize = 0.15
  ENDIF ELSE BEGIN 
      mmin = (-24) > mmin_act
      mmax = (-17) < mmax_act
      yrange=[0,0.014]
      IF n_elements(binsize) EQ 0 THEN binsize = 0.2
  ENDELSE 
  binner, struct.absmag[2], weight, binsize, xo, yo, sig, rev, hist, $
          max=mmax,min=mmin

;  nperbin = 3000
;  binner_bynum, struct.absmag[2], weight, nperbin, xo, yo, sig, hist, $
;                xmin = mmin, xmax = mmax

  ;; smooth a bit
  yo = smooth(yo, 3)

  ;; only use a subset
  nbin = n_elements(xo)
  wuse = lindgen(nbin)
  wuse = wuse[0:nbin-3]
  xo = xo[wuse]
  yo = yo[wuse]
  sig = sig[wuse]
  hist = hist[wuse]

  xtit = 'M!Dr!N - 5 log(h)'
  ytit = 'W(M!Dr!N)'
  aploterror, !gratio, xo, yo, sig, psym=8,xtit=xtit,ytit=ytit, $
         yrange=yrange,/ysty,xrange=[mmin,mmax],$
         xticklen=0.03

  oplot, [xo[0],xo], [0,hist*0.004/max(hist)],psym=10

  st = create_struct('absmag_r', xo, $
                     'weight', yo, $
                     'num', hist)

  print,'Writing file: ',asciifile
  colprint,xo,yo,hist,file=asciifile

  print,'Writing file: ',fitfile
  mwrfits, st, fitfile, /create
  
  IF keyword_set(dops) THEN endplot
  IF keyword_set(encap) THEN set_bbox, psfile, '%%BoundingBox: 20 15 480 345'

return

END 
