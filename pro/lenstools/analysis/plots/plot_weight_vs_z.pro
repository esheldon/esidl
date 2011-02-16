PRO plot_weight_vs_z, struct, weight, xo, yo, sig, hist, $
                      dops=dops, encap=encap, lrg=lrg

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: plot_weight_vs_z, struct, weight, xo, yo, sig, hist, /dops, /encap, /lrg'
      return
  ENDIF 

  IF keyword_set(lrg) THEN addstr='_lrg' ELSE addstr=''

  asciifile = '~/lensout/weights/weight_vs_z'+addstr+'.dat'
  fitfile = '~/lensout/weights/weight_vs_z'+addstr+'.fit'
  psfile = '~/plots/weights/weight_vs_z'+addstr
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

  IF keyword_set(lrg) THEN BEGIN 
      zbinsize = 0.01
      zmin = 0.15
      zmax = 0.6
      yrange = [0, 0.005]
  ENDIF ELSE BEGIN 
      zbinsize = 0.01
      zmin = 0.02
      zmax = 0.3
      yrange = [0,0.013]
  ENDELSE 

  binner, struct.z, weight, zbinsize, xo, yo, sig, rev, hist, $
          max=zmax

  yo = smooth(yo, 3)
  xtit = 'z'
  ytit = 'W(z)'
  aplot, !gratio, xo, yo, psym=8,xtit=xtit,ytit=ytit, $
         yrange=yrange,/ysty, xrange=[0,zmax],/xsty,$
         xticklen=0.03

  oplot, [xo[0],xo], [0,hist*0.004/max(hist)],psym=10

  st = create_struct('z', xo, $
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
