PRO combine_like, files, mean=mean, err=err, wmean=wmean, werr=werr, input_denscont=input_denscont

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: combine_like, files, mean=mean, err=err, wmean=wmean, werr=werr'
      return
  ENDIF 

  !p.multi=[0,0,2]

  setup_mystuff


  xtitle = !deltaytitle
  ytitle = 'Likelihood'


  nf = n_elements(files)


  IF !d.name EQ 'X' THEN BEGIN 
      R = 255L
      B = 0L
      G = long( arrscl( lindgen(nf), 0, 255 ) )
      colors = R + 256L*(G+256L*B)
      charsize = 1
  ENDIF ELSE BEGIN 
      colors = long( arrscl( lindgen(nf), !grey0, !grey90 ) )
      charsize = 0.7
  ENDELSE 
;  colors = long( arrscl( lindgen(nf), 255, 65280 ) )

  nd = 1000L

  FOR i=0L, nf-1 DO BEGIN 

      print,'reading file: ',files[i]
      t = mrdfits(files[i], 1,/silent)

      IF i EQ 0 THEN BEGIN 
          mind = min(t.denscont)
          maxd = max(t.denscont)
          denscont = arrscl(findgen(nd), mind, maxd)

          input_denscont = t.input_denscont

          like = interpol(t.denscont_likelihood, t.denscont, denscont)
          like = like/max(like)
          loglike = alog(like)

          plot, denscont, like,yrange=[0, 1.2], xtitle=xtitle, ytitle=ytitle
          oplot, denscont, like, color=colors[i]
 
          

          IF tag_exist(t, 'wmean_s') THEN BEGIN
              add_arrval, t.wmean_s, wmean_s
              add_arrval, t.werr_s, werr_s
          ENDIF 
      ENDIF ELSE BEGIN 
          IF t.input_denscont NE input_denscont THEN message,'Different input density contrasts'
          
          like = interpol(t.denscont_likelihood, t.denscont, denscont)
          like = like/max(like)
          loglike = loglike + alog(like)

          oplot, denscont, like, color=colors[i]

          IF tag_exist(t, 'wmean_s') THEN BEGIN
              add_arrval, t.wmean_s, wmean_s
              add_arrval, t.werr_s, werr_s
          ENDIF 
      ENDELSE 
  ENDFOR 

  like = (exp(1d))^(loglike - max(loglike))
  npts = 100
  norm = qgauss(like, denscont, npts)

  mean = qgauss(denscont*like, denscont, npts)/norm
  var  = qgauss((denscont-mean)^2*like, denscont, npts)/norm
  err = sqrt(var)

  nsig = 4.5
  xrange = [mean-nsig*err, mean+nsig*err]

  myusersym, 'fill_circle'
  plot, denscont, like, yrange=[0, 1.2], xrange=xrange, $
        xtitle=xtitle, ytitle=ytitle,xstyle=1
  oplot, denscont, like, psym=8

  nterms = 3
  guess = [1.0, mean, err]
  w=where(denscont GE xrange[0] AND denscont LE xrange[1])
  yfit = gaussfit(denscont[w], like[w], A, estimates=guess, nterms=3)
  oplot, denscont[w], yfit, color=!red

  oplot, [input_denscont, input_denscont], [0, 10000], color=!green

  print
  print,'Mean of likelihood: '+ntostr(mean)+' '+$
        !plusminus+' '+ntostr(err)+$
        '  (Nsig = '+ntostr( (mean-input_denscont)/err )+')'
  print,'Mean of Gausian:    '+ntostr(A[1])+' '+$
        !plusminus+' '+ntostr(A[2])+$
        '  (Nsig = '+ntostr( (A[1]-input_denscont)/A[2] )+')'

  IF n_elements(wmean_s) EQ nf THEN BEGIN 
      print
      print,'All had wmean_s, calculating mean'
      
      weights = 1d/werr_s^2
      wtot = total(weights)
      wmean = total(wmean_s*weights)/wtot
      werr = sqrt(1./wtot)

      print,'wmean =             '+$
            ntostr(wmean)+' '+!plusminus+' '+ntostr(werr)+$
        '  (Nsig = '+ntostr( (wmean-input_denscont)/werr )+')'
      
      wfit = gaussprob(denscont, wmean, werr)
      wfit = wfit/max(wfit)
      oplot, denscont, wfit, color=!cyan

      legend, ['Mean   = '+ntostr(mean,6,/round)+$
               !csym.plusminus+ntostr(err,6,/round), $
               'Wmean = '+ntostr(wmean,6,/round)+$
               !csym.plusminus+ntostr(werr,6,/round)],charsize=charsize
      legend,['Best Fit gaussian', 'Wmean'], $
             colors=[!red, !cyan], line=[0,0],/right,charsize=charsize

  ENDIF ELSE BEGIN 
      legend, 'Mean = '+ntostr(mean,6,/round)+$
              !csym.plusminus+ntostr(err,6,/round),charsize=charsize
      legend,['Best Fit gaussian'], $
             colors=[!red], line=[0],/right,charsize=charsize
  ENDELSE 

  !p.multi=0

END 
