PRO combine_like_wmean, files, mean=mean, err=err, input_denscont=input_denscont

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: combine_like_wmean, files, mean=mean, err=err, input_denscont=input_denscont'
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

          add_arrval, t.wmean_s, wmean_s
          add_arrval, t.werr_s, werr_s

          like = gaussprob(denscont, t.wmean_s, t.werr_s)
          like = like/max(like)

          plot, denscont, like,yrange=[0, 1.2], xtitle=xtitle, ytitle=ytitle
          oplot, denscont, like, color=colors[i]
 
      ENDIF ELSE BEGIN 
          IF t.input_denscont NE input_denscont THEN message,'Different input density contrasts'
          
          add_arrval, t.wmean_s, wmean_s
          add_arrval, t.werr_s, werr_s
          
          like = gaussprob(denscont, t.wmean_s, t.werr_s)
          like = like/max(like)

          oplot, denscont, like, color=colors[i]

      ENDELSE 
  ENDFOR 

  weights = 1d/werr_s^2
  wtot = total(weights)
  wmean = total(wmean_s*weights)/wtot
  werr = sqrt(1./wtot)
wmom,wmean_s,werr_s,wmean2,wsig2,werr2,/calcerr
  print,'wmean =             '+$
        ntostr(wmean)+' '+!plusminus+' '+ntostr(werr)+$
        '  (Nsig = '+ntostr( (wmean-input_denscont)/werr )+')'

  print,'wmean2 =             '+$
        ntostr(wmean2)+' '+!plusminus+' '+ntostr(werr2)+$
        '  (Nsig = '+ntostr( (wmean2-input_denscont)/werr2 )+')'

  nsig = 4.5
  xrange = [wmean-nsig*werr, wmean+nsig*werr]
  denscont = arrscl( findgen(1000), xrange[0], xrange[1] )
  like = gaussprob(denscont, wmean, werr)
  like = like/max(like)

  plot, denscont, like, yrange=[0, 1.2], xrange=xrange, $
        xtitle=xtitle, ytitle=ytitle,xstyle=1
  oplot, [input_denscont, input_denscont], [0, 1000], color=!green

  legend, 'Mean   = '+ntostr(wmean,6,/round)+$
          !csym.plusminus+ntostr(werr,6,/round),charsize=charsize

  !p.multi=0

END 
