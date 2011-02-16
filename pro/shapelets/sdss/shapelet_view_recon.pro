PRO shapelet_view_recon, image, recon, psFieldInfo, addnoise=addnoise

  IF n_params() LT 3 THEN BEGIN
      print,'-Syntax: shapelet_view_recon, image, recon, psFieldInfo, /addnoise'
      return
  ENDIF 

  pmulti_old = !p.multi
  !p.multi = [0,2,3]

  IF keyword_set(addnoise) THEN BEGIN 
      recon_backup = recon
      nrecon=n_elements(recon)
      recon[*] = $
        recon[*] + psFieldInfo.skysig*randomu(seed, nrecon, /normal)
  ENDIF 

  diff = image-recon

  siglevels = [3,5,10,20,50]
  c_lines=[0,0,0,0,0]

  IF (!d.name EQ 'Z') OR (!d.name EQ 'PS') THEN BEGIN 

      IF !d.name EQ 'Z' THEN BEGIN 
          setupplot,'Z'
          device, set_resolution=[770,890]
      ENDIF 
      loadct, 0
      c_colors = [!p.color,!p.color,!p.color,100,100]

  ENDIF ELSE IF !d.name EQ 'X' THEN BEGIN 
      setupplot,'X'
      c_colors = [!white, !green, !red, !blue, !darkGreen]
  ENDIF 

  IF !d.name EQ 'X' OR !d.name EQ 'Z' THEN BEGIN 
      charold = !p.charsize
      !p.charsize=2.0
  ENDIF 

  skysig = psFieldInfo.skysig

  levels = siglevels*skysig

  tvasinh, image, title='original'
  tvasinh, recon, title='recon'

  image_contour, image, levels=levels, c_colors=c_colors
  image_contour, recon, levels=levels, c_colors=c_colors

  legend,ntostr(siglevels)+' '+!csym.sigma,$
    /right,box=0,charsize=0.7, lines=c_lines, colors=c_colors

  IF keyword_set(addnoise) THEN BEGIN 
      recon = 0
      recon = temporary(recon_backup)
  ENDIF 

  image_contour, diff, levels=levels, c_colors=c_colors, title='difference'

  IF !d.name EQ 'X' OR !d.name EQ 'Z' THEN !p.charsize=charold

  !p.multi = pmulti_old 

END 
