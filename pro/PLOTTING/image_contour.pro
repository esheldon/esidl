PRO image_contour, image, x, y, proc=proc, levels=levels, $
        aspect=aspect, $    
        nocenter=nocenter, $
        nsmooth=nsmooth, $
        c_labels=c_labels, $
        c_colors=c_colors, $
        scale=scale, $
        xrange=xrange, yrange=yrange, $
        sky=sky, high=high, $
        ntick=ntick, $
        xtickn1=xtickn1, xtickn2=xtickn2, $
        ytickn1=ytickn1, ytickn2=ytickn2, $
        xtitle1=xtitle1, xtitle2=xtitle2, $
        ytitle1=ytitle1, ytitle2=ytitle2, $
        position=position, _extra=e

    IF n_params() LT 1 THEN BEGIN 
        print,'Syntax: image_contour, image, ...'
        return
    ENDIF 

    IF n_elements(ntick) EQ 0 THEN ntick = 0

    IF ntick THEN style = 4+1 ELSE style = 1

    IF n_elements(proc) EQ 0 THEN proc = 'tvasinh'

    noframe = 0
    IF n_elements(x) EQ 0 THEN BEGIN 
        call_procedure, proc, image, sky=sky, high=high, scale=scale, $
            aspect=aspect, nocenter=nocenter, $
            xrange=xrange, yrange=yrange
        contour, $
            image, $
            levels=levels, $
            c_labels=c_labels, c_colors=c_colors, $
            /overplot, $
            ystyle=style,xstyle=style, $
            _extra=e, position=position
    ENDIF ELSE BEGIN 
        call_procedure, proc, image, noframe=noframe, position=position, $
            aspect=aspect, nocenter=nocenter, $
            sky=sky, high=high, xrange=xrange, yrange=yrange
        if n_elements(nsmooth) ne 0 then begin
            contour, $
                smooth(image,nsmooth), x, y,$
                levels=levels, $
                c_labels=c_labels, c_colors=c_colors, $
                /overplot, $
                ystyle=style,xstyle=style, $
                _extra=e, position=position
 
        endif else begin
            contour, $
                image, x, y,$
                levels=levels, $
                c_labels=c_labels, c_colors=c_colors, $
                /overplot, $
                ystyle=style,xstyle=style, $
                _extra=e, position=position
        endelse
    ENDELSE 
;  tvasinh, image, /noframe, position=position, sky=sky, high=high
;  contour, image, x, y, levels=levels, /noerase, ystyle=style,xstyle=style,$
;    position=position, _extra=e

  return

  IF ntick THEN BEGIN 

      tmp = strarr(ntick)
      tmp[*] = ' '
      IF n_elements(xtickn1) EQ 0 THEN xtickn1 = tmp
      IF n_elements(xtickn2) EQ 0 THEN xtickn2 = tmp
      IF n_elements(ytickn1) EQ 0 THEN ytickn1 = tmp
      IF n_elements(ytickn2) EQ 0 THEN ytickn2 = tmp
      IF n_elements(xtitle1) EQ 0 THEN xtitle1 = ' '
      IF n_elements(xtitle2) EQ 0 THEN xtitle2 = ' '
      IF n_elements(ytitle1) EQ 0 THEN ytitle1 = ' '
      IF n_elements(ytitle2) EQ 0 THEN ytitle2 = ' '

      axis, xaxis=0, xticks=ntick-1, xtickn=xtickn1, xtitle=xtitle1
      axis, xaxis=1, xticks=ntick-1, xtickn=xtickn2, xtitle=xtitle2
      axis, yaxis=0, yticks=ntick-1, ytickn=ytickn1, ytitle=ytitle1
      axis, yaxis=1, yticks=ntick-1, ytickn=ytickn2, ytitle=ytitle2
  ENDIF 

  return

END 
