PRO plot_voronoi, cat

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: plot_voronoi, cat'
      return
  ENDIF 

  home = '/sdss3/usrdevel/esheldon/'
  htmldir = home+'WWW/Voronoi/'
  htmlfile = htmldir+'run752_756_foreground.html'
  
  ;openw, outfile, htmlfile, /get_lun
  ;printf, outfile, '<HTML>'
  ;printf, outfile, '    <HEAD><TITLE>'+htmlfile+'</TITLE></HEAD>'
  ;printf, outfile, '    <BODY bgcolor="#ffffff" link="#0066ff" vlink="#FF0000" text="#00000">'

;cat must be sorted by ra for fast results

  step = 9.0                    ;degrees
  
  x=cat.dec
  y=cat.ra

  indices = lindgen(n_elements(cat))

  
  y1=min(y)
  y2=y1+9.0
  w=where(y GE y1 AND y LT y2, nw)
  
  i=0
  outmessage='</PRE>'
  ;begplot,name=htmldir+'run752_756_foreground.ps'
  WHILE (nw NE 0) DO BEGIN

      IF i EQ 0 THEN BEGIN 
          voronoi_density, x[indices[w]], y[indices[w]], /plot, /true
          first=0
      ENDIF ELSE BEGIN
          voronoi_density, x[indices[w]], y[indices[w]], /plot, /true, xstyle=4
      ENDELSE 
      i=i+1
      gifname = 'test_N'+ntostr(i)+'.gif'
;      write_gif, htmldir+gifname, tvrd()
      outmessage='<img src="'+gifname+'" alt="Voronoi Tesselation"><br>'+outmessage

      key=get_kbrd(1)
      IF key EQ 'q' THEN return

      IF max(indices[w]) NE max(indices) THEN BEGIN 
          remove, w, indices

          y1 = min(y[indices])
          y2 = y1 + 9.0
          w=where(y[indices] GE y1 AND y[indices] LT y2, nw)
      ENDIF ELSE nw = 0
  ENDWHILE 
  ;ep
  ;printf, outfile, '    <PRE>'+outmessage
  ;printf, outfile, '    </BODY>'
  ;printf, outfile, '</HTML>'
  ;close, outfile
  ;free_lun, outfile

return
END 
