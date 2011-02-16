PRO shcorr_plt,e1psf,psfile=psfile

IF n_params() EQ 0 THEN BEGIN
  print,'-Syntax: shcorr_plt,e1psf,psfile=psfile'
  return
ENDIF

n=20
emeas = 2.0*(findgen(n)-n/2.0)/n
p = findgen(n)/n

title = 'Corrected e1        e1 psf = '+strtrim(string(e1psf),2)
xtitle='P smear polarizability'
ytitle = 'corrected e1'
yrange=[-1,1]

FOR i=0,n-1 DO BEGIN
  index = n/1.5
  egal = (emeas[i] - e1psf*p)/(1-p)
  IF i EQ 0 THEN BEGIN
    plot,p,egal,title=title,xtitle=xtitle,ytitle=ytitle,yrange=yrange
    IF abs(egal[index]) LT 1.0 THEN BEGIN 
      out='em = '+strmid( strtrim(string(emeas[i]),2), 0, 5)
      xyouts,p[index],egal[index]-.03,out,align=.5
    ENDIF
  ENDIF ELSE BEGIN
    oplot,p,egal
    IF abs(egal[index]) LT 1.0 THEN BEGIN
      out='em = '+strmid( strtrim(string(emeas[i]),2), 0, 5)
      xyouts,p[index],egal[index]-.03,out,align=.5
    ENDIF
  ENDELSE
ENDFOR


IF keyword_set(psfile) THEN BEGIN
  makeps,'/sdss4/data1/esheldon/correction.ps',/noland
  FOR i=0,n-1 DO BEGIN
    index = n/1.5
    egal = (emeas[i] - e1psf*p)/(1-p)
    IF i EQ 0 THEN BEGIN
      plot,p,egal,title=title,xtitle=xtitle,ytitle=ytitle,yrange=yrange
      IF abs(egal[index]) LT 1.0 THEN BEGIN 
        out='em = '+strmid( strtrim(string(emeas[i]),2), 0, 5)
        xyouts,p[index],egal[index]-.03,out,align=.5
      ENDIF 
    ENDIF ELSE BEGIN
      oplot,p,egal
      IF abs(egal[index]) LT 1.0 THEN BEGIN 
        out='em = '+strmid( strtrim(string(emeas[i]),2), 0, 5)
        xyouts,p[index],egal[index]-.03,out,align=.5
      ENDIF 
    ENDELSE
  ENDFOR
  ep
ENDIF


return
END
