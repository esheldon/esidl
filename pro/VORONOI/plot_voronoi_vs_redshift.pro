PRO plot_voronoi_vs_redshift, lcat

  w1=where(lcat.z1d LT 0.15,nw1)
  w2=where(lcat.z1d LT 0.15 AND lcat.voronoi_dens NE 0.0,nw2)
  DD = angdist_lambda(lcat.z1d, h=0.7, omega=0.3)

  plot,DD[w2], lcat[w2].ra*!dpi/180d, /polar, psym=3,/iso,/nodata,$
       xstyle=4,ystyle=4

  R = 255L
  ming = 20L
  maxg = 200L
  G = long(arrscl(findgen(nw2), ming, maxg))
  B = 20L
  oranges = R + 256L*(G+256L*B) 

  ;dark_orange = 255L + 256L*(50L + 256L*0L)
  ;oranges = arrscl(findgen(nw2), 

  myusersym,'fill_circle'
  s=sort(lcat[w2].voronoi_dens)
  FOR i=0L, nw2-1 DO BEGIN 
      ind=w2[s[i]]
      oplot, [DD[ind]], [lcat[ind].ra*!dpi/180d], $
             color=oranges[i],psym=8,/polar,symsize=0.3
  ENDFOR 

END 
