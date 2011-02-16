PRO compare_counts_model, str

  help,str
  print

  simpctable, rd, gr, bl

  thickold = !p.thick
  !p.thick = 2.
  xthickold = !x.thick
  ythickold = !y.thick
  !x.thick = 2.
  !y.thick = 2.

  xcharsizeold = !x.charsize
  ycharsizeold = !y.charsize
  fontold = !p.font

  IF (!d.flags AND 1) EQ 0 THEN doX = 1 ELSE doX=0

  backold = !p.background
  IF NOT doX THEN BEGIN
      regcolor = !black
      !p.background = !white
      !x.charsize = 1.5
      !y.charsize = 1.5
      !p.font = 0
  ENDIF ELSE BEGIN 
      regcolor = !black
      !p.background = !white
      !x.charsize = 1.3
      !y.charsize = 1.3
  ENDELSE 
  
  maxmag = 22
  wdev=where(str.dev_l[0] GT str.exp_l[0] AND $
             str.dev_l[1] GT str.exp_l[1] AND $
             str.dev_l[2] GT str.exp_l[2] AND $
             str.dev_l[3] GT str.exp_l[3] AND $
             str.dev_l[4] GT str.exp_l[4] AND $
             $
             abs(str.counts_dev[2]) LT maxmag, ndev)

  print,ndev
  print
  IF ndev NE 0 THEN BEGIN 

      diff = str[wdev].counts_model - str[wdev].counts_dev

      help,diff
      print

      print,max( abs(diff[0,*]) )
      print,max( abs(diff[1,*]) )
      print,max( abs(diff[2,*]) )
      print,max( abs(diff[3,*]) )
      print,max( abs(diff[4,*]) )

      title = 'dev_l > exp_l        counts_dev[2] < '+ntostr(maxmag)

      plothist, diff[1,*],bin=.05,xrange=[-1,1],$
        title=title,xtitle='counts_model - counts_dev', color=regcolor
      plothist, diff[1,*],bin=.05,/overplot,color=!green
      plothist, diff[0,*],bin=.05,/overplot,line=3, color=!blue
      plothist, diff[4,*],bin=.05,/overplot,line=4, color=!magenta
      plothist, diff[3,*],bin=.05,/overplot,line=5, color=!red

      mess = ['u','g','i','z']
      line = [3,0,5,4]
      bcolor = regcolor
      textcolor = [!blue, !green, !red, !magenta]
      thick=replicate(!p.thick,4)

      legend,mess,line=line,thick=thick, $
        bcolor=regcolor,textcolor=textcolor, $
        color=textcolor

      IF doX THEN write_gif, 'color_diff_dev.gif',tvrd(), rd, gr, bl

  ENDIF 
  
  key=get_kbrd(1)

  wexp=where(str.exp_l[0] GT str.dev_l[0] AND $
             str.exp_l[1] GT str.dev_l[1] AND $
             str.exp_l[2] GT str.dev_l[2] AND $
             str.exp_l[3] GT str.dev_l[3] AND $
             str.exp_l[4] GT str.dev_l[4] AND $
             $
             abs(str.counts_exp[2]) LT maxmag, nexp)

  print,nexp
  print
  IF nexp NE 0 THEN BEGIN 

      diff = str[wexp].counts_model - str[wexp].counts_exp

      help,diff
      print

      print,max( abs(diff[0,*]) )
      print,max( abs(diff[1,*]) )
      print,max( abs(diff[2,*]) )
      print,max( abs(diff[3,*]) )
      print,max( abs(diff[4,*]) )

      title = 'dev_l > exp_l        counts_exp[2] < '+ntostr(maxmag)

      plothist, diff[3,*],bin=.05,xrange=[-1,1],line=5, $
        title=title, xtitle='counts_model - counts_exp',color=regcolor
      plothist, diff[3,*],bin=.05,/overplot,color=!red,line=5
      plothist, diff[0,*],bin=.05,/overplot,line=3,color=!blue
      plothist, diff[4,*],bin=.05,/overplot,line=4,color=!magenta
      plothist, diff[1,*],bin=.05,/overplot,line=0,color=!green

      mess = ['u','g','i','z']
      line = [3,0,5,4]
      bcolor = regcolor
      textcolor = [!blue, !green, !red, !magenta]
      thick=replicate(!p.thick,4)
      legend,mess,line=line,thick=thick, $
        bcolor=regcolor,textcolor=textcolor, $
        color=textcolor

      IF doX THEN write_gif, 'color_diff_exp.gif',tvrd(), rd, gr, bl

  ENDIF 
  
  !p.background = backold
  !p.thick = thickold
  !x.charsize = xcharsizeold 
  !y.charsize = ycharsizeold 
  !p.font = fontold

  !x.thick = xthickold 
  !y.thick = ythickold
return
END 
