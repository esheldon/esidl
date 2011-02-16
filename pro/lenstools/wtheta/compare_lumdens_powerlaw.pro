PRO compare_lumdens_powerlaw, twobin=twobin

  indir = '/sdss5/data0/lensout/mass2light/'
  lall = mrdfits(indir+'lumdense_allplot.fit',1,/silent)
  lspiral = mrdfits(indir+'lumdense_spiral.fit',1,/silent)
  lellip = mrdfits(indir+'lumdense_ellip.fit',1,/silent)

  IF keyword_set(twobin) THEN ttstr='twobin' ELSE ttstr=''
  llum1 = mrdfits(indir+'lumdense_lumnum1'+ttstr+'.fit',1,/silent)
  llum2 = mrdfits(indir+'lumdense_lumnum2'+ttstr+'.fit',1,/silent)
  IF NOT keyword_set(twobin) THEN BEGIN 
      llum3 = mrdfits(indir+'lumdense_lumnum3'+ttstr+'.fit',1,/silent)
      llum4 = mrdfits(indir+'lumdense_lumnum4'+ttstr+'.fit',1,/silent)
  ENDIF 
  indir2 = '/sdss5/data0/lensout/stripe10/'

  all = mrdfits(indir2+'main_zgal_gal_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_fitparam_N1.fit',1,/silent)
  spiral=mrdfits(indir2+'spiral_zgal_gal_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_fitparam_N1.fit',1,/silent)
  ellip=mrdfits(indir2+'ellip_zgal_gal_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_fitparam_N1.fit',1,/silent)

  lum1=mrdfits(indir2+'sublum/i/lum1'+ttstr+'_zgal_gal_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_fitparam_N1.fit',1,/silent)
  lum2=mrdfits(indir2+'sublum/i/lum2'+ttstr+'_zgal_gal_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_fitparam_N1.fit',1,/silent)
  IF NOT keyword_set(twobin) THEN BEGIN 
      lum3=mrdfits(indir2+'sublum/i/lum3'+ttstr+'_zgal_gal_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_fitparam_N1.fit',1,/silent)
      lum4=mrdfits(indir2+'sublum/i/lum4'+ttstr+'_zgal_gal_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_fitparam_N1.fit',1,/silent)
  ENDIF 

  ;; make plots with alpha on x, norm on y

  erase & multiplot, [2,3]

  yrange=[0,10]
  xrange=[-2,0]
  IF keyword_set(twobin) THEN BEGIN 
      xrange=[-1.5,0]
      yrange=[0,8]
  ENDIF 
  FOR clr=0, 4 DO BEGIN 



      IF clr GT 0 THEN BEGIN 
          IF clr EQ 3 THEN multiplot,/doxaxis ELSE multiplot
          ;multiplot
      ENDIF 

      yy = [yrange[0], yrange[0], yrange[1], yrange[1]]
      
      plot, [0], /nodata, yrange=yrange, xrange=xrange, $
        xtitle=xtitle
      
      
      oplot, [all.powlow[0], all.powlow[0]], yrange,thick=!p.thick*2.
      oplot, [all.powhigh[0], all.powhigh[0]], yrange,thick=!p.thick*2.
      oplot, [spiral.powlow[0], spiral.powlow[0]], yrange,color=!blue
      oplot, [spiral.powhigh[0], spiral.powhigh[0]], yrange,color=!blue
      oplot, [ellip.powlow[0], ellip.powlow[0]], yrange,color=!red
      oplot, [ellip.powhigh[0], ellip.powhigh[0]], yrange,color=!red
      
      oplot, [lum1.powlow[0], lum1.powlow[0]], yrange,color=!cyan,line=2
      oplot, [lum1.powhigh[0], lum1.powhigh[0]], yrange,color=!cyan,line=2
      oplot, [lum2.powlow[0], lum2.powlow[0]], yrange,color=!green,line=2,$
        thick=!p.thick*2.
      oplot, [lum2.powhigh[0], lum2.powhigh[0]], yrange,color=!green,line=2,$
        thick=!p.thick*2.
      IF NOT keyword_set(twobin) THEN BEGIN 
          oplot, [lum3.powlow[0], lum3.powlow[0]], yrange,color=!magenta,line=2
          oplot, [lum3.powhigh[0], lum3.powhigh[0]], yrange,color=!magenta,line=2
          oplot, [lum4.powlow[0], lum4.powlow[0]], yrange,color=!yellow,line=2
          oplot, [lum4.powhigh[0], lum4.powhigh[0]], yrange,color=!yellow,line=2
      ENDIF 
      ;IF clr EQ 5 THEN legend,'Mass Power Laws',/left,box=0
      
      pname = !colors[clr]+'pow'
      nname = !colors[clr]+'norm'
      
      IF (clr EQ 0) OR (clr EQ 2) OR (clr EQ 4) THEN $
        ytitle = 'A (10!U10!N L'+sunsymbol()+' Mpc!U'+!tsym.minus+'2!N)' $
      ELSE ytitle=''
      IF clr EQ 3 OR clr EQ 4 THEN xtitle = !tsym.alpha ELSE xtitle=''
      plotderror,lall.pow[clr], lall.norm[clr], $
        lall.powlow[clr], lall.powhigh[clr],$
        lall.normlow[clr], lall.normhigh[clr], xrange=xrange,yrange=yrange,$
        xtitle=xtitle,ytitle=ytitle
      
      oplotderror, lspiral.pow[clr], lspiral.norm[clr], $
        lspiral.powlow[clr], lspiral.powhigh[clr],$
        lspiral.normlow[clr], lspiral.normhigh[clr], xrange=xrange,yrange=yrange,$
        color=!blue,errcolor=!blue
      oplotderror, lellip.pow[clr], lellip.norm[clr], $
        lellip.powlow[clr], lellip.powhigh[clr],$
        lellip.normlow[clr], lellip.normhigh[clr], xrange=xrange,yrange=yrange,$
        color=!red,errcolor=!red
      
      oplotderror, llum1.pow[clr], llum1.norm[clr], $
        llum1.powlow[clr], llum1.powhigh[clr],$
        llum1.normlow[clr], llum1.normhigh[clr], xrange=xrange,yrange=yrange,$
        color=!yellow,errcolor=!cyan
      oplotderror, llum2.pow[clr], llum2.norm[clr], $
        llum2.powlow[clr], llum2.powhigh[clr],$
        llum2.normlow[clr], llum2.normhigh[clr], xrange=xrange,yrange=yrange,$
        color=!green,errcolor=!green
      IF NOT keyword_set(twobin) THEN BEGIN 
          oplotderror, llum3.pow[clr], llum3.norm[clr], $
            llum3.powlow[clr], llum3.powhigh[clr],$
            llum3.normlow[clr], llum3.normhigh[clr], xrange=xrange,yrange=yrange,$
            color=!magenta,errcolor=!magenta
          oplotderror, llum4.pow[clr], llum4.norm[clr], $
            llum4.powlow[clr], llum4.powhigh[clr],$
            llum4.normlow[clr], llum4.normhigh[clr], xrange=xrange,yrange=yrange,$
            color=!cyan,errcolor=!yellow
      ENDIF 
      
      legend, !colorsp[clr], /right, box=0

      types =  ['all',$
                'spiral',$
                'ellip',$
                'lum1',$
                'lum2']
      colors= [!p.color,$
               !blue,$
               !red,$
               !cyan,$
               !green]
      IF NOT keyword_set(twobin) THEN BEGIN
          types = [types,'lum3','lum4']
          colors=[colors, !magenta, !yellow]
      ENDIF 
      nnn = n_elements(types)
      IF display_type() EQ 'X' THEN clear=0 ELSE clear=1
      IF clr EQ 0 THEN BEGIN 
          legend, types, $
            color=colors, $
            line=replicate(0,nnn),thick=replicate(!p.thick,nnn),clear=clear
      ENDIF 
      
      IF display_type() EQ 'X' THEN key=get_kbrd(1)

  ENDFOR 

  multiplot,/reset


END   
