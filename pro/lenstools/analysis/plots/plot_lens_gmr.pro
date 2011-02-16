PRO plot_lens_gmr, struct, dops=dops, encap=encap
  
  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: plot_lens_gmr, struct, /dops, /encap'
      return
  ENDIF 
  
  IF keyword_set(dops) THEN BEGIN 
      dir = '~/plots/'
      front = 'gmr_dist_bylum'
      IF keyword_set(encap) THEN BEGIN 
          file = dir+front+'.eps'
          begplot,name=file,/color,xsize=11,ysize=8.5,/encap
      ENDIF ELSE BEGIN 
          file = dir+front+'.ps'
          begplot,name=file,/color,/land
      ENDELSE 
  ENDIF 

  erase & multiplot, [3,2], /square
  
  read_cuts, range_struct
  minz = 0.02
  maxz = 0.3

  gmrmin = 0.1
  gmrmax = 1.1
  gmr = struct.abscounts[1] - struct.abscounts[2]

  weights = struct.weight

  FOR clr=0,4 DO BEGIN 

      lowcut = range_struct.main_minmag_threebin[clr,0]
      highcut = range_struct.main_maxmag_threebin[clr,0]
      w1 = where(struct.absmag[clr] le highcut AND $
                 struct.absmag[clr] gt lowcut  AND $
                 struct.z gt minz AND $
                 struct.z lt maxz AND $
                 gmr GT gmrmin AND $
                 gmr LT gmrmax)
      
      
      
      lowcut = range_struct.main_minmag_threebin[clr,1]
      highcut = range_struct.main_maxmag_threebin[clr,1]
      w2 = where(struct.absmag[clr] le highcut AND $
                 struct.absmag[clr] gt lowcut  AND $
                 struct.z gt minz AND $
                 struct.z lt maxz AND $
                 gmr GT gmrmin AND $
                 gmr LT gmrmax)
      
      
      
      lowcut = range_struct.main_minmag_threebin[clr,2]
      highcut = range_struct.main_maxmag_threebin[clr,2]
      w3 = where(struct.absmag[clr] le highcut AND $
                 struct.absmag[clr] gt lowcut  AND $
                 struct.z gt minz AND $
                 struct.z lt maxz AND $
                 gmr GT gmrmin AND $
                 gmr LT gmrmax)

      lowcut = range_struct.main_minmag_threebin[clr,2]
      highcut = range_struct.main_maxmag_threebin[clr,0]
      
      wall = where(struct.absmag[clr] le highcut AND $
                   struct.absmag[clr] gt lowcut  AND $
                   struct.z gt minz AND $
                   struct.z lt maxz AND $
                   gmr GT gmrmin AND $
                   gmr LT gmrmax)
      
      mgmr1 = total(weights[w1]*gmr[w1])/total(weights[w1])
      mgmr2 = total(weights[w2]*gmr[w2])/total(weights[w2])
      mgmr3 = total(weights[w3]*gmr[w3])/total(weights[w3])
      mgmr  = total(weights[wall]*gmr[wall])/total(weights[wall])
      
      bin = 0.01
      clr1 = !green
      clr2 = !blue
      clr3 = !red
      
      xtitle = 'g '+!csym.minus+' r'
      ytitle = 'dN/d('+xtitle+')'
      
      CASE clr OF 
          0: xtitle=''
          1: BEGIN 
              xtitle='' 
              ytitle='' 
          END 
          2: ytitle=''
          4: ytitle=''
          ELSE: 
      ENDCASE 
          

      norm=1
      fline=1
      IF !d.name EQ 'X' THEN BEGIN 
          thick=!p.thick*2 
      ENDIF ELSE BEGIN 
          thick=2
          !p.thick=2
      ENDELSE 
      
      plothist, gmr[wall], xhist, yhist, bin=bin, min=gmrmin, max=gmrmax,norm=norm, /noplot
      plothist, gmr[w1], xhist1, yhist1, bin=bin, min=gmrmin, max=gmrmax,norm=norm, /noplot
      plothist, gmr[w2], xhist2, yhist2, bin=bin, min=gmrmin, max=gmrmax,norm=norm, /noplot
      plothist, gmr[w3], xhist3, yhist3, bin=bin, min=gmrmin, max=gmrmax,norm=norm, /noplot
      
      maxy = max([max(yhist),max(yhist1),max(yhist2),max(yhist3)])
      yrange = [0, maxy]
      yrange = [0, 12.5]

;  plothist, gmr[wall], bin=bin, min=gmrmin, max=gmrmax,norm=norm,/fill,position=aspect(1./!gratio),$
;            xtitle=xtitle,ytitle=ytitle, yrange=yrange

      IF !d.name EQ 'PS' THEN fspacing = 0.1
      plothist, gmr[w1], bin=bin, min=gmrmin, max=gmrmax,norm=norm,/fill,$
                xtitle=xtitle,ytitle=ytitle, yrange=yrange,fline=fline,forientation=0,$
                xstyle=1+2,ystyle=1,charsize=1.2,xticklen=0.04,yticklen=0.04,fspacing=fspacing

      
      
      ;;oplot,[mgmr,mgmr],[0,1.e6]
;  if !d.name eq 'X' then key=get_kbrd(1)
;  plothist, gmr[w1], bin=bin, /overplot, color=clr1,norm=norm,$
;            thick=thick,/fill,fcolor=clr1,fline=fline,forientation=90
      ;;oplot,[mgmr1,mgmr1],[0,1.e6],color=clr1
      ;;if !d.name eq 'X' then key=get_kbrd(1)
      plothist, gmr[w2], bin=bin, /overplot, color=clr2,norm=norm,$
                thick=thick,/fill,fcolor=clr2,fline=fline,forientation=45,fspacing=fspacing
      ;;oplot,[mgmr2,mgmr2],[0,1.e6],color=clr2
      ;;if !d.name eq 'X' then  key=get_kbrd(1)
      plothist, gmr[w3], bin=bin, /overplot, color=clr3,norm=norm,$
                thick=thick,/fill,fcolor=clr3,fline=fline,forientation=135,fspacing=fspacing
      ;;oplot,[mgmr3,mgmr3],[0,1.e6],color=clr3

      legend,!colors[clr],box=0,/right,charsize=1.5

      CASE clr OF
          0: multiplot
          1: multiplot,/doxaxis
          2: multiplot
          3: multiplot
          4: multiplot,/reset
      ENDCASE 

      distfile = '~/lensout/gmrdist/'+!colors[clr]+'_gmr_dist.dat'
      colprint, xhist, yhist, yhist1, yhist2, yhist3, file=distfile

;      print
;      colprint,xhist,xhist1,xhist2,xhist3
;      print

  ENDFOR 

  IF keyword_set(dops) THEN BEGIN 
      endplot
      IF keyword_set(encap) THEN BEGIN 
          set_bbox, file, '%%BoundingBox: 35 15 770 515'
      ENDIF ELSE BEGIN 
          pslandfix, file
      ENDELSE 
  ENDIF 
  
END 
