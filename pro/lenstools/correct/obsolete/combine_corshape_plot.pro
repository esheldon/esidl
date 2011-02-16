PRO combine_corshape_plot_legend, intercept, intercepterr, slope, slopeerr, left=left, right=right

  sa=strmid(strtrim(string(intercept,format='(e8.1)'),2),0,8)
;  ssiga=strmid(strtrim(string(abs(intercepterr/intercept)),2),0,5)
  ssiga=strmid(strtrim(string(intercepterr,format='(e8.1)'),2),0,8)

  sb=strmid(strtrim(string(slope,format='(e8.1)'),2),0,8)
;  ssigb=strmid(strtrim(string(abs(slopeerr/slope)),2),0,5)
  ssigb=strmid(strtrim(string(slopeerr,format='(e8.1)'),2),0,8)

;  stta='a='+sa+' '+ssiga
;  sttb='b='+sb+' '+ssigb
  stta='a='+sa+!csym.plusminus+ssiga
  sttb='b='+sb+!csym.plusminus+ssigb

  ;;!p.background=!white
  legend, [stta,sttb], charsize = 0.7, /clear, left=left, right=right

END 

PRO combine_corshape_plot, struct, dolegend=dolegend, hirata=hirata, dops=dops


  rangefac = 0.5

  multold = !p.multi
  charold = !p.charsize
  symold = !p.symsize
  xmold = !x.margin
  ymold = !y.margin

;  !p.multi=[0,2,2,0,1]
  !p.charsize = 1.0
  symsize=0.75

  ;; defaults
;  !x.margin = [10.0,3.0]
;  !y.margin = [4.0, 2.0]

  !x.margin = [8.5, 0.5]
  !y.margin = [3.5, 0.15]

  yrange1 = [-0.055, 0.055]
  yrange2 = [-0.013, 0.013]
  xrange = [-0.14, 0.14]

;yrange1 = [-0.02,0.02]
;yrange2 = [-0.005, 0.005]
;xrange=[-0.02,0.02]

  xsize = 7
  ysize = 5

  xticklen=0.04
  yticklen=0.03

  IF keyword_set(hirata) THEN hirstr = '_h' ELSE hirstr=''
  IF keyword_set(dolegend) THEN legstr = '_leg' ELSE legstr=''
  
  file1 = '~/plots/e1gal_vs_e1psf'+hirstr+legstr+'.eps'
  file2 = '~/plots/e2gal_vs_e2psf'+hirstr+legstr+'.eps'

  IF keyword_set(dops) THEN begplot,name=file1, xsize=xsize,ysize=ysize,/encap

  maxx = max( abs(struct.tmean_e1psf) )
  xx=[-10.0*maxx,10.0*maxx]

  erase & multiplot,[1,2];, /square
  ;; e1-e1 plots
  xtitle = 'PSF e!D1!N'
  ytitle = 'Galaxy e!D1!N'
  fitlin, struct.tmean_e1psf, struct.tmean_e1gal, struct.terr_e1gal, $
          intercept, intercepterr, slope, slopeerr, /silent
  ploterror,struct.tmean_e1psf,struct.tmean_e1gal,struct.terr_e1gal,$
            psym=8, ytitle=ytitle, symsize=symsize, yrange=yrange1,/ystyle, xrange=xrange, /xsty, $
            xticklen=xticklen, yticklen=yticklen
  oplot, [-1,1],[0,0]
  oplot,xx,intercept + slope*xx, color=fcolor
  
  num = struct.hist_e1psf
  pnum = num*abs(yrange1[0])/max(num)*rangefac - abs(yrange1[0]) 
  s=sort(struct.tmean_e1psf)
  oplot, struct.tmean_e1psf[s], pnum[s] ,psym=10, color=!grey50

  IF keyword_set(dolegend) THEN BEGIN 
      combine_corshape_plot_legend, intercept, intercepterr, slope, slopeerr,/left
  ENDIF 

  multiplot
  ;; e1-e1 corrected plots
  xtitle = 'PSF e!D1!N'
  ytitle = 'Galaxy e!D1!N'
  fitlin, struct.tmean_e1psf, struct.tmean_e1galc, struct.terr_e1galc, $
          intercept, intercepterr, slope, slopeerr, /silent
  ploterror,struct.tmean_e1psf,struct.tmean_e1galc,struct.terr_e1galc,$
             psym=8, xtitle=xtitle, ytitle=ytitle, symsize=symsize, yrange=yrange2,/ystyle, xrange=xrange, xstyle=1, $
            xticklen=xticklen, yticklen=yticklen
  oplot, [-1,1],[0,0]
  oplot,xx,intercept + slope*xx, color=fcolor

  IF keyword_set(dolegend) THEN BEGIN 
      combine_corshape_plot_legend, intercept, intercepterr, slope, slopeerr,/right
  ENDIF 
  
  IF keyword_set(dops) THEN endplot 

  IF !d.name EQ 'X' THEN key=get_kbrd(1)

  IF keyword_set(dops) THEN begplot,name=file2, xsize=xsize,ysize=ysize,/encap
  multiplot,/reset
  erase & multiplot,[1,2];, /square
  ;; e2-e2 plots
  xtitle = 'PSF e!D2!N'
  ytitle = 'Galaxy e!D2!N'
  fitlin, struct.tmean_e2psf, struct.tmean_e2gal, struct.terr_e2gal, $
          intercept, intercepterr, slope, slopeerr, /silent
  ploterror,struct.tmean_e2psf,struct.tmean_e2gal,struct.terr_e2gal,$
            psym=8, ytitle=ytitle, symsize=symsize, yrange=yrange1,/ystyle, xrange=xrange, /xsty, $
            xticklen=xticklen, yticklen=yticklen
  oplot, [-1,1],[0,0]
  oplot,xx,intercept + slope*xx, color=fcolor

  s=sort(struct.tmean_e2psf)
  num = struct.hist_e2psf
  pnum = num*abs(yrange1[0])/max(num)*rangefac - abs(yrange1[0]) 
  oplot, struct.tmean_e2psf[s], pnum[s] ,psym=10, color=!grey50

  IF keyword_set(dolegend) THEN BEGIN 
      combine_corshape_plot_legend, intercept, intercepterr, slope, slopeerr,/left
  ENDIF 

  multiplot
  ;; e2-e2 corrected plots
  xtitle = 'PSF e!D2!N'
  ytitle = 'Galaxy e!D2!N'
  fitlin, struct.tmean_e2psf, struct.tmean_e2galc, struct.terr_e2galc, $
          intercept, intercepterr, slope, slopeerr, /silent
  ploterror,struct.tmean_e2psf,struct.tmean_e2galc,struct.terr_e2galc,$
            psym=8, xtitle=xtitle, ytitle=ytitle, symsize=symsize, yrange=yrange2,/ystyle, xrange=xrange, /xsty, $
            xticklen=xticklen, yticklen=yticklen
  oplot, [-1,1],[0,0]
  oplot,xx,intercept + slope*xx, color=fcolor

  IF keyword_set(dolegend) THEN BEGIN 
      combine_corshape_plot_legend, intercept, intercepterr, slope, slopeerr,/right
  ENDIF 
  multiplot,/reset

  IF keyword_set(dops) THEN endplot

  ;;!p.multi=multold
  !p.charsize = charold
  !x.margin = xmold
  !y.margin = ymold

END 
