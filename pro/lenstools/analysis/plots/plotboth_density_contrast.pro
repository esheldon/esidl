PRO plotboth_density_contrast, struct, struct2, $
                               logbin=logbin, $
                               $
                               miny=miny, maxy=maxy, $
                               center=center, $
                               wuse=wuse, $
                               shear=shear, $
                               clr=clr, $
                               $
                               position1=position1, position2=position2,$
                               xwindow1=xwindow1,ywindow1=ywindow1,$
                               xwindow2=xwindow2,ywindow2=ywindow2, $
                               botyrange=botyrange, yrange=yrange, $
                               $
                               xfit=xfit, yfit=yfit, fitsym=fitsym, $
                               xoplot=xoplot, yoplot=yoplot, $
                               oploterr=oploterr, oplotcolor=oplotcolor, $
                               oplotsym=oplotsym, rebin=rebin, _extra=_extra

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: plotboth_density_contrast, struct, struct2, $'
      print,'                  /logbin, $'
      print,'                  $'
      print,'                  miny=miny, maxy=maxy, $'
      print,'                  /center, $'
      print,'                  wuse=wuse, $'
      print,'                  shear=shear, $'
      print,'                  clr=clr, $'
      print,'                  $'
      print,'                  position1=position1, position2=position2,$'
      print,'                  xwindow1=xwindow1,ywindow1=ywindow1,$'
      print,'                  xwindow2=xwindow2,ywindow2=ywindow2, $'
      print,'                  botyrange=botyrange, yrange=yrange, $'
      print,'                  $'
      print,'                  xfit=xfit, yfit=yfit, fitsym=fitsym, $'
      print,'                  xoplot=xoplot, yoplot=yoplot, $'
      print,'                  oploterr=oploterr, oplotcolor=oplotcolor, $'
      print,'                  oplotsym=oplotsym, rebin=rebin, _extra=_extra'
      return
  ENDIF 

  ;; This one for plotting two at once, second either random or ortho
  ;; depending on if struct2 is sent
  ;; fixed to square on top
  ;; lower plot is always linear on y-axis

  IF n_elements(struct2) NE 0 THEN random=1 ELSE random=0

  psym=1
  minfac=0.5
  maxfac=1.5
  IF keyword_set(noxtitle) THEN xtitle='' ELSE xtitle=!kpcxtitle2
  IF keyword_set(noytitle) THEN BEGIN
      sigytitle=''
      sigytitle2=''
  ENDIF ELSE BEGIN 
          
      IF NOT keyword_set(shear) THEN BEGIN 
          sigytitle=!deltaytitle
          IF keyword_set(random) THEN BEGIN
              sigytitle2=!csym.delta_cap+!csym.sigma_cap+$
                '!S!D'+!csym.plus+'!N!R!Urand!N'
          ENDIF ELSE BEGIN 
              sigytitle2=!csym.delta_cap+!csym.sigma_cap+'!D'+!csym.times+'!N'
          ENDELSE 
      ENDIF ELSE BEGIN 
          sigytitle=!shytitle
          IF keyword_set(random) THEN BEGIN
              sigytitle2=!csym.gamma+'!S!D'+!csym.plus+'!N!R!Urand!N'
          ENDIF ELSE BEGIN 
              sigytitle2=!orthoytitle
          ENDELSE 
      ENDELSE 

  ENDELSE 
  IF keyword_set(shear) THEN BEGIN 
      tmp=tag_exist(struct, 'shear', ind=tag)
      tmp=tag_exist(struct, 'shearerr', ind=errtag)
      IF random THEN BEGIN 
          tmp=tag_exist(struct2, 'shear', ind=tag2)
          tmp=tag_exist(struct2, 'shearerr', ind=errtag2)
      ENDIF ELSE BEGIN 
          tmp=tag_exist(struct, 'ortho', ind=tag2)
          tmp=tag_exist(struct, 'orthoerr', ind=errtag2)
      ENDELSE 

      mfac = mean(struct.shear/struct.sigma)
      shfac = 1.e3
  ENDIF ELSE BEGIN
      mfac=1.0
      shfac=1.0

      IF keyword_set(rebin) THEN restr = '_rebin' ELSE restr=''
      name = 'sigma'+restr
      ename = 'sigmaerr'+restr
      oname = 'orthosig'+restr
      oename = 'orthosigerr'+restr
      rname = 'meanr'+restr

      tmp=tag_exist(struct, name, ind=tag)
      tmp=tag_exist(struct, ename, ind=errtag)
      tmp=tag_exist(struct, rname, ind=rtag)
      IF random THEN BEGIN 
          tmp=tag_exist(struct2, name, ind=tag2)
          tmp=tag_exist(struct2, ename, ind=errtag2)
          tmp=tag_exist(struct2, rname, ind=rtag2)
      ENDIF ELSE BEGIN 
          tmp=tag_exist(struct, oname, ind=tag2)
          tmp=tag_exist(struct, oename, ind=errtag2)
          rtag2=rtag
      ENDELSE 

  ENDELSE 

  IF n_elements(yfit) NE 0 THEN BEGIN 
      yfitsend = yfit*shfac
  ENDIF 
  IF n_elements(miny) EQ 0 THEN miny=0.01*mfac*shfac ;Msolar/pc^2
  IF n_elements(maxy) EQ 0 THEN maxy=10000.*mfac*shfac

  IF n_elements(wuse) EQ 0 THEN BEGIN
      IF keyword_set(logbin) THEN w=where(struct.(rtag) GT 0.0) $
      ELSE w=lindgen(n_elements(struct.(rtag)))
  ENDIF ELSE w=wuse

  topaspect = 1
  fract = 0.8
  IF keyword_set(logbin) THEN BEGIN 
      xrange=[minfac*min(struct.(rtag)[w]), maxfac*max(struct.(rtag)[w])]

      IF xrange[0] LT 100 THEN BEGIN 
          IF xrange[0] GT 10 THEN xrange[0] = 10.0
      ENDIF 

      xticklen=0.04 & yticklen=0.04
      xlog=1 & ylog=1
      square=1

      w2=where(struct.(tag)[w] GT 0.0,nw)

      yr=prange(shfac*struct.(tag)[w[w2]], $
                    shfac*struct.(errtag)[w[w2]],slack=0.5)
      yr[0] = yr[0] > miny
      yr[1] = yr[1] < maxy

      IF random THEN BEGIN 
          yr2=prange(shfac*struct2.(tag2)[w], $
                         shfac*struct2.(errtag2)[w],/slack,/symmetric)
      ENDIF ELSE BEGIN 
          yr2=prange(shfac*struct.(tag2)[w],$
                         shfac*struct.(errtag2)[w],/slack,/symmetric)
      ENDELSE 

      xtickf = 'loglabels'
      ytickf = 'loglabels'

  ENDIF ELSE BEGIN 
      yr=prange(shfac*struct.(tag)[w], $
                shfac*struct.(errtag)[w],/slack)

      yr[0] = yr[0] < 0

      IF random THEN BEGIN 
          yr2=prange(shfac*struct2.(tag2)[w], $
                         shfac*struct2.(errtag2)[w],/slack,/symmetric)
      ENDIF ELSE BEGIN 
          yr2=prange(shfac*struct.(tag2)[w],$
                         shfac*struct.(errtag2)[w],/slack,/symmetric)
      ENDELSE 


      xrange=[0, 1.1*max(struct.(rtag)[w])]

  ENDELSE 
  
  IF n_elements(xfit) NE 0 THEN BEGIN 

      xoplot = xfit
      yoplot = yfitsend
      IF n_elements(fitsym) NE 0 THEN oplotsym = fitsym

  ENDIF

  IF n_elements(yrange) EQ 0 THEN yrange=yr
  IF n_elements(botyrange) EQ 0 THEN botyrange=yr2

  myusersym, 'fill_circle'
  psym=8
  IF random THEN BEGIN 
      plottwoerr, struct.(rtag)[w], $
                  struct.(tag)[w]*shfac, struct.(errtag)[w]*shfac,$
                  struct2.(rtag2)[w], $
                  struct2.(tag2)[w]*shfac,struct2.(errtag2)[w]*shfac,$
                  xtitle=xtitle, $
                  topytitle=sigytitle, $
                  botytitle=sigytitle2,$
                  botyminor=5,$
                  topyrange=yrange, botyrange=botyrange,$
                  xrange=xrange,$
                  xlog=xlog,topylog=ylog,$
                  toppsym=psym, $
                  botpsym=8,$
                  ystyle=1, xstyle=1, center=center,$
                  xticklen=xticklen, yticklen=yticklen, frac=0.8, topaspect=1,$
                  xoplot=xoplot, yoplot=yoplot,oploterr=oploterr,oplotcolor=oplotcolor,$
                  oplotsym=oplotsym, $
                  position1=position1, position2=position2,$
                  xwindow1=xwindow1,ywindow1=ywindow1,$
                  xwindow2=xwindow2,ywindow2=ywindow2, xtickformat=xtickf,ytickformat1=ytickf,$
                  _extra=_extra

      IF n_elements(clr) NE 0 THEN BEGIN 
          xwold = !x.window & !x.window=xwindow1
          ywold = !y.window & !y.window=ywindow1
          
          legend, !colorsp[clr], /right, box=0, charsize=!p.charsize*1.3
          !x.window=xwold & !y.window=ywold
      ENDIF 
  ENDIF ELSE BEGIN 
      plottwoerr, struct.(rtag)[w], $
                  struct.(tag)[w]*shfac, struct.(errtag)[w]*shfac,$
                  struct.(rtag)[w], $
                  struct.(tag2)[w]*shfac,struct.(errtag2)[w]*shfac,$
                  xtitle=xtitle, $
                  topytitle=sigytitle, $
                  botytitle=sigytitle2,$
                  botyminor=5,$
                  topyrange=yrange, botyrange=botyrange,$
                  xrange=xrange,$
                  xlog=xlog,topylog=ylog,$
                  toppsym=psym, $
                  botpsym=8,$
                  ystyle=1, xstyle=1, center=center,$
                  xticklen=xticklen, yticklen=yticklen, frac=0.8, topaspect=1,$
                  xoplot=xoplot, yoplot=yoplot,oploterr=oploterr,oplotcolor=oplotcolor,$
                  oplotsym=oplotsym, $
                  position1=position1, position2=position2,$
                  xwindow1=xwindow1,ywindow1=ywindow1,$
                  xwindow2=xwindow2,ywindow2=ywindow2, xtickformat=xtickf,ytickformat1=ytickf,$
                  _extra=_extra

      IF n_elements(clr) NE 0 THEN BEGIN 
          xwold = !x.window & !x.window=xwindow1
          ywold = !y.window & !y.window=ywindow1
          
          legend, !colorsp[clr], /right, box=0, charsize=!p.charsize*1.3
          !x.window=xwold & !y.window=ywold
      ENDIF 
  ENDELSE 
  oplot,[0.01,1.e6], [0,0]

END 
