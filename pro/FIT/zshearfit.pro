PRO zshearfit, files, shguess, siguess, yrange=yrange, $
               names=names, munit=munit, sigma=sigma, sigerr=sigerr, $
               mass=mass, merr=merr, meanr=meanr, str=str, $
               nowait=nowait, useall=useall, _extra=extra
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    ZSHEARFIT
;       
; PURPOSE:
;    Fits up to three different bandpasses to a model.  
;    Different from shearfit in that the data it reads is vs. physial radius
;
; CALLING SEQUENCE:
;    
;
; INPUTS: 
;    
;
; OPTIONAL INPUTS:
;    
;
; KEYWORD PARAMETERS:
;    
;       
; OUTPUTS: 
;    
;
; OPTIONAL OUTPUTS:
;    
;
; CALLED ROUTINES:
;    
; 
; PROCEDURE: 
;    
;	
;
; REVISION HISTORY:
;    
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF N_params() LT 1 THEN BEGIN 
     print,'-Syntax: zshearfit, files [, shguess, siguess, yrange=yrange,'
     print,'   names=names, munit=munit, sigma=sigma, sigerr=sigerr,'
     print,'   mass=mass, merr=merr, meanr=meanr, str=str, '
     print,'   nowait=nowait, _extra=extra]'
     print,''
     print,'shguess~[.005, -.7]         siguess~[40.,-.7]'
     print,'Use doc_library,""  for more help.'  
     return
  ENDIF 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Some parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(shguess) EQ 0 THEN BEGIN
      shguess=[.005, -.7]
      siguess=[40.,  -.7]
  ENDIF 

  IF keyword_set(useall) THEN useall = 1 ELSE useall = 0
  nfile = n_elements(files)
  IF nfile gt 3 THEN BEGIN
      print,'At most 3 bandpasses'
      return
  ENDIF 
  IF n_elements(munit) EQ 0 THEN munit = 1.e15
  IF n_elements(names) EQ 0 THEN names = strarr(nfile)
  IF NOT keyword_set(nowait) THEN nowait=0

  waittime = 10                  ;seconds
  xmarg = !x.margin
  !x.margin = [10,10]
  pold=!p.multi
  !p.multi=0

  ;; number of elements in strings
  nstring = 7
  ;; number from last x position to print xyouts
  xn = 4
  ;; Shear polarizeability correction
  Ssh = 1.15
  ;; convert from Msun/pc^2 to g/cm^2
  sigconvert = 2.2e-4

  ;; Conver to Mpc to prevent underflow
  fac = 1000.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; First do shear
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; find plot range
  FOR i=0, nfile-1 DO BEGIN
      rdzobjshear, files[i], str,/silent
      str.shear = str.shear*Ssh
      trange = prange(str.shear, str.shearerr)
      IF i EQ 0 THEN yrange = trange ELSE BEGIN
          yrange[0] = min([trange[0], yrange[0]])
          yrange[1] = max([trange[1], yrange[1]])
      ENDELSE 
  ENDFOR 

  xtitle = 'Projected Radius (1/h kpc)'
  ytitle = 'Tangential Shear'
  erase & multiplot, [1, nfile]
  FOR i=0, nfile-1 DO BEGIN 

      rdzobjshear, files[i], str,/silent
      n=n_elements(str.meanr)
      IF useall THEN nuse = n ELSE nuse = n-1
      str = str[0:nuse-1]          ;last bin is usually bad

      str.shear = str.shear*Ssh

      fitpower, str.meanr/fac,str.shear,str.shearerr,shguess,$
        yfit,a,asigma,P,/silent

      ;; a[1] (the power) will be negative.   Convert to positive
      ;; since we fit model a[0]*r^-beta
      a[1] = -a[1]

      IF nfile EQ 1 THEN BEGIN 
          ploterror,str.meanr,str.shear,str.shearerr,psym=1, $
            xtitle=xtitle,ytitle=ytitle,yrange=yrange,_extra=extra
      ENDIF ELSE BEGIN 
          CASE i OF 
              nfile-1: ploterror,str.meanr,str.shear,str.shearerr,psym=1, $
                xtitle=xtitle,ytitle=ytitle,yrange=yrange
              0:       ploterror, str.meanr, str.shear,str.shearerr,psym=1, $
                yrange=yrange, _extra=extra
              ELSE:    ploterror, str.meanr,str.shear,str.shearerr,psym=1,$
                yrange=yrange
          ENDCASE 
      ENDELSE 

      ;; overplot the model fit
      oplot, str.meanr, yfit
      oplot,[0,10000],[0,0]

      ;; annotation
      mess = strarr(3)
      mess[0] = names[i]
      mess[1] = 'a = '+ntostr(a[0],nstring) + ' !M+ '+ntostr(asigma[0],nstring)
      mess[2] = 'b = '+ntostr(a[1],nstring) + ' !M+ '+ntostr(asigma[1],nstring)
      legend, mess, /right

      IF i NE nfile-1 THEN multiplot
  ENDFOR 
  multiplot, /reset

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Now do sig_crit*shear
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; find plot range
  FOR i=0, nfile-1 DO BEGIN
      rdzobjshear, files[i], str,/silent
      str.sigma = str.sigma*Ssh
      trange = prange(str.sigma, str.sigmaerr)
      IF i EQ 0 THEN yrange = trange ELSE BEGIN
          yrange[0] = min([trange[0], yrange[0]])
          yrange[1] = max([trange[1], yrange[1]])
      ENDELSE 
  ENDFOR 

  IF NOT nowait THEN wait,waittime
  ytitle='!S!7R!3!R!A-!N(!Ml!3r) - !S!7R!3!R!A-!N(r) (h M!M!Ln!N!3 pc!U-2!N)!X'
  ytitle2='gm cm!U-2'
  erase & multiplot, [1, nfile]
  FOR i=0, nfile-1 DO BEGIN 

      rdzobjshear, files[i], str,/silent
      n=n_elements(str.meanr)
      IF useall THEN nuse = n ELSE nuse = n-1
      str = str[0:nuse-1]          ;last bin is usually bad

      str.shear = str.shear*Ssh
      str.sigma = str.sigma*Ssh

      fitpower, str.meanr/fac,str.sigma,str.sigmaerr,siguess,$
        yfit,a,asigma,P,/silent

      ;; a[1] (the power) will be negative.   Convert to positive
      ;; since we fit model a[0]*r^-beta
      a[1] = -a[1]

      IF nfile EQ 1 THEN BEGIN
          plot, str.meanr, str.sigma, /nodata, ytick_get=v, $
            xtitle=xtitle, ytitle=ytitle, yrange=yrange, _extra=extra
      ENDIF ELSE BEGIN 
          CASE i OF
              nfile-1:  plot, str.meanr, str.sigma, /nodata, ytick_get=v, $
                xtitle=xtitle, ytitle=ytitle, yrange=yrange,ystyle=8
              0:        plot, str.meanr, str.sigma,/nodata, $
                ytick_get=v, yrange=yrange, _extra=extra,ystyle=8
              ELSE:     plot, str.meanr, str.sigma,/nodata, $
                ytick_get=v, yrange=yrange,ystyle=8
          ENDCASE 
      ENDELSE 
      ;; plot only axis first to get tick_values.

      ;; do other axis in different units.
      oploterror, str.meanr,str.sigma,str.sigmaerr,psym=1
      vstr=ntostr(v*sigconvert,5)
      IF i EQ nfile-1 THEN axis,yaxis=1,ytickname=vstr,ytitle=ytitle2 $
      ELSE axis,yaxis=1,ytickname=vstr

      ;; Overplot the model fit
      oplot, str.meanr, yfit
      oplot,[0,10000],[0,0]

      ;; annotation
      mess = strarr(3)
      mess[0] = names[i]
      mess[1] = 'a = '+ntostr(a[0],nstring) + ' !M+ '+ntostr(asigma[0],nstring)
      mess[2] = 'b = '+ntostr(a[1],nstring) + ' !M+ '+ntostr(asigma[1],nstring)
      legend, mess, /right

      IF i NE nfile-1 THEN multiplot
  ENDFOR 
  multiplot, /reset

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Now sigma(r) and mass(r)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF NOT nowait THEN wait,waittime

  erase & multiplot, [1, nfile]

  ;; First do Sigma then M(r)
  ytitle = ['!S!7R!3!R!A-!N(r)   (h M!M!Ln!N!3 pc!U-2!N)!X', $
            'M(r)  (10!U15!N 1/h M!M!Ln!N)!X' ]
  ls = [0]

  yrange=[-50,200]
  FOR i=0, nfile-1 DO BEGIN

      rdzobjshear, files[i], str,/silent
      n=n_elements(str.meanr)
      IF useall THEN nuse = n ELSE nuse = n-1
      str = str[0:nuse-1]          ;last bin is usually bad
      
      str.shear = str.shear*Ssh
      str.sigma = str.sigma*Ssh
      
      fitpower, str.meanr/fac,str.sigma,str.sigmaerr,siguess,$
        yfit,a,P,/silent
      a[1] = -a[1]
      
      betafac = (2. - a[1])/a[1]
      yfit = yfit*betafac       ;Convert to sigma
      yerr = str.sigmaerr*betafac
      yr=[yrange[0], yrange[1]]
                                ;yr = prange(yfit,yerr)*1.1
      xt='' 
      sigma=yfit & sigerr=yerr
       
      IF i EQ nfile-1 THEN BEGIN
          yt=ytitle[0]
          yt2 = ytitle2
      ENDIF ELSE BEGIN
          yt=''
          yt2 = ''
      ENDELSE 
      plot, str.meanr, yfit, /nodata, ytick_get=v, $
        xtitle=xt,ytitle=yt,yrange=yr, $
        _extra=extra, ystyle=8
      oploterror, str.meanr, yfit, yerr ;, linestyle=i
      
      ;; do other axis in different units for sigma
          
      vstr=ntostr(v*sigconvert,5)
      axis,yaxis=1,ytickname=vstr,ytitle=yt2
      oplot,[0,10000],[0,0]

      IF i NE nfile-1 THEN multiplot

      IF i GT 0 THEN ls = [ls,i]
      mess = names[i]
      legend, mess, /right
  ENDFOR 
  multiplot,/reset

  erase & multiplot, [1, nfile]

  ;; M(r)
  yr=[0,2]
  FOR i=0, nfile-1 DO BEGIN

      rdzobjshear, files[i], str,/silent
      n=n_elements(str.meanr)
      IF useall THEN nuse = n ELSE nuse = n-1

      str = str[0:n-2]          ;last bin is usually bad
      
      str.sigma = str.sigma*Ssh
      
      fitpower, str.meanr/fac,str.sigma,str.sigmaerr,siguess,$
        yfit,a,P,/silent
      a[1] = -a[1]
      
      betafac = 2.*!pi/a[1]*(str.meanr*1000.)^2/munit
      yfit = yfit*betafac       ;convert to Mass, also kpc->pc
      yerr = str.sigmaerr*betafac
      ;yr = prange(yfit,yerr)*1.1
      xt=xtitle 
      mass=yfit & merr=yerr

      IF i EQ nfile-1 THEN yt=ytitle[1] ELSE yt=''
      plot, str.meanr, yfit, /nodata, ytick_get=v, $
        xtitle=xt,ytitle=yt,yrange=yr
      oploterror, str.meanr, yfit, yerr ;, linestyle=i
      
      ;; do other axis in different units for sigma
      
      oplot,[0,10000],[0,0]
      
      IF i NE nfile-1 THEN multiplot

      IF i GT 0 THEN ls = [ls,i]
      mess = names[i]
      legend, mess, /right
  ENDFOR 

  multiplot,/reset
  !p.noerase = 0

  meanr = str.meanr
  print
  print,'This used MULTIPLOT, you may need to use ERASE'
  !p.multi=pold
  !x.margin = xmarg

  return 
END 
