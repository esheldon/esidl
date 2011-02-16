
PRO pow_chisq_conf, datax, data, covariance, powvals, normvals, $
                    chisq_surf, $
                    bestp, bestn, $
                    powlow, powhigh, $
                    normlow, normhigh, $
                    perrlow=perrlow,perrhigh=perrhigh,$
                    nerrlow=nerrlow,nerrhigh=nerrhigh, $
                    $
                    likelihood=likelihood, $
                    pow_like=pow_like, norm_like=norm_like, $
                    $
                    wuse=wuse,$
                    $
                    diagonals=diagonals,$
                    $
                    yfit=yfit,$
                    yallow_low=yallow_low, yallow_high=yallow_high,$
                    minchisq=minchisq,$
                    degfree=degfree,$
                    $
                    chisq_diff=chisq_diff, $
                    $
                    powallow = powallow, normallow=normallow, $
                    $
                    getpath=getpath, info=info, xy=xy, $
                    $
                    nodisplay=nodisplay, $
                    aspect=aspect, center=center, $
                    xtitle=xtitle, ytitle=ytitle, $
                    noplotmin=noplotmin,$
                    xtick_get=xtick_get, ytick_get=ytick_get, $
                    dolegend=dolegend, names=names, nkeep=nkeep, $
                    project=project, domean=domean, $
                    noprompt=noprompt, $
                    _extra=extra

  IF n_params() LT 5 THEN BEGIN
      print,'-Syntax: pow_chisq_conf, datax, data, covariance, powvals, normvals, $'
      print,'       chisq_surf, $'
      print,'       bestp, bestn, $'
      print,'       powlow, powhigh, $'
      print,'       normlow, normhigh, $'
      print,'       perrlow=perrlow,perrhigh=perrhigh,$'
      print,'       nerrlow=nerrlow,nerrhigh=nerrhigh, $'
      print,'       $'
      print,'       likelihood=likelihood, $'
      print,'       pow_like=pow_like, norm_like=norm_like, $'
      print,'       $'
      print,'       wuse=wuse,$'
      print,'       $'
      print,'       diagonals=diagonals,$'
      print,'       $'
      print,'       yfit=yfit,$'
      print,'       yallow_low=yallow_low, yallow_high=yallow_high,$'
      print,'       minchisq=minchisq,$'
      print,'       degfree=degfree,$'
      print,'       $'
      print,'       chisq_diff=chisq_diff, $'
      print,'       $'
      print,'       powallow = powallow, normallow=normallow, $'
      print,'       $'
      print,'       getpath=getpath, info=info, xy=xy, $'
      print,'       $'
      print,'       nodisplay=nodisplay, $'
      print,'       aspect=aspect, center=center, $'
      print,'       noplotmin=noplotmin,$'
      print,'       xtick_get=xtick_get, ytick_get=ytick_get, $'
      print,'       dolegend=dolegend, names=names, nkeep=nkeep, $'
      print,'       project=project, domean=domean, $'
      print,'       noprompt=noprompt, $'
      print,'       _extra=extra'
      print
      return
  ENDIF 
      
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Define conf. levels
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nlevels = 3                   ;68, 95.5, 99.7%
  c_line=[0,0,0]

  errcode=['Successful completion',$
           'Singular Covariance Matrix', $
           'Small pivot element in Covariance Matrix: Accuracy Compromised']

  IF n_elements(names) EQ 0 THEN BEGIN 
      names=['Index', 'Norm']
      IF n_elements(xtitle) NE 0 THEN names[0] = xtitle
      IF n_elements(ytitle) NE 0 THEN names[1] = ytitle
  ENDIF 

  npow = n_elements(powvals)
  nnorm = n_elements(normvals)

  ;; projected values
  default = 1.e10
  p_powlow = replicate(default,nlevels)
  p_powhigh = replicate(default,nlevels)
  p_normlow = replicate(default,nlevels)
  p_normhigh = replicate(default,nlevels)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; index for x,y
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  index = lindgen(long(npow)*long(nnorm))
  x = index MOD npow
  y = index/npow

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; check form of covariance. Is it
  ;; cov matrix or just errors?
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  ndata = n_elements(data)
  cs = size(covariance)
  IF (cs[0] EQ 2) THEN BEGIN ;; 2-D?
      IF cs[1] NE ndata THEN BEGIN
          message,'Covariance matrix must be NdataXNdata',/inf
          return
      ENDIF 
      docov = 1
  ENDIF ELSE BEGIN ;; its 1-D
      docov = 0
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; which data points do we use?
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(wuse) EQ 0 THEN BEGIN
      ndata = n_elements(data)
      wuse = lindgen(ndata)
  ENDIF ELSE BEGIN 
      ndata = n_elements(wuse)
  ENDELSE 
  degfree = ndata-2

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; calculate chi squared surface
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  chisq_surf = dblarr( npow, nnorm )
  modfunc = data

  IF NOT docov THEN BEGIN 
      ;; Input covariance is really just dataerr
      errsquared=covariance^2
      FOR in=0L, nnorm-1 DO BEGIN ; loop over norms
          norm = normvals[in]
          FOR ip=0L, npow-1 DO BEGIN ;loop over powers
              pow = powvals[ip]
              
              modfunc[wuse] = norm*(datax[wuse])^pow
              
              ff= (data[wuse]-modfunc[wuse])^2/errsquared[wuse]
              chisq_surf[ip, in] = total(ff)
              
          ENDFOR 
      ENDFOR 

  ENDIF ELSE BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; copy in block corresponding to wuse indices
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      tmpcov = replicate(covariance[0], ndata, ndata)
      FOR i=0L, ndata-1 DO BEGIN 
          FOR j=0L, ndata-1 DO BEGIN 
              tmpcov[i,j] = covariance[wuse[i], wuse[j]]
          ENDFOR 
      ENDFOR 
        
      cinv = invert(tmpcov, stat, /double)
      IF stat NE 0 THEN message,errcode[stat]

      FOR in=0L, nnorm-1 DO BEGIN ; loop over norms
          norm = normvals[in]
          FOR ip=0L, npow-1 DO BEGIN ;loop over powers
              pow = powvals[ip]
              
              modfunc[wuse] = norm*(datax[wuse])^pow
              diff = data[wuse] - modfunc[wuse]

              chisq_surf[ip, in] = transpose(diff)#(reform(cinv##diff))

          ENDFOR 
      ENDFOR 

  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; find min chi squared
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  minchisq = min(chisq_surf, /nan)
  
  if minchisq lt 0.0 then print,'WARNING: minchisq is < 0'
  chisq_diff = chisq_surf - minchisq

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; get pow/val of min
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  w=where(chisq_surf EQ minchisq, nw)

  bestp = (powvals[x[w]])[0]
  bestn = (normvals[y[w]])[0]
  pltmin=1
  IF keyword_set(noplotmin) THEN pltmin=0
  IF keyword_set(nodisplay) THEN pltmin=0

  print,'Min Chisq: ',minchisq,'/',ntostr(long(degfree)),' = ',ntostr(minchisq/degfree)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; display?
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF NOT keyword_set(nodisplay) THEN BEGIN 

      IF n_elements(xtitle) EQ 0 THEN xtit=names[0] ELSE xtit=xtitle
      IF n_elements(ytitle) EQ 0 THEN ytit=names[1] ELSE ytit=ytitle

      IF n_elements(aspect) EQ 0 THEN aspect = 1

      IF keyword_set(getpath) THEN BEGIN 
          acontour, aspect, chisq_diff, powvals, normvals, $
                    levels=!siglevels2,c_line=c_line, $
                    _extra=extra, xtick_get=xtick_get, $
                    ytick_get=ytick_get, path_xy=xy, path_info=info, $
                    /path_data_coords
      ENDIF 
      
      acontour, aspect, chisq_diff, powvals, normvals, $
                levels=!siglevels2,c_line=c_line, $
                _extra=extra, center=center, xtick_get=xtick_get, $
                ytick_get=ytick_get, $
                xtitle=xtit, ytitle=ytit
      
      IF pltmin THEN oplot,[bestp],[bestn],psym=7

  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; find confidence region
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; joint allowed pow/norm

  FOR i=0L, nlevels-1 DO BEGIN
          
      w = where(chisq_diff LE !siglevels2[i], nw)
      
      IF nw NE 0 THEN BEGIN 
          p_powlow[i] = min( powvals[ x[w] ] )
          p_powhigh[i] = max( powvals[ x[w] ] )
          
          p_normlow[i] = min( normvals[ y[w] ] )
          p_normhigh[i] = max( normvals[ y[w] ] )
      ENDIF 
      IF i EQ 0 THEN BEGIN ;; save 1-sigma region
          powallow = powvals[ x[w] ]
          normallow = normvals[ y[w] ]
      ENDIF 
  ENDFOR 

  IF keyword_set(project) THEN BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Do projection to get error regions
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      powlow = p_powlow
      powhigh = p_powhigh
      normlow = p_normlow
      normhigh = p_normhigh

      range2error, powlow[0], bestp, powhigh[0], perrlow0, perrhigh0
      range2error, powlow[1], bestp, powhigh[1], perrlow1, perrhigh1
      range2error, powlow[2], bestp, powhigh[2], perrlow2, perrhigh2
      perrlow=[perrlow0,perrlow1,perrlow2]
      perrhigh=[perrhigh0,perrhigh1,perrhigh2]

      range2error, normlow[0], bestn, normhigh[0], nerrlow0, nerrhigh0
      range2error, normlow[1], bestn, normhigh[1], nerrlow1, nerrhigh1
      range2error, normlow[2], bestn, normhigh[2], nerrlow2, nerrhigh2
      nerrlow=[nerrlow0,nerrlow1,nerrlow2]
      nerrhigh=[nerrhigh0,nerrhigh1,nerrhigh2]
      
  ENDIF ELSE BEGIN 
      
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Marginalize to get error regions
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ;; create likelihood surface
      likelihood = exp( - (0.5d*chisq_diff < 20d) )

      marginalize_like2d, powvals, normvals, likelihood, pow_like, norm_like

      ;; find 1,2,3-sigma regions

      get_like_conf, normvals, norm_like, $
                     bestn, normlow, normhigh, nerrlow, nerrhigh, /silent, domean=domean
      get_like_conf, powvals, pow_like, $
                     bestp, powlow, powhigh, perrlow, perrhigh, /silent, domean=domean

; This was a test, where I generated values from the 2-d likelihood 
; surface.  The results agreed with the marginalized likelihood.

;      nrand = 100000L
;      genrand2d, likelihood, powvals, normvals, nrand, randpow, randnorm

;      print,mean(randpow),mean(randnorm)
;      print,bestp,bestn

;      !p.multi = [0,0,2]

;      plothist, randpow, bin=0.005, /norm
;      oplot, powvals, pow_like

;      plothist, randnorm, bin=0.01, /norm
;      oplot, normvals, norm_like

;      !p.multi=0

  ENDELSE 

  ;; Y values for best fit and 1-sigma allowed region
  yfit = fltarr(n_elements(datax))
  yallow_low = yfit
  yallow_high = yfit

  ;; best fit
  yfit[wuse] = bestn*datax[wuse]^bestp

  ;; Allowed 1-sigma region of y values
  yallow_low = yfit
  yallow_high = yfit

  FOR i=0L, ndata-1 DO BEGIN 
      mody = normallow*(datax[wuse[i]])^powallow

      yallow_low[wuse[i]] = min(mody)
      yallow_high[wuse[i]] = max(mody)

  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; some more plotting stuff
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(dolegend) AND NOT keyword_set(nodisplay) THEN BEGIN 
      mean_error_legend,names,[bestp,bestn],[perrlow[0],nerrlow[0]],$
                        [perrhigh[0],nerrhigh[0]],nkeep=nkeep,/right,$
                        charsize=1, box=0


      mess2 = !csym.chi+'!U2!N/'+!csym.nu+' = '+ntostr(minchisq,4,/round)+'/'+ntostr(degfree)+$
        ' = '+ntostr(minchisq/degfree,4,/round)
      legend, mess2,charsize=1, box=0

  ENDIF 

  IF NOT keyword_set(project) AND NOT keyword_set(nodisplay) THEN BEGIN 

      IF !d.name EQ 'X' THEN BEGIN
          errcolor = !green

          IF NOT keyword_set(noprompt) THEN BEGIN 
              print
              print,'Hit a key'
              key=get_kbrd(1)
          ENDIF 
      ENDIF ELSE BEGIN 
          errcolor = !blue
      ENDELSE 

      pold = !p.multi
      !p.multi=[0,0,2]

      ytit='L/L!Dmax!N'

      IF n_elements(xtitle) EQ 0 THEN xtit=names[0] ELSE xtit=xtitle
      aplot, 1, powvals, pow_like/max(pow_like), xtitle=xtit, ytitle=ytit,$
             yrange=[0,1.1], xstyle=2, ystyle=1, $
             charsize=1, center=center

      oplot, [bestp, bestp], [0, 5], color=!red
      oplot, [powlow[0], powlow[0]], [0,5], line=2, color=errcolor
      oplot, [powhigh[0], powhigh[0]], [0,5], line=2, color=errcolor

      IF keyword_set(dolegend) THEN BEGIN 
          IF n_elements(nkeep) NE 0 THEN BEGIN 
              mean_error_legend,names[0],bestp,perrlow[0], perrhigh[0],$
                                nkeep=nkeep[0],/right,$
                                charsize=1, box=0
          ENDIF ELSE BEGIN
              mean_error_legend,names[0],bestp,perrlow[0], perrhigh[0],$
                                /right,$
                                charsize=1, box=0
          ENDELSE 
                    
      ENDIF 

      IF n_elements(ytitle) EQ 0 THEN xtit=names[1] ELSE xtit=ytitle
      aplot, 1, normvals, norm_like/max(norm_like), $
             xtitle=xtit, ytitle=ytit,yrange=[0,1.1], xstyle=2, ystyle=1, $
             charsize=1, center=center

      oplot, [bestn, bestn], [0, 5], color=!red
      oplot, [normlow[0], normlow[0]], [0,5], line=2, color=errcolor
      oplot, [normhigh[0], normhigh[0]], [0,5], line=2, color=errcolor

      IF keyword_set(dolegend) THEN BEGIN 
          IF n_elements(nkeep) NE 0 THEN BEGIN 
              mean_error_legend,names[1],bestn,nerrlow[0],nerrhigh[0],$
                                nkeep=nkeep[1],/right,$
                                charsize=1, box=0
          ENDIF ELSE BEGIN 
              mean_error_legend,names[1],bestn,nerrlow[0],nerrhigh[0],$
                                /right,$
                                charsize=1, box=0
          ENDELSE 
          
      ENDIF 


      !p.multi=pold

  ENDIF 


END 

