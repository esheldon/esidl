PRO chisq_conf, datax, data, covariance, modelx, model, param1, param2, $
                chisq_surf, $
                best1, best2, $
                low1, high1, $
                low2, high2, $
                errlow1=errlow1,errhigh1=errhigh1,$
                errlow2=errlow2,errhigh2=errhigh2, $
                $
                likelihood=likelihood, $
                like1=like1, like2=like2, $
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
                allow1 = allow1, allow2=allow2, $
                $
                getpath=getpath, info=info, xy=xy, $
                $
                nodisplay=nodisplay, $
                aspect=aspect, center=center, $
                xtitle=xtitle, ytitle=ytitle, $
                noplotmin=noplotmin,$
                xtick_get=xtick_get, ytick_get=ytick_get, $
                dolegend=dolegend, names=names, nkeep=nkeep, $
                project=project, $
                domean=domean, input_best1=input_best1, input_best2=input_best2, $
                noprompt=noprompt, psfile_front=psfile_front, $
                _extra=extra

  IF n_params() LT 7 THEN BEGIN
      print,'-Syntax: chisq_conf, datax, data, covariance, modelx, model, param1, param2, $'
      print,'             chisq_surf, $'
      print,'             best1, best2, $'
      print,'             low1, high1, $'
      print,'             low2, high2, $'
      print,'             errlow1=errlow1,errhigh1=errhigh1,$'
      print,'             errlow2=errlow2,errhigh2=errhigh2, $'
      print,'             $'
      print,'             likelihood=likelihood, $'
      print,'             like1=like1, like2=like2, $'
      print,'             wuse=wuse,$'
      print,'             $'
      print,'             diagonals=diagonals,$'
      print,'             $'
      print,'             yfit=yfit,$'
      print,'             yallow_low=yallow_low, yallow_high=yallow_high,$'
      print,'             minchisq=minchisq,$'
      print,'             degfree=degfree,$'
      print,'             $'
      print,'             chisq_diff=chisq_diff, $'
      print,'             $'
      print,'             allow1 = allow1, allow2=allow2, $'
      print,'             $'
      print,'             getpath=getpath, info=info, xy=xy, $'
      print,'             $'
      print,'             nodisplay=nodisplay, $'
      print,'             aspect=aspect, center=center, $'
      print,'             noplotmin=noplotmin,$'
      print,'             xtick_get=xtick_get, ytick_get=ytick_get, $'
      print,'             dolegend=dolegend, names=names, nkeep=nkeep, $'
      print,'             project=project, $'
      print,'             domean=domean, input_best1=input_best1, input_best2=input_best2, $'
      print,'             noprompt=noprompt, $'
      print,'             _extra=extra'
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
      names=['P1', 'P2']
      IF n_elements(xtitle) NE 0 THEN names[0] = xtitle
      IF n_elements(ytitle) NE 0 THEN names[1] = ytitle
  ENDIF 

  ;; projected values
  default = 1.e10
  p_low1 = replicate(default,nlevels)
  p_high1 = replicate(default,nlevels)
  p_low2 = replicate(default,nlevels)
  p_high2 = replicate(default,nlevels)

  s1=n_elements(param1)
  s2=n_elements(param2)
  
  sz_model = size(model)
  sz_modelx = n_elements(modelx)

  IF sz_model[0] NE 3 THEN BEGIN
      print,'bad1'
      return
  ENDIF 
  IF sz_model[1] NE s1 THEN BEGIN
      print,'bad2'
      return
  ENDIF 
  IF sz_model[2] NE s2 THEN BEGIN
      print,'bad3'
      return
  ENDIF 
  IF sz_model[3] NE sz_modelx THEN BEGIN
      print,'bad4'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; index for x,y
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  index = lindgen(long(s1)*long(s2))
  x = index MOD s1
  y = index/s1

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

  chisq_surf = dblarr( s1, s2 )

  IF NOT docov THEN BEGIN 
      ;; Input covariance is really just dataerr
      errsquared=covariance^2
      FOR iy=0L, s2-1 DO BEGIN ; loop over par2
          
          FOR ix=0L, s1-1 DO BEGIN ;loop over par1
              
              tmpmod = model[ix, iy, *]
          
              modfunc=interpol(tmpmod, modelx, datax)
                            
              ff= (data[wuse]-modfunc[wuse])^2/errsquared[wuse]
              chisq_surf[ix, iy] = total(ff)
              
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

      FOR iy=0L, s2-1 DO BEGIN ; loop over par2
          
          FOR ix=0L, s1-1 DO BEGIN ;loop over par1
              
              tmpmod = reform( model[ix, iy, *] )
          
              modfunc=interpol(tmpmod, modelx, datax)
              
              diff = data[wuse] - modfunc[wuse]

              chisq_surf[ix, iy] = transpose(diff)#(reform(cinv##diff))

          ENDFOR 
      ENDFOR 

  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; find min chi squared
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  minchisq = min(chisq_surf, /nan)
  
  chisq_diff = chisq_surf - minchisq

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; get pow/val of min
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  w=where(chisq_surf EQ minchisq, nw)

  IF n_elements(input_best1) NE 0 THEN best1=input_best1 ELSE best1 = (param1[x[w]])[0]
  IF n_elements(input_best2) NE 0 THEN best2=input_best2 ELSE best2 = (param2[y[w]])[0]

  pltmin=1
  IF keyword_set(noplotmin) THEN pltmin=0
  IF keyword_set(nodisplay) THEN pltmin=0

  print,'Min Chisq: ',minchisq,'/',ntostr(long(degfree)),' = ',ntostr(minchisq/degfree)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; display?
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF NOT keyword_set(nodisplay) THEN BEGIN 

      IF n_elements(psfile_front) NE 0 THEN BEGIN 
          begplot,name=psfile_front+'_surf.eps',/encap
      ENDIF 

      IF n_elements(xtitle) EQ 0 THEN xtit=names[0] ELSE xtit=xtitle
      IF n_elements(ytitle) EQ 0 THEN ytit=names[1] ELSE ytit=ytitle

      IF n_elements(aspect) EQ 0 THEN aspect = 1

      IF keyword_set(getpath) THEN BEGIN 
          acontour, aspect, chisq_diff, param1, param2, $
                    levels=!siglevels2,c_line=c_line, $
                    _extra=extra, xtick_get=xtick_get, ytick_get=ytick_get, $
                    path_xy=xy, path_info=info, /path_data_coords
      ENDIF 
      
      acontour, aspect, chisq_diff, param1, param2, $
                levels=!siglevels2,c_line=c_line, $
                _extra=extra, center=center, $
                xtick_get=xtick_get, ytick_get=ytick_get, $
                xtitle=xtit, ytitle=ytit
      
      IF pltmin THEN oplot,[best1],[best2],psym=7

  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; find confidence region
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; joint allowed par1,par2

  FOR i=0L, nlevels-1 DO BEGIN
          
      w = where(chisq_diff LE !siglevels2[i], nw)
      
      IF nw NE 0 THEN BEGIN 
          p_low1[i] = min( param1[ x[w] ] )
          p_high1[i] = max( param1[ x[w] ] )
          
          p_low2[i] = min( param2[ y[w] ] )
          p_high2[i] = max( param2[ y[w] ] )
      ENDIF 
      IF i EQ 0 THEN BEGIN ;; save 1-sigma region
          allow1 = param1[ x[w] ]
          allow2 = param2[ y[w] ]
      ENDIF 
  ENDFOR 

  IF keyword_set(project) THEN BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Do projection to get error regions
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      low1 = p_low1
      high1 = p_high1
      low2 = p_low2
      high2 = p_high2

      range2error, low1[0], best1, high1[0], errlow10, errhigh10
      range2error, low1[1], best1, high1[1], errlow11, errhigh11
      range2error, low1[2], best1, high1[2], errlow12, errhigh12
      errlow1=[errlow10,errlow11,errlow12]
      errhigh1=[errhigh10,errhigh11,errhigh12]

      range2error, low2[0], best2, high2[0], errlow20, errhigh20
      range2error, low2[1], best2, high2[1], errlow21, errhigh21
      range2error, low2[2], best2, high2[2], errlow22, errhigh22
      errlow2=[errlow20,errlow21,errlow22]
      errhigh2=[errhigh20,errhigh21,errhigh22]
      
  ENDIF ELSE BEGIN 
      
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Marginalize to get error regions
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ;; create likelihood surface
      likelihood = exp( - (0.5d*chisq_diff < 20d) )

      marginalize_like2d, param1, param2, likelihood, like1, like2

      ;; find 1,2,3-sigma regions

      get_like_conf, param1, like1, $
                     best1, low1, high1, errlow1, errhigh1, /silent, domean=domean, input_vmax=input_best1

      get_like_conf, param2, like2, $
                     best2, low2, high2, errlow2, errhigh2, /silent, domean=domean, input_vmax=input_best2

  ENDELSE 

  ;; Y values for best fit and 1-sigma allowed region
  IF keyword_set(domean) THEN BEGIN 
      p1diff = abs(param1-best1)
      p2diff = abs(param2-best2)
      twx = where(p1diff EQ min(p1diff))
      twy = where(p2diff EQ min(p2diff))
  ENDIF ELSE BEGIN 
      twx = where(param1 EQ best1)
      twy = where(param2 EQ best2)
  ENDELSE 
  tmpmod = reform(model[ twx, twy, *])

  yfit = interpol(tmpmod, modelx, datax[wuse])

  yallow_low = fltarr(ndata)
  yallow_high = yallow_low

  ;; Allowed 1-sigma region of y values
  ;; don't know how to do this yet
;  FOR i=0L, ndata-1 DO BEGIN 
;      mody = allow2*(datax[wuse[i]])^allow1

;      yallow_low[wuse[i]] = min(mody)
;      yallow_high[wuse[i]] = max(mody)

;  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; some more plotting stuff
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(dolegend) AND NOT keyword_set(nodisplay) THEN BEGIN 

      mean_error_legend,names,[best1,best2],[errlow1[0],errlow2[0]],$
                        [errhigh1[0],errhigh2[0]],nkeep=nkeep,/right,$
                        charsize=1, box=0

      mess2 = !csym.chi+'!U2!N/'+!csym.nu+' = '+ntostr(minchisq,4,/round)+'/'+ntostr(degfree)+$
        ' = '+ntostr(minchisq/degfree,4,/round)
      legend, mess2,charsize=1, box=0


;      mean_error_legend,names,[best1,best2],[errlow1[0],errlow2[0]],$
;                        [errhigh1[0],errlow2[0]],nkeep=nkeep,/right

      IF n_elements(psfile_front) NE 0 THEN endplot
  ENDIF 

  IF NOT keyword_set(project) AND NOT keyword_set(nodisplay) THEN BEGIN 

      IF n_elements(psfile_front) NE 0 THEN BEGIN 
          begplot,name=psfile_front+'_marg.eps',/encap, /color
      ENDIF 

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
      aplot, 1, param1, like1/max(like1),$
             yrange=[0,1.1], xstyle=2, ystyle=1, $
             charsize=1, center=center, xtitle=xtit, ytitle=ytit

      oplot, [best1, best1], [0, 5], color=!red
      oplot, [low1[0], low1[0]], [0,5], line=2, color=errcolor
      oplot, [high1[0], high1[0]], [0,5], line=2, color=errcolor

      IF keyword_set(dolegend) THEN BEGIN 

          IF n_elements(nkeep) NE 0 THEN BEGIN 
              mean_error_legend,names[0],best1,errlow1[0],errhigh1[0],$
                                nkeep=nkeep[0],/right,$
                                charsize=1, box=0
          ENDIF ELSE BEGIN 
              mean_error_legend,names[0],best1,errlow1[0],errhigh1[0],$
                                /right,$
                                charsize=1, box=0
          ENDELSE 
      ENDIF 

      IF n_elements(ytitle) EQ 0 THEN xtit=names[1] ELSE xtit=ytitle
      aplot, 1, param2, like2/max(like2), $
             yrange=[0,1.1], xstyle=2, ystyle=1, $
             charsize=1, center=center, xtit=xtit,ytitle=ytit

      oplot, [best2, best2], [0, 5], color=!red
      oplot, [low2[0], low2[0]], [0,5], line=2, color=errcolor
      oplot, [high2[0], high2[0]], [0,5], line=2, color=errcolor

      IF keyword_set(dolegend) THEN BEGIN 

          IF n_elements(nkeep) NE 0 THEN BEGIN 
              mean_error_legend,names[1],best2,errlow2[0],errhigh2[0],$
                                nkeep=nkeep[1],/right,$
                                charsize=1, box=0
          ENDIF ELSE BEGIN 
              mean_error_legend,names[1],best2,errlow2[0],errhigh2[0],$
                                /right,$
                                charsize=1, box=0
          ENDELSE 
      ENDIF 

      !p.multi=pold

      IF n_elements(psfile_front) NE 0 THEN endplot

  ENDIF 




END 
