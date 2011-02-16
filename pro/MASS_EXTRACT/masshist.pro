PRO masshist, kcat, ncat, bin, khist, khist_err, nhist, nhist_err, xhist, name=name, med=med, _extra=extra, maxx=maxx, minx=minx, maxy=maxy, miny=miny, step=step, nstep=nstep, slength=slength

  IF N_params() LT 3 THEN BEGIN 
     print,'-Syntax: masshist, kcat, ncat, bin, khist, khist_err, nhist, nhist_err, xhist, name=name,_extra=extra'
     print,''
     print,'Use doc_library,"masshist"  for more help.'  
     return
  ENDIF 

  IF NOT keyword_set(med) THEN med=0

  IF n_elements(maxx) NE 0 THEN BEGIN
      ;; Build title

      binarea = ntostr( (maxx - minx)^2*step, 6)
      nfields = ntostr(long(nstep))
      sstr = ntostr(long(slength))
      title='Field Size: '+binarea+' deg^2    Nfields: '+nfields+'   S='+sstr
  ENDIF ELSE title=''

  yrange=[0., .3]

  min = -5.0
  max = 5.0

  IF n_elements(psym) EQ 0 THEN psym=0

  xtitle = 'S/N'
  ytitle = 'pdf'
  items = ['Kappa', 'Noise']
  linestyle = [0,2]

  wild = kcat[ rem_dup(kcat.wild) ].wild

  nwild = n_elements(wild)
  IF nwild EQ 1 THEN BEGIN 

      nhist, ncat.s2n, bin, linestyle = linestyle[1], $
        xtitle=xtitle, ytitle=ytitle, psym=psym, _extra=extra
      nhist, kcat.s2n, bin, linestyle = linestyle[0], /overplot, psym=psym
      legend, items, linestyle=linestyle
  ENDIF ELSE BEGIN 
  
      khist0 = histogram(kcat.s2n, binsize=bin, min=min, max=max)
      nhist0 = histogram(ncat.s2n, binsize=bin, min=min, max=max)

      nbin = n_elements(khist0)

      xhist = arrscl(findgen(nbin), min, max)

      khist_tmp = fltarr(nwild, nbin)
      nhist_tmp = khist_tmp
      khist = fltarr(nbin) & nhist = khist
      khist_err = khist & nhist_err = khist
      
      ;; Estimate the error bars.
      FOR i=0, nwild-1 DO BEGIN 
          
          w = where(kcat.wild EQ wild[i])
          
          k_tmp = histogram(kcat[w].s2n, binsize=bin, min=min, max=max)
          n_tmp = histogram(ncat[w].s2n, binsize=bin, min=min, max=max)
          k_tmp = normalize(float(k_tmp), xhist)
          n_tmp = normalize(float(n_tmp), xhist)
      
          khist_tmp[i, *] = k_tmp
          nhist_tmp[i, *] = n_tmp
      ENDFOR 

      FOR i=0, nbin-1 DO BEGIN
          
          mtmp = moment( khist_tmp[*, i] )

          IF NOT med THEN BEGIN 
              khist[i] = mtmp[0]
              khist_err[i] = sqrt( !pi/2. * mtmp[1]/nwild )
          ENDIF ELSE BEGIN 
              khist[i] = median(khist_tmp[*, i])
              khist_err[i] = sqrt(mtmp[1]/nwild)
          ENDELSE 

          mtmp = moment( nhist_tmp[*, i] )
          IF NOT med THEN BEGIN 
              nhist[i] = mtmp[0]
              nhist_err[i] = sqrt( !pi/2. * mtmp[1]/nwild )
          ENDIF ELSE BEGIN 
              nhist[i] = median(nhist_tmp[*, i])
              nhist_err[i] = sqrt(mtmp[1]/nwild)
          ENDELSE 
      ENDFOR 

      ;; Individual plots
;      ploterr, xhist, nhist, nhist_err, psym=psym, $
;               linestyle=linestyle[1], xtitle=xtitle, ytitle=ytitle, $
;               title = 'Noise',ystyle=1, yrange=[0., .3]
;      key=get_kbrd(1)
;      ploterr, xhist, khist, khist_err, psym=psym, $
;               linestyle=linestyle[1], xtitle=xtitle, ytitle=ytitle, $
;               title='Kappa',ystyle=1, yrange=[0., .3]
;      key=get_kbrd(1)

      ;; Over plot without error bars
;      plot, xhist, nhist, psym=psym, $
;               linestyle=linestyle[1], xtitle=xtitle, ytitle=ytitle, $
;               ystyle=1, yrange=[0., .3], _extra=extra
;      oplot, xhist, khist, psym=psym, $
;               linestyle=linestyle[0]
      ;; Over plot with error bars
      ploterr, xhist, nhist, nhist_err, psym=psym, $
               linestyle=linestyle[1], xtitle=xtitle, ytitle=ytitle, $
               ystyle=1, yrange=[0., .3], title=title, _extra=extra
      oploterr, xhist, khist, khist_err, psym=psym, $
               linestyle=linestyle[0]
      legend, items, linestyle=linestyle

      IF n_elements(name) NE 0 THEN BEGIN 
          openw, lun, name+'.txt', /get_lun
          message = '    s/n bin     noise      noise error      kappa    kappa error    '

          printf, lun, message
          !textunit = lun
          forprint,text=5, /silent, $
                   xhist, nhist, nhist_err, khist, khist_err
          close, lun
          free_lun, lun
      ENDIF 

  ENDELSE     
  return
END 
