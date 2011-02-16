
PRO stripe_range_plots, stripe, dops=dops

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: stripe_range, stripe, [/dops]'
      return
  ENDIF 

  !p.multi=[0,1,2]

  radec_search_dir = sdssidl_config('radec_search_dir')
  outdir=radec_search_dir+'Radecplots/'

  psfile = outdir+'stripe'+ntostr(stripe)+'_range.ps'
  pngfile = outdir+'stripe'+ntostr(stripe)+'_range.png'
  print
  print,'Doing stripe = '+ntostr(stripe)
  print

  get_good_lenstripe, stripe, run, rerun, stripes, strips

  IF run[0] EQ -1 THEN return

  ;; get rid of duplicates
  rmd = rem_dup(run)
  run=run[rmd]
  rerun = rerun[rmd]
  strips = strips[rmd]

  ;; sort by strip
  s=sort(strips)
  run = run[s]
  rerun = rerun[s]
  strips = strips[s]

  stripestr = ntostr(stripe)

  IF keyword_set(dops) THEN begin
      begplot, name=psfile, /color , /landscape
      simpctable
      !p.charsize=1
  ENDIF 
  colors=[!black, !red, !blue, !green, !magenta, !orange, !lightblue,$
          !grey50, !hotpink, !darkred, !seagreen, !purple, !salmon, !turquoise]
  nclr = n_elements(colors)
  nrun=n_elements(run)

  check=1
  FOR i=0L, nrun-1 DO BEGIN 
        
      runstr = run2string(run[i])

      FOR camcol=1, 6 DO BEGIN 
          file=radec_search_dir+$
            'lameta-range-'+runstr+'-'+ntostr(camcol)+'.fit'
          IF fexist(file) THEN BEGIN 
              lametatmp=mrdfits(file,1,/silent,status=status)
              IF status LT 0 THEN print,'What!'
              add_tag, lametatmp, 'run', run[i], lameta
              
              IF check THEN BEGIN 
                  lametarange = lameta
                  check=0
              ENDIF ELSE BEGIN
                  concat_structs, lametarange, lameta, tmp
                  lametarange=temporary(tmp)
              ENDELSE 
          ENDIF ELSE BEGIN
              IF camcol EQ 6 THEN run[i] = -1
          ENDELSE 
      ENDFOR 
  ENDFOR 

  maxeta = max(lametarange.etamax)
  mineta = min(lametarange.etamin)
  maxlam = max(lametarange.lambdamax)
  minlam = min(lametarange.lambdamin)

  IF minlam LT 0 THEN BEGIN
      minlam = 1.5*minlam
  ENDIF ELSE BEGIN
      minlam = 0.5*minlam
  ENDELSE 

  IF stripe GT 45 THEN xrange=[-300,200] ELSE xrange=[-150,100]

  plot,[0],xrange=xrange,yrange=[mineta,maxeta],/nodata, $
       xtitle='lambda', ytitle='eta', title='Stripe '+ntostr(stripe),$
       background=!white,color=!black

  clri=0
  line=0
  FOR i=0, nrun-1 DO BEGIN 

      IF run[i] NE -1 THEN BEGIN 
          print,stripestr+' ',strips[i]+' ',run[i]
          prunstr=stripestr+'-'+strips[i]+'-'+run2string(run[i])
          add_arrval, run[i], goodruns, /front
          add_arrval, strips[i], goodstrips, /front
          add_arrval, prunstr, comm, /front
          add_arrval, colors[clri], runclr, /front
          add_arrval, !black, textcolor, /front
          add_arrval, line, linestyle, /front
          add_arrval, !p.thick, thick, /front

          w=where(lametarange.run EQ run[i],nw)

          tminlam = min(lametarange[w].lambdamin, max=tmaxlam)
          add_arrval, tminlam, rminlam, /front
          add_arrval, tmaxlam, rmaxlam, /front

          FOR j=0L, nw-1 DO BEGIN 
              ind=w[j]
              IF (lametarange[ind].lambdamin LT -170) AND (lametarange[ind].lambdamax GT 170.) THEN BEGIN
                  cross=1
              ENDIF ELSE BEGIN 
                  cross=0
              ENDELSE 
              IF NOT cross THEN BEGIN
                  x = [lametarange[ind].lambdamin,$
                       lametarange[ind].lambdamax,$
                       lametarange[ind].lambdamax,$
                       lametarange[ind].lambdamin]
                  y = [lametarange[ind].etamin,$
                       lametarange[ind].etamin,$
                       lametarange[ind].etamax,$
                       lametarange[ind].etamax]

                  polyfill, x, y, color=colors[clri]
              ENDIF
          ENDFOR 
          clri=clri+1
          IF clri GE nclr THEN BEGIN 
              clri=0
              line=line+2
          ENDIF 
      ENDIF 
  ENDFOR 

  lchar = 0.75
  legend, reverse(comm), linestyle=linestyle, $
          colors=reverse(runclr), textcolor=textcolor,$
          thick=thick, /left,charsize=lchar,box=0
;  legend, comm, linestyle=linestyle, $
;          colors=runclr, textcolor=textcolor,$
;          thick=thick, /left,charsize=lchar,box=0

  ;;IF !d.name EQ 'X' THEN key=get_kbrd(1)
  
  ;; plot in not-overlapping way
  plot,[0],xrange=xrange,yrange=[0,1],/nodata, $
       xtitle='lambda', ytitle='run', title='Stripe '+ntostr(stripe),$
       background=!white,color=!black
  ngood = n_elements(goodruns)
  miny = 0.1
  maxy = 0.6
  yvals = arrscl(findgen(ngood), miny, maxy)

  myusersym,'fill_circle'
  FOR i=0L, ngood-1 DO BEGIN 

      w=where(lametarange.run EQ goodruns[i],nw)
      ysend = replicate(yvals[i], nw)
      oplot, lametarange[w].lambdamin, ysend, psym=8, symsize=0.3, $
             color=runclr[i]
;      oplot, [rminlam[i], rmaxlam[i]], [yvals[i], yvals[i]], $
;             color=runclr[i], thick=4

  ENDFOR 
  legend, reverse(comm), linestyle=linestyle, $
          colors=reverse(runclr), textcolor=textcolor,$
          thick=thick, /left,charsize=lchar,box=0
;  legend, comm, linestyle=linestyle, $
;          colors=runclr, textcolor=textcolor,$
;          thick=thick, /left,charsize=lchar,box=0


  IF keyword_set(dops) THEN BEGIN 
      endplot,/noprint
      pslandfix,psfile
      spawn,'pstoimg -antialias -density 75 -crop a -flip cw -out '+pngfile+' '+psfile
  ENDIF 

  !p.multi=0

  return
END 
