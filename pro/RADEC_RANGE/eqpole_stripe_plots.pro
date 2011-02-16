PRO eqpole_stripe_plots, dops=dops, dopng=dopng

  radec_search_dir = sdssidl_config('radec_search_dir')

  outdir = radec_search_dir + 'Radecplots/'
  IF keyword_set(dops) THEN BEGIN 
      outfile = outdir + 'stripe_range.ps'
      begplot,name=outfile,/color,xsize=8,ysize=8
      position = [250,20000]
      device=1
  ENDIF ELSE IF keyword_set(dopng) THEN BEGIN 
      device_old = !d.name
      outfile = outdir + 'stripe_range.png'
      setupplot,'Z'
      device, set_resolution=[640,640]
      position = [25,625]
      device=1
  ENDIF ELSE BEGIN 
      window,xsize=640,ysize=640
      position = [25,625]
      device=1
  ENDELSE 

  run_status = sdss_runstatus()
  w=where(run_status.stripe NE -1,nw)
  IF nw EQ 0 THEN message,'No good stripes'

  w2=rem_dup(run_status[w].stripe)
  w=w[w2]

  allstripes = run_status[w].stripe

  nstripe=n_elements(allstripes)

  simpctable, rmap, gmap, bmap

  IF !d.name EQ 'PS' THEN BEGIN 
      colors=[!black, !red, !blue, !green, !magenta, !orange, !lightblue,$
              !grey50, !hotpink, !darkred, !seagreen, !purple, !salmon, $
              !turquoise]
      !p.charsize=1      
  ENDIF ELSE BEGIN 
      colors=[!white, !red, !blue, !green, !magenta, !orange, !lightblue,$
              !grey50, !hotpink, !darkred, !seagreen, !purple, !salmon, $
              !turquoise]
  ENDELSE 

  label=1
  gcolor=!p.color

  nclr = n_elements(colors)

  pold = !p.color
  !p.color=gcolor
  eqpole_grid,/new,color=gcolor,label=label

  FOR is=0L, nstripe-1 DO BEGIN 
      
      stripe = allstripes[is]
      
      print
      print,'Stripe = '+ntostr(stripe)
      get_good_lenstripe, stripe, run, rerun, stripes, strips

      IF run[0] NE -1 THEN BEGIN 

          ;; get rid of duplicates
          rmd = rem_dup(run)
          run=run[rmd]
          rerun = rerun[rmd]
          strips = strips[rmd]
          
          nrun = n_elements(run)

          ;; sort by strip
          s=sort(strips)
          run = run[s]
          rerun = rerun[s]
          strips = strips[s]

          stripestr = ntostr(stripe)
          add_arrval, stripestr, sstripestr
          scolor = is MOD nclr
          add_arrval, colors[scolor], scolors

          FOR i=0L, nrun-1 DO BEGIN 

              runstr = run2string(run[i])

              FOR camcol=1, 6 DO BEGIN 
                  file=radec_search_dir+$
                    'lameta-range-'+runstr+'-'+ntostr(camcol)+'.fit'
                  IF fexist(file) THEN BEGIN 
                      lameta = mrdfits(file,1,/silent,status=status)
                      IF status LT 0 THEN message,'What!'

                      eta = [lameta.etamin, lameta.etamax]
                      lambda = [lameta.lambdamin, lameta.lambdamax]

                      survey2eq, lambda, eta, ra, dec

                      eqpole, ra, dec, x,y  

                      plots,x,y,psym=3,color=colors[scolor]

                  ENDIF 
              ENDFOR 
          ENDFOR 
      ENDIF 

  ENDFOR 

  ;;eqpole_grid, label=label

  !p.color=pold

  ngs=n_elements(sstripestr)
;  legend, sstripestr, line=replicate(0, ngs), $
;          thick=replicate(!p.thick, ngs),$
;          colors=scolors,/clear
  legend, reverse(sstripestr), line=replicate(0, ngs), $
          thick=replicate(!p.thick, ngs),$
          colors=reverse(scolors),/clear, position=position, device=device

  IF keyword_set(dops) THEN BEGIN 
      endplot
  ENDIF ELSE IF keyword_set(dopng) THEN BEGIN 
      print,'Making png file: ',outfile
      IF float(!version.release) LT 5.4 THEN BEGIN 
          tmp = tvrd()
          tmp = rotate( rotate(tmp,1), 4)
          write_png, outfile, tmp, rmap, gmap, bmap
      ENDIF ELSE BEGIN 
          write_png, outfile, tvrd(), rmap, gmap, bmap
      ENDELSE 
      setupplot,device_old
  ENDIF 

END 
