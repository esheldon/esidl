PRO des_plot_psf_sampling, setnum, dops=dops

  ;; Files

  IF sdssidl_config('hostname') EQ 'treebeard' THEN BEGIN 
      inDir = '/data0/des_sim/sampling_errors/'
  ENDIF ELSE BEGIN 
      inDir = '/net/cheops1/data0/esheldon/des_sim/sampling_errors/'
  ENDELSE 

  front = 'sampling'+strn(setnum, length=2, padchar='0')

  files = file_search(inDir+front+'_psf*.st')

  t = read_idlstruct(files[0],/silent)
  nellip = n_elements(t.e1in)
  e1vals = t.e1in
  e2vals = t.e2in


  print
  print,'Sampling set: '+front
  print,'Number of ellipticities: ',nellip


  ;; For /all, we only need to input the seeing set
  seeing_set = 1
  des_psf_sampling_getsets, seeing_set, dummy2, seeing, apix, /all

  nseeing = n_elements(seeing)
  napix   = n_elements(apix)

  IF nseeing GT 1 THEN BEGIN 
      message,$
        'This program no longer supports more than one '+$
        'value of seeing per setnum'
  ENDIF 

  e1_slopes = dblarr(napix)
  e1_offsets = e1_slopes
  e2_slopes = e1_slopes
  e2_offsets = e1_slopes

  e1_sdev = fltarr(napix, nellip)
  e2_sdev = e1_sdev

  sz_in = fltarr(napix, nellip)
  sz_meas = fltarr(napix, nellip)

  fracuse = e1_sdev

  ;; Units for (e1meas-e1in)
  ediffunit = 1.e-3

  !p.multi=[0,2,4]
  !p.charsize = 2
  e1_xtitle = 'e!S!D1!N!R!UInput!N'
  e2_xtitle = 'e!S!D2!N!R!UInput!N'
  ytitle = 'difference/10!U-3!N'
  yrange = [-50, 50]
  irun = 0

  seeing_string = ntostr(seeing, 4, /round)
  FOR ip=0L, napix-1 DO BEGIN 

      arcperpix_string = ntostr(apix[ip], 5, /round)
          
      file = $
        indir + $
        front+'_psf'+seeing_string+'_pixel'+arcperpix_string+'.st'
      print,file

      pixscale_over_seeing = apix[ip]/seeing
      title = $
        'pixscale/seeing = '+ntostr(pixscale_over_seeing,5,/round)
      
      IF fexist(file) THEN BEGIN 
          
          t=read_idlstruct(file,/silent)
          
          IF n_elements(t.e1in) NE nellip THEN BEGIN 
              message,'n_elements(e1in) mismatch'
          ENDIF 
          
          s=sort(t.e1in)
          

          e1in = t[s].e1in
          e1diff    = (t[s].e1meas-e1in)/ediffunit
          e1differr = t[s].e1sdev/sqrt( t[s].npsf_use )/ediffunit


          s=sort(t.e2in)
          e2in = t[s].e2in
          e2diff    = (t[s].e2meas-e2in)/ediffunit
          e2differr = t[s].e2sdev/sqrt( t[s].npsf_use )/ediffunit


          IF nellip GT 1 THEN BEGIN 

              fitlin, e1in, e1diff, e1differr, e1a, e1erra, e1b, e1errb, $
                /silent
              fitlin, e2in, e2diff, e2differr, e2a, e2erra, e2b, e2errb, $
                /silent

              e1_offsets[ip] = e1a*ediffunit
              e1_slopes[ip] = e1b*ediffunit
              e2_offsets[ip] = e2a*ediffunit
              e2_slopes[ip] = e2b*ediffunit

              ploterror, e1in, e1diff, e1differr, $
                psym=8, $
                title=title, xtitle=e1_xtitle, ytitle=ytitle, $
                xstyle=1+2
              oplot, e1in, e1a+e1b*e1in

              ploterror, e2in, e2diff, e2differr, $
                psym=8, $
                title=title, xtitle=e2_xtitle, ytitle=ytitle, $
                xstyle=1+2
              oplot, e2in, e2a+e2b*e2in
          ENDIF 
          
          FOR ie=0L,nellip-1 DO BEGIN 
              e1_sdev[ip,ie] = t[ie].e1sdev
              e2_sdev[ip,ie] = t[ie].e2sdev

              findabtheta, t[ie].e1in, t[ie].e2in, aratio, theta
              sz = $
                t[ie].psf_sigma_pixels^2 + $
                (aratio*t[ie].psf_sigma_pixels)^2

              sz_in[ip, ie] = sz
              sz_meas[ip, ie] = t[ie].ixxmeas + t[ie].iyymeas

              fracuse[ip,ie] = float(t[ie].npsf_use)/t[ie].npsf
          ENDFOR 
                        
      ENDIF ELSE BEGIN 
          print,'FILE IS MISSING: '+file
      ENDELSE 


      irun = irun + 1
;      IF ( ( ( irun MOD 4) EQ 0 ) AND $
;           ( NOT keyword_set(dops) ) ) THEN key=prompt_kbrd()      
                                ;print,irun

  ENDFOR 

  ;; go back to defaults
  !p.multi=0

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Plot the noise in e1/e2 as a function of apix/seeing
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  psFile = inDir + front+'_noise_pixoverseeing.eps'
  IF keyword_set(dops) THEN BEGIN 
      begplot, name=psFile, /encap, /color
      myusersym,'fill_circle', symsize=0.7
  ENDIF

  setupplot

  colors=[!p.color, $
          !red, !darkGreen, !blue, $
          !magenta, !orange,$
          !grey50, !hotpink, !dodgerBlue, $
          !darkred, !seagreen, !purple, $
          !salmon, !turquoise, !sienna, $
          !ForestGreen, !FireBrick]

  ytitle = 'Noise e!D1!N'
  xtitle = 'pixel size (arcsec)'
  xrange = [0.24, 0.33]
  yrange = [1.e-4, 1.0]

  fracline = 1
  fraccolor = !p.color

  !p.multi=[0,0,2]

  se1 = sort(e1vals)
  se2 = sort(e2vals)
  e1str = strarr(nellip)
  e2str = strarr(nellip)

  xrange = [0.4, 1.6]/2.35
  xtitle = 'pixel size / seeing'
  ytitle = 'Noise e!D1!N'
  ;; e1 noise
  plot, apix/seeing, e1_sdev[*, se1[0]], $
    /ylog, $
    yrange=yrange, ystyle=1+2, $
    ytickf='loglabels', $
    xtitle=xtitle, ytitle=ytitle, $
    xrange=xrange, xstyle=1+2, $
    psym=8
  oplot, apix/seeing, e1_sdev[*, se1[0]]
  oplot, apix/seeing, fracuse[*, se1[0]], line=fracline, color=fraccolor

  e1str[0] = ntostr(e1vals[se1[0]],7,/round)
  iclr = 1
  FOR ie=1, nellip-1 DO BEGIN 
      oplot, apix/seeing, e1_sdev[*, se1[ie]], psym=8, $
        color=colors[iclr]
      oplot, apix/seeing, e1_sdev[*, se1[ie]], color=colors[iclr]

      e1str[ie] = ntostr(e1vals[se1[ie]],7,/round)
      iclr = iclr+1
  ENDFOR 

  legend,'Fraction used',$
    line=fracline,color=fraccolor,box=0,thick=!p.thick, $
    charsize=1

  legend, 'e!D1!N = '+e1str, colors=colors, line=replicate(0, nellip), $
          box=0, $
          charsize=1, thick=replicate(!p.thick, nellip),/right,/bottom

  ;; e2 noise
  ytitle = 'Noise e!D2!N'
  plot, apix/seeing, e2_sdev[*, se2[0]], $
    /ylog, $
    yrange=yrange, ystyle=1+2, $
    ytickf='loglabels', $
    xtitle=xtitle, ytitle=ytitle, $
    xrange=xrange, xstyle=1+2, $
    psym=8

  oplot, apix/seeing, e2_sdev[*, se2[0]]
  oplot, apix/seeing, fracuse[*, se2[0]], line=fracline, color=fraccolor

  e2str[0] = ntostr(e2vals[se2[0]],7,/round)
  iclr = 1
  FOR ie=1, nellip-1 DO BEGIN 
      oplot, apix/seeing, e2_sdev[*, se2[ie]], psym=8, $
        color=colors[iclr]
      oplot, apix/seeing, e2_sdev[*, se2[ie]], color=colors[iclr]

      e2str[ie] = ntostr(e2vals[se2[ie]],7,/round)
      iclr = iclr+1
  ENDFOR 

  legend, 'e!D2!N = '+e2str, colors=colors, line=replicate(0, nellip), $
          box=0, $
          charsize=1, thick=replicate(!p.thick, nellip)


  IF keyword_set(dops) THEN endplot, /trim_bbox

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Plot bias info
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF nellip GT 1 THEN BEGIN 

      psFile = inDir + front+'_bias_pixoverseeing.eps'
      IF keyword_set(dops) THEN BEGIN 
          begplot, name=psFile, /encap, /color
          myusersym, 'fill_circle',symsize=0.7
      ENDIF       
      setupplot


      e1color = !darkGreen
      e2color = !DodgerBlue
      
;      slope_yrange = [-0.01,0.01]
      slope_yrange = [-1,1]

      offset_yrange = [-1.e-4,1.e-4]
      
      slope_ytitle = 'b'
      offset_ytitle = 'a'

      key=prompt_kbrd()

      plot, $
        apix/seeing, e1_slopes, psym=8, $
        yrange=slope_yrange, ystyle=1+2, $
        xtitle=xtitle, ytitle=slope_ytitle

      oplot, apix/seeing, e1_slopes, psym=8, color=e1color
      oplot, apix/seeing, e1_slopes, color=e1color

      oplot, $
        apix/seeing, e2_slopes, psym=8, $
        color=e2color
      oplot, apix/seeing, e2_slopes, color=e2color
      oplot, apix/seeing, fracuse[*, se1[0]], line=fracline, color=fraccolor

      legend,'e!UMeas!N - e!UInput!N = a + b '+!csym.times+' e!UInput!N',$
        box=0,charsize=1
      legend,['e!D1!N', 'e!D2!N'], $
        line=[0,0],color=[e1color,e2color],box=0, $
        thick=[!p.thick,!p.thick],/right
      legend,'Fraction used',$
        line=fracline,color=fraccolor,box=0,thick=!p.thick, $
        charsize=1,/bottom

      plot, $
        apix/seeing, e1_offsets, psym=8, $
        yrange=offset_yrange, ystyle=1+2, $
        xtitle=xtitle, ytitle=offset_ytitle

      oplot, apix/seeing, e1_offsets, psym=8, color=e1color
      oplot, apix/seeing, e1_offsets, color=e1color

      oplot, $
        apix/seeing, e2_offsets, psym=8, $
        color=e2color
      oplot, apix/seeing, e2_offsets, color=e2color

      IF keyword_set(dops) THEN endplot, /trim_bbox

  ENDIF 

  !p.multi=0

  ;; Turns out, the size measurement is mostly independent of
  ;; the ellipticity. Just plot the first

  psFile = inDir + front+'_size_pixoverseeing.eps'
  IF keyword_set(dops) THEN BEGIN 
      begplot, name=psFile, /encap, /color
      myusersym, 'fill_circle',symsize=0.7
  ENDIF       
  setupplot

  key=prompt_kbrd()

  aplot, !gratio, $
    apix/seeing, sz_meas[*,0]/sz_in[*,0], $
    psym=8, $
    xtitle=xtitle, ytitle='size!UMeas!N / size!UInput!N'
  
  fraccolor=!blue
  fracline=0
  oplot, apix/seeing, fracuse, color=fraccolor, line=fracline
  legend,'Fraction used',$
    line=fracline,color=fraccolor,box=0,thick=!p.thick, $
    charsize=1,/bottom

  IF keyword_set(dops) THEN endplot, /trim_bbox  

  myusersym, 'fill_circle'
END 
