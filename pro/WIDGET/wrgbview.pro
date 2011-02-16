
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    RGBVIEW
;       
; PURPOSE:
;    Create and display a RGB image from red, green, and blue input images.
;
;    SDSS should use red=i, grn=r, blue=g although this will have false
;    color.  The result is similar to Steve Kent's images.
;    If the device is postscript, then a color postscript is made.
;    NOTE: The color map is NOT inverted so there is often large amounts
;          of black space.
;
;    Works best on displays with 16 million+ colors.
;    At the least the program requires 256 colors to display properly.  
;    To insure the you get all the colors, you can request a private color map.  
;    This is done when the first window is opened.  You must use this command:
;         IDL> window, colors=256
;
;    It won't work except on the first window so you might consider putting
;    it in your .idl.startup file.  Note that using a private color map
;    may cause "flashing" when you point in the window.
;
;    Also note that the number of colors on your display does not affect
;    the optional jpeg/gif files that RGBVIEW outputs because it uses full
;    color range when producing them.
;
;    WARNING: this program can use lots of memory. It makes 
;    2 copies of each image for speed: one that converted to float and
;    sky subtracted and then yet another that is scaled. If your images
;    are already sky-subtracted and float then you should modify this program
;    so it doesn't make 1 of the extra copies.
;
; CALLING SEQUENCE:
;     rgbview4, rimage, gimage, bimage, sdss=sdss, 
;             gamma=gamma, saturation=saturation, 
;             contrast=contrast, low_cut=low_cut, 
;             nodisplay=nodisplay, 
;             jpegfile=jpegfile, giffile=giffile, 
;             rsig=rsig, gsig=gsig, bsig=bsig, 
;             rsky=rsky, gsky=gsky, bsky=bsky, 
;             rmap=rmap, gmap=gmap, bmap=bmap, 
;             color_im=color_im, ar=ar, ag=ag, ab=ab, imtot=imtot, 
;             xrange=xrange, yrange=yrange, noprompt=noprompt, $
;             title=title, xtitle=xtitle, ytitle=ytitle, subtitle=subtitle, 
;             noframe=noframe, nolabels=nolabels, 
;             _extra=extra_key
;
; INPUTS: 
;    red, grn, blue: The red, green and blue images.  Images must be same size.
;
; OPTIONAL INPUTS:
;    gamma: gamma value.  typical values are 1.5-2.5, default is 2.5.  
;           rescales images to (image)^(1./gamma) to bring out low surface
;           brightness features.
;    saturation: color saturation. Bring out the red and blue by scaling 
;                images: red = red * (red/green)^(saturation)
;                        blue = blue * (blue/green)^(saturation)
;           default value is zero.  Be careful with this, it only works 
;           well for very high S/N images, and tends to make things
;           look a little unrealistic (see Frei and Gunn..)
;    contrast:  The number of sigma above the mean to use in images.
;            The default is 15 but good results depend on the image.
;            For example, a really bright galaxy will require contrast
;            to be set _many_ sigma above mean.
;            You will be prompted for a change in contrast unless /noprompt
;            is set (good for batch processing)
;    low_cut: Lowest value in image to use in # of sigma above mean.  Images
;            are sky subtracted, so this number should be greater than zero.
;            Default is 0.7sigma
;    jpegfile: If a string is sent in this keyword, rgbview will write a 
;          jpeg file containing the image with this filename.  If you have
;          8-bit display, then this is better than
;          the display because it is compressed 24-bit. This file is
;             created directly from the images, so no axes or labels will 
;             appear.  To save axes, read from the screen: 
;             write_gif, file, tvrd(), rmap, gmap, bmap
;    giffile: same as above except writes gif (8-bit) image file.  This is
;             created directly from the images, not read from the screen, so 
;             no axes or labels will appear.  
;    rsig, gsig, bsig, : input value for sky noise in r image. Saves time by 
;          avoiding sigma-clipping the images.  If not input, these values
;          can be returned through these keywords.
;    rsky, gsky, bsky: input sky value of r,g,b images.  Avoides running 
;          sigma-clipping to find sky.  If not input, these values
;          can be returned through these keywords.
;    xrange, yrange, noframe, nolabels: see tvim2
;    title,xtitle,ytitle,subtitle: Plot labels.
;    _extra=extra_key:  Other plotting keywords.
;
; KEYWORD PARAMETERS:
;    sdss: if set, rescales images by energy based on sdss filters.  Doesn't
;       always result in better _looking_ images.
;    nodisplay: if set, no display is made.  You might use this if you are
;          just outputting the jpeg files, maybe in a batch job if used in
;          conjunction with /noprompt
;    noprompt: if set then don't ask for a change of contrast, etc
;       
; OPTIONAL OUTPUTS: 
;    rmap,gmap,bmap: color map vectors.  These are the vectors used to display 
;           this image.  If using an 8-bit display, they can be sent to 
;
;                IDL> WRITE_GIF, filename, TVRD(), rmap, gmap, bmap
;           
;           If on 8-bit display, you might need to go back to BW linear display.
;           The color map can be reset to BW linear with loadct,0
; 
;    color_im: a byte 2-d image containing the image used for display with the 
;        color maps above.
;    ar, ag, ab: red, green, blue byte scaled images used to make full 24-bit
;         image (as used to output jpeg images)
;    imtot: image containing ar, ag, ab in the form bytarr(3, n, m).  Can be 
;         input directly to write_jpeg to produce 24-bit image.
;  
; CALLED ROUTINES: (lower case are built in IDL procedures)
;    DCENTER
;    SIGMA_CLIP
;    color_quan
;    bytscl
;    tv
;    tvlct
;    (device)
;    (write_gif)
;    (write_jpeg)
;
; PROCEDURE: 
;    
;
; REVISION HISTORY:
;    Author: Erin Scott Sheldon  UofMich  11/28/99
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO setup_plot, nx, ny, aspect, black, xsize, ysize, px, py, xrng, yrng

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set up plot
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  plot, [0,1],[0,1],/nodata,xstyle=4,ystyle=4, background=black
  px=!x.window*!d.x_vsize
  py=!y.window*!d.y_vsize
  xsize=px[1]-px[0]
  ysize=py[1]-py[0]
  
  IF xsize GT ysize*aspect THEN xsize=ysize*aspect ELSE ysize=xsize/aspect 
  px[1]=px[0]+xsize
  py[1]=py[0]+ysize

  nxm=nx-1
  nym=ny-1

  IF n_elements(xrange) EQ 0 THEN BEGIN
      xrng=[ -0.5, nxm+0.5]
  ENDIF ELSE BEGIN
      xrng=[xrange(0), xrange(n_elements(xrange)-1)]
  ENDELSE 

  IF n_elements(yrange) EQ 0 THEN BEGIN
      yrng = [-0.5,nym+0.5]
  ENDIF ELSE BEGIN
      yrng = [yrange(0), yrange(n_elements(yrange)-1)]
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; center up the display
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  dcenter, xsize, ysize, px, py, /silent

  ;; PS device?
  IF !d.name EQ 'PS' THEN BEGIN 
      ;; this is white on black, make sure we use small
      ;; thickness
      !p.thick=1
      !x.thick=1
      !y.thick=1
      !p.charsize=1
      !p.charthick=1
  ENDIF 


END 

PRO wrgbview, rr, gg, bb, sdss=sdss, gunn=gunn, $
              gamma=gamma, saturation=saturation, $
              contrast=contrast, low_cut=low_cut, $
              nodisplay=nodisplay, $
              rsig=rsig, gsig=gsig, bsig=bsig, $
              rsky=rsky, gsky=gsky, bsky=bsky, $
              rmap=rmap, gmap=gmap, bmap=bmap, $
              color_im=color_im, ar=ar, ag=ag, ab=ab, imtot=imtot, $
              xrange=xrange, yrange=yrange, noprompt=noprompt, $
              title=title, xtitle=xtitle, ytitle=ytitle, subtitle=subtitle, $
              noframe=noframe, nolabels=nolabels, $
              jpegfile=jpegfile, giffile=giffile, pngfile=pngfile, $
              tvread=tvread, $
              _extra=extra_key


  IF n_params() EQ 0 THEN BEGIN
      print,'-Syntax: rgbview, rimage, gimage, bimage, '
      print,' sdss=sdss, gunn=gunn, '
      print,' gamma=gamma, saturation=saturation, '
      print,' contrast=contrast, low_cut=low_cut, '
      print,' /nodisplay, '
      print,' jpegfile=jpegfile, giffile=giffile, '
      print,' /tvread, '
      print,' rsig=rsig, gsig=gsig, bsig=bsig, '
      print,' rsky=rsky, gsky=gsky, bsky=bsky, '
      print,' /noprompt, '
      print,' /noframe, /nolabels, '
      print,' rmap=rmap, gmap=gmap, bmap=bmap, '
      print,' color_im=color_im, ar=ar, ag=ag, ab=ab, imtot=imtot, '
      print,' _extra=extra_key'
      print
      print,' -For sloan, use red=i, grn=r, blue=g'
      print,' -contrast=nsigma above mean'
      return
  ENDIF 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set up parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; titles and stuff for axes.
  IF n_elements(title) EQ 0 THEN title = ''
  IF n_elements(xtitle) EQ 0 THEN xtitle=''
  IF n_elements(ytitle) EQ 0 THEN ytitle=''
  IF n_elements(subtitle) EQ 0 THEN subtitle=''

  ;; see how many colors we got
  max_color=!d.n_colors-1
  IF max_color LT 255 THEN BEGIN
      print
      print,'Only got ',ntostr(max_color+1),' colors for display'
      print
  ENDIF 

  IF max_color GT 10000000 THEN true_color=1 ELSE true_color=0
  IF ( (n_elements(jpegfile) NE 0) AND keyword_set(tvread) AND $
       (NOT true_color) ) THEN BEGIN 
      message,'You cannot read a jpeg from the display in 8-bit color mode',/inf
      message,'Try using a gif and /tvread'
  ENDIF 
  
  IF ( ( (n_elements(giffile) NE 0) OR (n_elements(jpegfile) NE 0) ) AND $
       keyword_set(tvread) AND keyword_set(nodisplay) ) THEN BEGIN 
      message,'Cannot read from display if /nodisplay is set'
  ENDIF 

  ;; gamma level
  IF n_elements(gamma) EQ 0 THEN gam = 2.2 ELSE gam = gamma > 1.

  ;; saturation level
  IF n_elements(saturation) EQ 0 THEN sat=0. ELSE sat=saturation > 0.

  ;; lower contrast level
  IF n_elements(low_cut) EQ 0 THEN lcut = .7 ELSE lcut = low_cut > .1 

  ;; upper contrast level
  IF n_elements(contrast) EQ 0 THEN cont = 15. ELSE BEGIN 
      cont = contrast > lcut
  ENDELSE 


  cont = float(cont)
  lcut = float(lcut)

  print,'Gamma = ',ntostr(gam,5)
  print,'Saturation = ',ntostr(sat,5)
  print,'Low level = ',ntostr(lcut,5),' sigma'
  print,'Beginning contrast = ',ntostr(cont,5),' sigma'

  ;; prompting
  IF NOT keyword_set(noprompt) THEN noprompt = 0

  ;; Check current device
  IF (!d.flags AND 1) EQ 0 THEN doX = 1 ELSE doX = 0

  ;; Check size of arrays
  szr = size(rr)
  szg = size(gg)
  szb = size(bb)

  IF (szr[4] NE szg[4]) OR (szr[4] NE szb[4]) THEN BEGIN
      print,'Arrays must be of same size'
      return
  ENDIF 

  nx = szr[1]
  ny = szr[2]
  aspect = float(nx)/ny

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Determine Sky and Skysig
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF (n_elements(rsig) EQ 0) OR (n_elements(rsky) EQ 0) $
  THEN dosig=1 ELSE dosig=0

  time=systime(1)
  IF dosig THEN BEGIN 
      print,'Sigma Clipping'
      niter = 2
      nsig = 3.5
      sigma_clip, rr, sky_r, sigmar, niter=niter, nsig=nsig, /silent
      sigma_clip, gg, sky_g, sigmag, niter=niter, nsig=nsig, /silent
      sigma_clip, bb, sky_b, sigmab, niter=niter, nsig=nsig, /silent
  ENDIF

  IF n_elements(rsig) EQ 0 THEN BEGIN
      rsig = sigmar
      gsig = sigmag
      bsig = sigmab
  ENDIF 
  IF n_elements(rsky) EQ 0 THEN BEGIN
      rsky = sky_r
      gsky = sky_g
      bsky = sky_b
  ENDIF 

  print
  print,'rsky ',rsky,' sigma(rimage) ',rsig
  print,'gsky ',gsky,' sigma(gimage) ',gsig
  print,'bsky ',bsky,' sigma(bimage) ',bsig

  ggain = 1.
  rgain = 1.
  igain = 1.
  gfreq = 1.
  rfreq = 1.
  ifreq = 1.

  bgreat = 0.
  rless = 2.
  IF keyword_set(sdss) THEN BEGIN 
      ;; scale to sloan colors
      print
      print,'Scaling to energy using SDSS bands'
      gfreq = 1./4700.
      rfreq = 1./6250.
      ifreq = 1./7700.
  ENDIF ELSE IF keyword_set(gunn) THEN BEGIN 
      ;; Gunn stuff
      print
      print,'Using Gunn bands'

      bgreat = 1.
      gfreq = 1./4930.
      rfreq = 1./6550.
      ifreq = 1./8200.
  ENDIF 

  ;; sky subtract
  red = float(rr) - rsky
  grn = float(gg) - gsky
  blue = float(bb) - bsky

  ;; scale by gain.  Also, scale to relative energy
  IF (ggain NE 1.) OR (ifreq NE 1.) THEN BEGIN 
      red  = red*ifreq/rfreq/igain
      IF rgain NE 1. THEN grn  = grn/rgain
      blue = blue*gfreq/rfreq/ggain
  ENDIF 

  ;; Convert so color of sun is white (or grey)

  IF keyword_set(sdss) OR keyword_set(gunn) THEN BEGIN 
      fac1=10.^(-0.09/2.5)
      fac2=10.^(-0.02/2.5)
      grn = grn*fac1
      red = red*fac2
  ENDIF 
  maxr = max(red)
  maxg = max(grn)
  maxb = max(blue)
  minr = min(red)
  ming = min(grn)
  minb = min(blue)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; byte scale each of the images
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  continue = 1

  WHILE continue DO BEGIN 

      power = 1./gam

      _b = blue
      _g = grn
      _r = red

      ;; highest value
      highr = cont*rsig < maxr/2.
      highg = cont*gsig < maxg/2.
      highb = cont*bsig < maxb/2.

      ;; lowest value
      lowr = lcut*rsig
      lowg = lcut*gsig
      lowb = lcut*bsig

      high = max([highr, highg, highb]) 
      low  = max([lowr, lowg, lowb])

      ;; divide so high stuff is relatively lowered by gamma correction
      ;; and low stuff is raised.

      mid = (high + low)/2.

      high = high/mid
      low = low/mid
      _b = _b/mid
      _g = _g/mid
      _r = _r/mid

      IF (gam GT 1.) OR (sat GT 0.) THEN BEGIN 
          wr = where(_r GT low AND _r LT high, nwr)
          wg = where(_g GT low AND _g LT high, nwg)
          wb = where(_b GT low AND _b LT high, nwb)
          IF nwr NE 0 THEN BEGIN
              IF sat NE 0. THEN BEGIN
                  _r[wr]=_r[wr]*(_r[wr]/(_g[wr]>low<high) < rless > 0.)^sat
              ENDIF 
              IF gam GT 1. THEN _r[wr] = _r[wr]^power
          ENDIF 
          
          IF nwg NE 0 THEN BEGIN 
              IF gam GT 1. THEN _g[wg] = _g[wg]^power
          ENDIF 
          
          IF nwb NE 0 THEN BEGIN 
              IF sat NE 0. THEN BEGIN
                  _b[wb] = _b[wb]*(_b[wb]/(_g[wb]>low<high) > bgreat )^sat
              ENDIF 
              IF gam GT 1. THEN _b[wb] = _b[wb]^power
          ENDIF 
      ENDIF 

      low  = low^power
      high = high^power

      ar = bytscl(_r, min=low, max=high, top=255)
      ag = bytscl(_g, min=low, max=high, top=255)
      ab = bytscl(_b, min=low, max=high, top=255)

      imtot = bytarr(3, nx, ny )
      
      imtot[0,*,*] = ar
      imtot[1,*,*] = ag
      imtot[2,*,*] = ab

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Create pseudo color image for display
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      IF (max_color LE 255) OR (n_elements(giffile) NE 0) THEN BEGIN 
          color_im = color_quan(ar, ag, ab, rmap, gmap, bmap)
      ENDIF 
      IF ( (n_elements(giffile) NE 0) AND (max_color LT 255) AND $
           (NOT keyword_set(tvread) ) )  THEN BEGIN 
          ;; gif not limited by available colors
          altcolor_im = color_quan(ar, ag, ab, $
                                   altrmap, altgmap, altbmap, $
                                   colors=256)
      ENDIF 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output to proper device (unless nodisplay)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      IF NOT keyword_set(nodisplay) THEN BEGIN 

          ;; special situation if 8-bit X display
          IF (max_color LE 255) AND doX THEN BEGIN 

              ww = where(ar EQ 0, nww)
              IF nww EQ 0 THEN BEGIN 
                  ww = where(ag eq 0, nww)
                  IF nww EQ 0 THEN BEGIN
                      ww = where(ab EQ 0, nww)
                      IF nww EQ 0 THEN BEGIN
                          ww = where(ar EQ min(ar), nww)
                      ENDIF 
                  ENDIF 
              ENDIF 
              black = color_im[ww[0]]
          
              ww = where(ar EQ 255 AND ag EQ 255 AND ab EQ 255, nww)
              IF nww EQ 0 THEN BEGIN 
                  print,'redo'
                  ww = where(rr EQ max(rr))
              ENDIF 
              white = color_im[ww[0]] 
          ENDIF ELSE BEGIN 
              defsysv, '!white', exist=wexist
              IF wexist THEN BEGIN 
                  white=!white
                  black=!black
              ENDIF ELSE BEGIN
                  ;; try something else
                  white=!d.n_colors-1
                  black=0L
              ENDELSE 
          ENDELSE 

          setup_plot, nx, ny, aspect, black, xsize, ysize, px, py, xrng, yrng
          IF max_color GT 255 THEN BEGIN 
              delvarx,black,white
          ENDIF 

          IF doX THEN BEGIN     ;X window
              IF max_color LE 255 THEN BEGIN ;; 8-bit display
                  tv, congrid(color_im, xsize, ysize), px[0],py[0]
                  tvlct, rmap, gmap, bmap
                  pos = [px[0], py[0], px[1], py[1]]
              ENDIF ELSE BEGIN 
                  tv, congrid(imtot, 3, xsize, ysize),true=1,px[0],py[0]
                  pos = [px[0], py[0], px[1], py[1]]
              ENDELSE 
          ENDIF ELSE BEGIN      ;Postscript
              device,/color
              tvlct,indgen(256),indgen(256),indgen(256)

              ;; add black background
              POLYFILL, [1,1,0,0,1], [1,0,0,1,1], /NORMAL, COLOR=0
              pos = [px[0], py[0], px[1], py[1]]
              tv, imtot,true=1,px[0],py[0], xsize=xsize, ysize=ysize, /device
          ENDELSE 

          IF keyword_set(noframe) OR keyword_set(nolabels) THEN BEGIN 
              plot, [0,0], [0,0], xstyle=5, ystyle=5, $
                title=title,xtitle=xtitle,ytitle=ytitle, subtitle=subtitle, $
                xrange=xrng, yrange=yrng, position=pos, color=white, $
                /noerase, /device, /nodata, background=black
          ENDIF ELSE BEGIN 
              plot, [0,0], [0,0], xstyle=1, ystyle=1, $
                    title=title, xtitle=xtitle, ytitle=ytitle, subtitle=subtitle, $
                    xrange=xrng, yrange=yrng, position=pos, color=white,$
                    /noerase, /device, /nodata, background=black
          ENDELSE 
          
          IF (NOT keyword_set(noframe)) AND keyword_set(nolabels) THEN BEGIN 
              axis,xaxis=1,xtickname=strarr(10)+" ",color=white
              axis,xaxis=0,xtickname=strarr(10)+" ",color=white
              axis,yaxis=1,ytickname=strarr(10)+" ",color=white
              axis,yaxis=0,ytickname=strarr(10)+" ",color=white
          ENDIF  

          IF (NOT noprompt) AND doX THEN BEGIN 
              print
              ans = ' '
              print,format='($, "Command: (g gamma) (s saturation) (c contrast) (l lower cut)")'
              read,ans
              ans = STRLOWCASE(ntostr(ans,1,0))
              CASE ans OF 
                  'g': BEGIN
                          print,format='($, "New gamma")'
                          read, gam
                          gam = gam > 1.
                       END 
                  's': BEGIN
                          print,format='($, "New saturation")'
                          read, sat
                          sat = sat > 0.
                       END
                  'c': BEGIN
                          print,format='($, "New contrast")'
                          read, cont
                          cont = cont > lcut
                       END
                  'l': BEGIN
                          print,format='($, "New low cut")'
                          read, lcut
                          lcut = lcut > .1
                       END 
                  ELSE: continue = 0
              ENDCASE 
          ENDIF ELSE continue = 0 
      ENDIF ELSE continue = 0; nodisplay

  ENDWHILE 
  print

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; output the image if requested
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(jpegfile) NE 0 THEN BEGIN 

      IF n_elements(quality) EQ 0 THEN qual=75 ELSE BEGIN
          qual = quality < 100 > 0
      ENDELSE 

      print
      print,'Writing jpeg file: ',jpegfile
      print,'Using quality: ',ntostr(qual)

      IF keyword_set(tvread) THEN BEGIN
          write_jpeg,jpegfile,tvrd(true=1),true=1,quality=qual 
      ENDIF ELSE BEGIN
          write_jpeg,jpegfile,imtot,true=1,quality=qual
      ENDELSE 
  ENDIF 

  IF n_elements(pngfile) NE 0 THEN BEGIN 

      print,'Writing png file: ',pngfile

      ;; get image from display?
      IF keyword_set(tvread) THEN BEGIN 
          IF true_color THEN BEGIN 
              tmp = tvrd(true=1)
              
              IF float(!version.release) LT 5.4 THEN BEGIN 
                  ;; In old version png is written out flipped
                  tmp[0,*,*] = rotate(rotate(reform(tmp[0,*,*]),1),4)
                  tmp[1,*,*] = rotate(rotate(reform(tmp[1,*,*]),1),4)
                  tmp[2,*,*] = rotate(rotate(reform(tmp[2,*,*]),1),4)
              ENDIF 

              write_png, pngfile, tmp
          ENDIF ELSE BEGIN 
              write_png, pngfile, tvrd(), rmap, gmap, bmap
          ENDELSE 
      ENDIF ELSE BEGIN 

          IF float(!version.release) LT 5.4 THEN BEGIN 
              ;; In old version png is written out flipped
              imtot[0,*,*] = rotate(rotate(reform(imtot[0,*,*]),1),4)
              imtot[1,*,*] = rotate(rotate(reform(imtot[1,*,*]),1),4)
              imtot[2,*,*] = rotate(rotate(reform(imtot[2,*,*]),1),4)
          ENDIF 
          
          write_png, pngfile, imtot

          IF float(!version.release) LT 5.4 THEN BEGIN 
              ;; In old version png is written out flipped
              imtot[0,*,*] = rotate(rotate(reform(imtot[0,*,*]),1),4)
              imtot[1,*,*] = rotate(rotate(reform(imtot[1,*,*]),1),4)
              imtot[2,*,*] = rotate(rotate(reform(imtot[2,*,*]),1),4)
          ENDIF 
          
      ENDELSE 

  ENDIF 

  IF n_elements(giffile) NE 0 THEN BEGIN 
      print,'Writing gif file: ',giffile

      ;; get image from display
      IF keyword_set(tvread) THEN BEGIN 
          IF true_color THEN BEGIN 
              tmp = tvrd(true=1)
              colorim = reform( color_quan(tmp[0,*,*], tmp[1,*,*], tmp[2,*,*], $
                                           ttrmap, ttgmap, ttbmap, $
                                           colors=256) )
              write_gif, giffile, colorim, ttrmap, ttgmap, ttbmap
          ENDIF ELSE BEGIN 
              write_gif, giffile, tvrd(), rmap, gmap, bmap
          ENDELSE 
      ENDIF ELSE BEGIN 
          IF n_elements(altcolor_im) EQ 0 THEN BEGIN 
              write_gif, giffile, color_im, rmap, gmap, bmap
          ENDIF ELSE BEGIN
              write_gif, giffile, altcolor_im, altrmap, altgmap, altbmap
          ENDELSE 
      ENDELSE 
  ENDIF 

return
END 
