
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    OBJECT_VIEW
;       
; PURPOSE:
;    Display the atlas images and a color composite for an object or list
;    of objects in sequence.
;
; CALLING SEQUENCE:
;    object_view, struct, contrast=contrast, saturation=saturation, gamma=gamma,$
;                addu=addu, imtot=imtot
;
; INPUTS: 
;    struct: PHOTO structure containing the tags required by GET_ATLAS.
;
; OPTIONAL INPUTS:
;    contrast: The number of sigma above the mean to use in images.
;            Sent to RGBVIEW.
;            The default is 15 but good results depend on the image.
;            For example, a really bright galaxy will require contrast
;            to be set _many_ sigma above mean.
;            You will be prompted for a change in contrast unless /noprompt
;            is set (good for batch processing)
;    saturation: color saturation. Sent to RGBVIEW.
;                Bring out the red and blue by scaling 
;                images: red = red * (red/green)^(saturation)
;                        blue = blue * (blue/green)^(saturation)
;           default value is zero.  Be careful with this, it only works 
;           well for very high S/N images, and tends to make things
;           look a little unrealistic (see Frei and Gunn..)
;    gamma: gamma value sent to RGBVIEW.  
;           typical values are 1.5-2.5, default is 2.5.  
;           rescales images to (image)^(1./gamma) to bring out low surface
;           brightness features.
;
; KEYWORD PARAMETERS:
;    /addu: add the u-band to the color image.  
;       
; OUTPUTS: 
;    NONE
;
; OPTIONAL OUTPUTS:
;    imtot: The true color image of the last object displayed.
;
; CALLED ROUTINES:
;    FETCH_DIR
;    GET_ATLAS
;    DISPLAY_ATLAS
;    RGBVIEW
; 
; PROCEDURE: 
;    
;	
;
; REVISION HISTORY:
;    Created 3-Jul-2002, Erin S. Sheldon  UofMich
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO object_view_fixim, imu, img, imr, imi          

  ;; these are ints.  Max is 32767
  maxval = 32767L
  
  imu = long(imu)
  img = long(img)
  imr = long(imr)
  imi = long(imi)
  
  wu = where(imu LT 0, nwu)
  IF nwu NE 0 THEN BEGIN 
      imu[wu] = maxval + (imu[wu] - (-32767) )
  ENDIF 
  wg = where(img LT 0, nwg)
  IF nwg NE 0 THEN BEGIN 
      img[wg] = maxval + (img[wu] - (-32767) )
  ENDIF 
  wr = where(imr LT 0, nwr)
  IF nwr NE 0 THEN BEGIN 
      imr[wr] = maxval + (imr[wu] - (-32767) )
  ENDIF 
  wi = where(imi LT 0, nwi)
  IF nwi NE 0 THEN BEGIN 
      imi[wi] = maxval + (imi[wu] - (-32767) )
  ENDIF 
END 

; This copied from get_atlas_display
pro object_view_display, objstruct, imtot, i, ret, zoom, run, camcol, $
                   field, obj_id, cnum, rerun=rerun, sep=sep, pa=pa, $
                   noprompt=noprompt,silent=silent,$
                   hideradec=hideradec,_extra=extra


   IF NOT keyword_set(silent) THEN silent=0

   ;; find a magnitude array

   domag = 1b
   tags = tag_names(objstruct)
   wmag = where(tags EQ 'COUNTS_MODEL',nmod)
   IF nmod EQ 0 THEN BEGIN 
       wmag = where(tags EQ 'PETROCOUNTS', npet)
       IF npet EQ 0 THEN BEGIN 
           wmag = where(tags EQ 'FIBERCOUNTS',nfib)
           IF nfib EQ 0 THEN BEGIN 
               wmag = where(tags EQ 'PSFCOUNTS',npsf)
               IF npsf EQ 0 THEN BEGIN 
                   ;;print
                   ;;print,'No magnitude arrays found. Not displaying colors.'
                   domag = 0b
               ENDIF 
           ENDIF 
       ENDIF 
   ENDIF 

   ;; check for ra/dec
   wra = where(tags EQ 'RA')
   wdec = where(tags EQ 'DEC')

   IF wra[0] EQ -1 OR wdec[0] EQ -1 THEN BEGIN
       ;;print
       ;;print,'No RA/DEC found. Not displaying position'
       dora=0b 
   ENDIF ELSE dora=1b

;   IF (NOT silent) AND (wmag[0] NE -1) THEN BEGIN
;       print,'GET_ATLAS_DISPLAY: Using '+tags[wmag]
;       print,'------------------------------------------------'
;   ENDIF 

   IF n_elements(rerun) NE 0 THEN rerstr = '-'+ntostr(rerun) ELSE rerstr=''
   IF n_elements(sep) EQ 0 THEN sep = -1
   IF n_elements(pa) EQ 0 THEN pa = -1
   
   ;; We have had problems fitting all the info on the screen! 
   ;; Need an adaptive character size.
   char_old = !p.charsize

   px = !p.multi[1]
   py = !p.multi[2]

   CASE px OF 
       0.: !p.charsize= 1.
       1.: !p.charsize= 1.
       2.: !p.charsize= 1.
       ELSE: !p.charsize = .7
   ENDCASE 

   CASE py OF
       3: !p.charsize=1.75*!p.charsize
       4: !p.charsize=1.4*!p.charsize
       5: !p.charsize=1.4*!p.charsize
       ELSE: 
   ENDCASE

   xtitle = ''
   subtitle = ''
   IF domag THEN BEGIN 
       c = objstruct.(wmag[0])
       mag=['','','','','']
       color = ['u','g','r','i','z']

       FOR jj=0,n_elements(cnum)-1 DO BEGIN 
           mag[cnum[jj]]=strmid( strtrim(string(c[cnum[jj]]),2), 0, 5)
       ENDFOR 

       FOR jj=0,4 DO BEGIN 
           IF (mag[jj] NE '') THEN BEGIN 
               IF (xtitle NE '') THEN xtitle=xtitle+'  '
               xtitle=xtitle + color[jj]+'='+mag[jj]
           ENDIF 
       ENDFOR 
   
       FOR jj=0,3 DO BEGIN 
           IF (mag[jj] NE '' AND mag[jj+1] NE '') THEN BEGIN 
               diff=strmid(strtrim(string(c[jj] - c[jj+1]),2), 0, 5)
               subtitle=subtitle+'  '+color[jj]+'-'+color[jj+1]+'='+diff
           ENDIF 
       ENDFOR 
       xtitle=xtitle + '  ' + subtitle

   ENDIF 

   ;; changed order to run-rerun-camcol-field-id
   ;;title = run2string(run)+'-'+ntostr(camcol)+rerstr
   title = run2string(run)+rerstr+'-'+ntostr(camcol)
   title = title+ '-'+field2string(field)+'-'+ntostr(obj_id)
   
   if (not keyword_set(hideradec)) AND dora then begin
       radecstr, objstruct.ra, objstruct.dec, ra, dec
       title=title+ '  '+ra+'  '+dec
   endif
       
   IF sep NE -1 THEN BEGIN
       sepst=strtrim(sep,2)     ;take only one digit after decimal place
       title=title+' sep: '+strmid(sepst,0,strpos(sepst,'.')+2)
       title=title+' pa: '+strtrim(string(fix(pa)),2)
   ENDIF

   zoom_old = zoom
   xtitle_old=xtitle
   sub_old = subtitle

   sigma_clip, imtot, mean, sigma, nsig=3.5, niter=2,/silent
   rdis_setup, imtot, pls
   IF sigma LT 5. THEN BEGIN
       low = mean
       high = mean + 25.
   ENDIF ELSE BEGIN
       low = mean
       high = mean + 10.*sigma
   ENDELSE 
   pls.low = low
   pls.high = high
   ;; Display and prompt user.  Zoom if requested.
   REPEAT BEGIN

       rdis, imtot, pls, zoom=zoom, silent=1, $
         xtitle=xtitle,title=title,_extra=extra

       zoom = zoom_old

       IF (noprompt EQ 0) THEN BEGIN 
           print,''
           print,'(n for next) (p for previous) (r to redisplay) (q to exit)'
           redo = strlowcase( get_kbrd(1) )
       ENDIF ELSE redo='n'    

       IF (redo EQ 'r') THEN BEGIN 
           print,'------------------------------------------------'
           print,'Want to zoom(y/n)?  n to display original.'
           zm =  strlowcase( get_kbrd(1) )
           print,''
           IF (zm EQ 'y') THEN BEGIN            ;Zooming
               zoom = 1
               xtitle=''
               subtitle=''
           ENDIF ELSE BEGIN                     ;Redisplay original
               rdis_setup, imtot, pls  &  pls.low = low  &  pls.high = high
               xtitle=xtitle_old
               subtitle=sub_old
           ENDELSE 
       ENDIF ELSE IF (redo EQ 'p') THEN begin
           IF (i EQ 0) THEN BEGIN               ;At the first object.
               print,'------------------------------------------------'
               print,'No such object. Starting over'
               print,'------------------------------------------------'
               tmp = get_kbrd(1)
               i = i-1
           ENDIF ELSE BEGIN 
               i = i-2                          ;Going to previous object.
           ENDELSE 
           redo ='n'
       ENDIF ELSE IF (redo EQ 'q') THEN BEGIN 
           print,'Quitting get_atlas'
           ret = 'q'
           !p.charsize=char_old
           return
       ENDIF ELSE BEGIN 
           redo = 'n'
           xtitle=xtitle_old
           subtitle=sub_old
       ENDELSE 
    
   ENDREP UNTIL (redo EQ 'n')
  
   ret='next'
   zoom = zoom_old

   !p.charsize=char_old
   return 
END 

PRO object_view, struct, $
                 contrast=contrast, saturation=saturation, gamma=gamma,$
                 noprompt=noprompt, $
                 addu=addu, imtot=imtot, dojpeg=dojpeg, outdir=outdir

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: object_view, struct, contrast=contrast, saturation=saturation, gamma=gamma, addu=addu'
      print,' struct should be a PHOTO object struct or array of structures.'
      return
  ENDIF 

  IF n_elements(outdir) EQ 0 THEN outdir=''

  ;; just choose some feducial numbers for
  ;; skysig since much of atlas images are exactly 1000

  isig   =    8.07816
  rsig   =    6.71133
  gsig   =    6.03404
  isky = 1000.
  rsky = 1000.
  gsky = 1000.


  nstruct=n_elements(struct)

  !p.multi = [0,1,2]
  
  cont=1
  i=0L
  WHILE cont DO BEGIN 
;  FOR i=0L, nstruct-1 DO BEGIN 

      fetch_dir, struct[i].run, struct[i].camcol, $
        struct[i].rerun, dir, atldir,$
        /check
      
      IF dir[0] NE '' THEN BEGIN 
                                
          get_atlas, struct, i, dir=atldir, $
            imu=imu, img=img, imr=imr, imi=imi,$
            imtot=aimtot, $
            /noprompt,/nodisplay, row0=row0, col0=col0

          print,row0
          print,col0
          

          IF n_elements(imi) NE 0 THEN BEGIN 

              object_view_fixim, imu, img, imr, imi

;              help,imu,img,imr,imi
              
;              sz = size(imu, /dim)
;              nx = sz[0]
;              ny = sz[1]
;              imu = rebin(imu, nx/2, ny/2) 
;              img = rebin(img, nx/2, ny/2) 
;              imr = rebin(imr, nx/2, ny/2) 
;              imi = rebin(imi, nx/2, ny/2) 
;              help,imu,img,imr,imi

              IF keyword_set(addu) THEN BEGIN 
                  img = (imu + img) - 1000
              ENDIF 

              rgbview_lup, imi, imr, img, /noprompt, $
                contrast=contrast, saturation=saturation, gamma=gamma,$
                imtot=imtot, $
                rsky=isky, gsky=rsky, bsky=gsky, jpegfile='~/tmp/test.jpg',$
                rsig=isig, gsig=rsig, bsig=gsig, /sdss;,/silent
                            
              object_view_display, struct[i], aimtot, i,ret,0,$
                struct[i].run,struct[i].camcol,$
                struct[i].field,struct[i].id,[0,1,2,3,4],$
                rerun=struct[i].rerun,/noprompt,/silent;,$
                ;color=!white,invbw=0

              IF keyword_set(dojpeg) AND !d.name EQ 'X' THEN BEGIN 
                  jpegfile = outdir + $
                    'obj-'+run2string(struct[i].run)+'-'+$
                    ntostr(struct[i].rerun)+'-'+$
                    ntostr(struct[i].camcol)+'-'+$
                    field2string(struct[i].field)+'-'+$
                    ntostr(struct[i].id)+'.jpg'

                  write_jpeg, jpegfile, tvrd(true=1), true=1
                  
              ENDIF 

          ENDIF 
      ENDIF ELSE BEGIN 
          print,'Cannot view object '+ntostr(i)
      ENDELSE 
      IF !d.name EQ  'X' THEN BEGIN

          IF keyword_set(noprompt) THEN BEGIN 
              i = i+1
          ENDIF ELSE BEGIN 
              IF i EQ (nstruct-1) THEN print,'This is the last object'
              print,'(q to quit) (p for previous) (any key for next)'
              key=get_kbrd(1)
              CASE strlowcase(key) OF
                  'q': return
                  'p': i = i-1
                  ELSE: i = i+1
              ENDCASE 
          ENDELSE 

      ENDIF ELSE i=i+1

      IF i GT (nstruct-1) THEN cont=0
      IF i LT 0 THEN BEGIN 
          print,'Index less than zero. Viewing first object'
          i=0L
      ENDIF 
      delvarx, img, imr, imi, imu
;  ENDFOR 
  ENDWHILE 
  !p.multi=0

END 
