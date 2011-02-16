PRO snefchart_print_header, lun, title

  printf,lun,'<html>'
  printf,lun,'<head>'
  printf,lun,'  <title>'+title+'</title>'
  printf,lun,'</head>'
  
  printf,lun,'<BODY bgcolor="#ffffff" link="#0066ff" vlink="#FF0000" text="#000000">'

END 

PRO snefchart_print_footer, lun

  printf,lun,'  <hr>'
  printf,lun,'  <b>Email</b>: esheldon at cfcp.uchicago.edu<br>'
  printf,lun,'  Last modified: '+systime()
  printf,lun,'</body>'
  printf,lun,'</html>'

END 

PRO snefchart_names, spobjname, objname, $
                     family_psfile, png_atlasfile,$
                     fchartfile, fchartzoomfile, $
                     rawjpgfile, rawzoomjpgfile, $
                     fchartpsfile,$
                     specfile, specpsfile, htmlfile, html_multirun_file, $
                     html_multirun_AtlasFile, $
                     html_multirun_FchartFile, $
                     html_multirun_FchartZoomFile, $
                     multi=multi

  IF keyword_set(multi) THEN dir='multi_obs/' ELSE dir=''

  front = dir + 'sneCand-'+spobjname+'-run'+objname
  family_psfile  = front + '_Family.ps'
  png_atlasfile  = front + '_Atlas.png'
  
  fchartfile     = front + '_Fchart.jpg'
  fchartzoomfile = front + '_FchartZoom.jpg'
  
  rawjpgfile     = front + '_Image.jpg'
  rawzoomjpgfile = front + '_ImageZoom.jpg'
  
  fchartpsfile   = front + '_Fchart.ps'
  
  specfile       = front + '_Spectra.png'
  specpsfile     = front + '_Spectra.ps'
  
  htmlfile       = front + '.html'
  html_multirun_file = front + '_multirun.html'
  html_multirun_AtlasFile = front + '_multirun_Atlas.html'
  html_multirun_FchartFile = front + '_multirun_Fchart.html'
  html_multirun_FchartZoomFile = front + '_multirun_FchartZoom.html'


;  family_psfile = dir + $
;    'sneFamily-'+spobjname+'-run'+objname+'.ps'
;  png_atlasfile = dir + $
;    'sneAtlas-' +spobjname+'-run'+objname+'.png'
  
;  fchartfile  = dir + $
;    'sneFchart-'+spobjname+'-run'+objname+'.jpg'
;  fchartzoomfile  = dir + $
;    'sneFchartZoom-'  +spobjname+'-run'+objname+'.jpg'
  
;  rawjpgfile     = dir + $
;    'sneImage-'+spobjname+'-run'+objname+'.jpg'
;  rawzoomjpgfile = dir + $
;    'sneImageZoom-'  +spobjname+'-run'+objname+'.jpg'
  
;  fchartpsfile    = dir + $
;    'sneFchart-'+spobjname+'-run'+objname+'.ps'
  
;  specfile    = dir + $
;    'sneSpectra-'+spobjname+'-run'+objname+'.png'
;  specpsfile      = dir + $
;    'sneSpectra-'+spobjname+'-run'+objname+'.ps'
  
;  htmlfile = dir + $
;    'sneCand-'+spobjname+'-run'+objname+'.html'
;  html_multirun_file = dir + $
;    'sneCand-'+spobjname+'-run'+objname+'_multirun.html'


END 

PRO snefchart_makehtml, spobjname, $
                        htmlfile, html_multirun_file, $
                        specpsfile, specfile, $
                        fchartpsfile, $
                        fchartfile, fchartzoomfile, $
                        png_atlasfile, table_entry, $
                        std_table_entries=std_table_entries

  ;; Only create the html file if we were
  ;; successful

  print
  print,'Opening HTML file: ',htmlfile
  openw, lun, htmlfile, /get_lun
  printf,lun,'<html>'
  printf,lun,'<head>'
  printf,lun,'  <title>SNe Candidate '+spobjname+'</title>'
  printf,lun,'</head>'
  
  printf,lun,'<BODY bgcolor="#ffffff" link="#0066ff" vlink="#FF0000" text="#000000">'
  printf,lun,'  <h1>Finding charts and spectra for SNe candidate '+spobjname+'</h1>'
  printf,lun,'  The SN fit page is <a href="../CANDpage.php?cand='+spobjname+'">here</a>'
  printf,lun,'  <p>'
  IF fexist(html_multirun_file) THEN printf,lun,'  <h2>There are multiple observations of this object described <a href="./'+html_multirun_file+'">here</a></h2>'

  ;; Some photometric data for this candidate
  printf,lun,'<table border =" 1">'
  printf,lun,'<CAPTION><EM>Photometric Informaton for Candidate</EM></CAPTION>'
  printf,lun,'  <tr>'
  printf,lun,'    <th>Obj name</th><th>Observation Date</th><th>fiber r mag</th><th>magerr</th><th>Photo Vers.</th>'
  printf,lun,'  </tr>'
  printf,lun,table_entry
  printf,lun,'</table>'

  ;; Standard stars from the plate
  printf,lun,'<p>'
  printf,lun,'<table border =" 1">'
  printf,lun,'<CAPTION><EM>Standard Stars on this Plate</EM></CAPTION>'
  printf,lun,'  <tr>'
  printf,lun,'    <th>Standard Star Name</th><th>RA</th><th>Dec</th><th>fiber r mag</th><th>Separation (arcmin)</th>'
  printf,lun,'  </tr>'
  nstd = n_elements(std_table_entries)
  IF nstd GT 0 THEN BEGIN 
      FOR i=0L, nstd-1 DO printf, lun, std_table_entries[i]
  ENDIF 

  printf,lun,'</table>'


  printf,lun,'  <p>'
  
  printf,lun,'  <h2>The spectrum.</h2><br>'
  printf,lun,'  A postscript version is <a href="'+specpsfile+'">here</a><br>'
  printf,lun,'  <img src="'+specfile+'">'
  
  printf,lun,'  <p>'
  printf,lun,'  <h2>A finding chart with the object circled.<br></h2>'
  printf,lun,'  If matched by ra/dec, circle is nearest object, corresponding to the quoted magnitude<br>'
  printf,lun,'  A greyscale postscript version is <a href="'+fchartpsfile+'">here</a><br>'
  printf,lun,'  <img src="'+fchartfile+'">'
  printf,lun,'  <p>'
  printf,lun,'  <h2>Zoomed in a bit. Circle is the size of an SDSS fiber<br></h2>'
  printf,lun,'  If matched by ra/dec, circle is nearest object, corresponding to the quoted magnitude<br>'
  printf,lun,'  <img src="'+fchartzoomfile+'">'
  printf,lun,'  <p>'
  printf,lun,'  <h2>Atlas Images<br></h2>'
  printf,lun,'  <img src="'+png_atlasfile+'">'
  
  printf,lun,'  <hr>'
  printf,lun,'  <b>Email</b>: esheldon at cfcp.uchicago.edu<br>'
  printf,lun,'  Last modified: '+systime()
  printf,lun,'</body>'
  printf,lun,'</html>'
  free_lun,lun

END 

PRO snefchart_makeimages, struct, photoid, markStruct, $
                          fchartfile, rawjpgfile, $
                          fchartzoomfile, rawzoomjpgfile, $
                          family_psfile, png_atlasfile, $
                          extra_markStruct=extra_markStruct, $
                          fchartpsfile, $
                          overwrite=overwrite

  COMMON snefchart_block, wr, radarcmin1, radarcmin2, fibrad

  files_exist = $
    fexist(fchartfile) + $
    fexist(rawjpgfile) + $
    fexist(fchartzoomfile) + $
    fexist(rawzoomjpgfile) + $
    fexist(family_psfile) + $
    fexist(png_atlasfile) + $
    fexist(fchartpsfile)

  IF keyword_set(overwrite) OR (files_exist NE 7) THEN BEGIN 

      maguse = 'fibercounts'

      IF n_elements(extra_markstruct) NE 0 THEN BEGIN 
          nextra_mark=n_elements(extra_markstruct.ra)
          IF nextra_mark GT 1 THEN message,'Only want one extra mark'
      ENDIF 

      ;; default is arcminutes
      rgbfchart, $
        struct[photoid].run, $
        struct[photoid].rerun, $
        struct[photoid].camcol, $
        struct[photoid].field, $
        struct[photoid].id, $
        markStruct=markStruct, $
        extra_markStruct=extra_markStruct, $
        fr=fr, objx=objx, objy=objy, $
        jpegfchart=fchartfile, jpegfile=rawjpgfile, $
        radius=radarcmin1, /directions, $
        /nodisplay,maguse=maguse, order=1
      
      rgbfchart, $
        struct[photoid].run, $
        struct[photoid].rerun, $
        struct[photoid].camcol, $
        struct[photoid].field, $
        struct[photoid].id, $
        markStruct=markStruct, $
        extra_markStruct=extra_markStruct, $
        jpegfchart=fchartzoomfile, jpegfile=rawzoomjpgfile, $
        expand=2, $
        radius=radarcmin2, /directions, $
        /nodisplay,maguse=maguse, order=1
  
      begplot, name=family_psfile
      get_family, struct, photoid, parent, children, siblings, $
        clr=2, pmulti=[0,2,4],maguse=maguse
      endplot 
      
      nsib = n_elements(siblings)
      nsibstr = ntostr(nsib)
      setupplot,'Z'
      
      IF parent[0] EQ -1 THEN BEGIN 
          ;; This is either the parent or an orphan
          ;; Just plot it
          
          device, set_resolution=[735,265]
          get_atlas, struct, photoid,/silent,$
            maguse=maguse
          
      ENDIF ELSE IF nsib EQ 2 THEN BEGIN 
          
          ;; Just this object and sibling
          
          device, set_resolution=[735,675]
          
          !p.multi=[0,0,3]
          
          get_atlas, struct, parent,/silent,maguse=maguse
          legend,'Parent',/right,box=0
          
          w=where(siblings NE photoid)
          get_atlas, struct, siblings[w],/silent,maguse=maguse
          legend,'Sibling',/right,box=0
          
          get_atlas, struct, photoid,/silent,maguse=maguse
          legend,'Candidate',/right,box=0
          !p.multi=0
      ENDIF ELSE IF nsib GT 2 THEN BEGIN 
          ;; Find the closest sibling
          
          device, set_resolution=[735,675]
          
          w=where(struct[siblings].field EQ struct[photoid].field AND $
                  struct[siblings].id    NE struct[photoid].id)
          
          w=siblings[w]
          
          dis = $
            (struct[w].colc[2] - struct[photoid].colc[2])^2 + $
            (struct[w].rowc[2] - struct[photoid].rowc[2])^2
          
          w2=where(dis EQ min(dis))
          w=w[w2]
          
          !p.multi=[0,0,3]
          get_atlas, struct, parent,/silent,maguse=maguse
          legend,'Parent',/right,box=0
          
          get_atlas, struct, w[0],/silent,maguse=maguse
          legend,'Closest of '+nsibstr+' Siblings',/right,box=0
          
          get_atlas, struct, photoid,/silent,maguse=maguse
          legend,'Candidate',/right,box=0
          !p.multi=0
      ENDIF ELSE IF nsib EQ 1 THEN BEGIN 
          ;; Faint version of parent (which has
          ;; no atlas image)
          
          device, set_resolution=[735,265]
          
          get_atlas, struct, photoid,/silent,maguse=maguse
      ENDIF ELSE message,'What to do?'
      
      print,'Writing atlas image: ',png_atlasfile
      write_png, png_atlasfile, tvrd()
      
      
;      setupplot,'X'
;      
;      display_fchart, fr, struct[photoid], objx, objy, 2, $
;        /ps, fnameps=fchartpsfile, $
;        /nodisplay, $
;        /directions,maguse=maguse, order=1

  ENDIF 

END 

PRO snefchart_makepages, ra, dec, run, rerun, camcol, field, id,$
                         spobjname, objname, astrans=astrans, $
                         std_table_entries=std_table_entries, $
                         overwrite=overwrite

    COMMON snefchart_block, wr, radarcmin1, radarcmin2, fibrad

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; Is the input one of the available runs/reruns with atlas, etc?
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


    run_status = sdss_runstatus()
    select_struct = {tsobj_exist:'y',fpatlas_exist:'y'}
    rsi = sdss_flag_select(run_status.flags, 'runstatus', select_struct, nrsi) 

    wgood = where( $
        run_status[rsi].run EQ run AND $
        run_status[rsi].rerun EQ rerun,ngood)
    if ngood eq 0 then begin
        print,'  * No with tsObj/fpAtlas matched ra/dec'
        return
    endif
    wgood = rsi[wgood]

    find_radec, ra, dec, frun, fcamcol, ffield, astrans=astrans
    match, [RUN_STATUS[wgood].run], [frun], mst, mfrun
    IF mst[0] EQ -1 THEN BEGIN 
        print,'  * No with tsObj/fpAtlas matched ra/dec'
        return
    ENDIF ELSE BEGIN 
        frun = frun[mfrun]
        fcamcol = fcamcol[mfrun]
        ffield = ffield[mfrun]

        frerun = RUN_STATUS[wgood[mst]].rerun
    ENDELSE 

  nfound = n_elements(frun)
  fid = replicate(-1,nfound)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Check for the input requested run in our matches
  ;; If not found, we match to one of the found runs by
  ;; RA/DEC
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  wsamerun=where(frun EQ run, nsamerun)
  IF nsamerun NE 0 AND ngood NE 0 THEN BEGIN 
      ;; Replace this occurrence with the input 

      frun[wsamerun] = run
      frerun[wsamerun] = rerun
      fcamcol[wsamerun] = camcol
      ffield[wsamerun] = field

      
      fid[wsamerun] = id

      IF nsamerun GT 1 THEN message,'Duplicate runs here'
  ENDIF ELSE BEGIN 
      print,'--------------------------'
      print,' Matching by RA/DEC'
      print,'--------------------------'
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Make the file names
  ;;;;;;;;;;;;;;;;;;;;;;;;

  snefchart_names, spobjname, objname, $
    family_psfile, png_atlasfile,$
    fchartfile, fchartzoomfile, $
    rawjpgfile, rawzoomjpgfile, $
    fchartpsfile,$
    specfile, specpsfile, $
    htmlfile, html_multirun_file, $
    html_multirun_AtlasFile, $
    html_multirun_FchartFile, $
    html_multirun_FchartZoomFile



  IF nfound GT 1 THEN BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Make a multiple-observations page
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      openw, mlun, html_multirun_file, /get_lun

      snefchart_print_header, mlun, 'SNe Candidate '+spobjname

      printf,mlun,'<h1>Multiple observations for candidate '+spobjname+'</h1>'
      
      printf,mlun,'<p>'
      printf,mlun,'We have '+ntostr(nfound)+$
        ' runs with imaging for this object:'

      printf,mlun,'<table border =" 1">'
      printf,mlun,'<CAPTION><EM>Photometric Informaton for Candidate</EM></CAPTION>'
      printf,mlun,'  <tr>'
      printf,mlun,'    <th>Obj name</th><th>Observation Date</th><th>fiber r mag</th><th>magerr</th><th>Photo Vers.</th>'
      printf,mlun,'  </tr>'

  ENDIF 

  months = ['Jan','Feb','Mar','Apr','May','Jun','Jul',$
            'Aug','Sep','Oct','Nov','Dec']

  FOR i=0L, nfound-1 DO BEGIN 

      fetch_dir, frun[i], fcamcol[i], frerun[i], $
        dir, atldir, tsObjFile, field=ffield[i]
;      struct = mrdfits(tsObjFile,1)
      read_tsobj,[frun[i],frerun[i],fcamcol[i]],struct,$
        start=ffield[i],verbose=0

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Sometimes there are more fields in asTrans than in tsObj.
      ;; This is one of the stupidest things I've ever encountered.
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      IF n_elements(struct) NE 0 THEN BEGIN 

          ;; The style for marking the object.  All default settings, 
          ;; except using red

          markStruct = {color:"red"}

          ;; The id is input if we have the run/rerun of the object on disk
          ;; Otherwise, we match up
          IF fid[i] NE -1 THEN BEGIN 
              ;; This is input run/rerun/etc.
              photoid = where(struct.id EQ fid[i],np)
              IF np EQ 0 THEN BEGIN 
                  message,'Object not found!'
              ENDIF  
              photoid=photoid[0]
          ENDIF ELSE BEGIN 
              
              mygcirc, ra, dec, struct.ra, struct.dec, distance
              tmin = min(distance, photoid)
              
              print,'Distance = '+ntostr(distance[photoid]*3600.)+' arcsec'
              
              IF photoid[0] EQ -1 THEN BEGIN
                  message,'No matches found'
              ENDIF 

              ;; We will circle the input ra/dec separately, since the
              ;; object by run/rerun/... may be at a different position

              extra_markStruct = {ra:ra, dec:dec, radius:fibrad, type:"circle"}

          ENDELSE 
          
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; For main object, we just put everything here.  Multirun we
          ;; Put things into the subdirectory
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          tobjname = $
            sdss_objname(frun[i],frerun[i],fcamcol[i],ffield[i],photoid)

          IF i GT 0 THEN BEGIN 
              snefchart_names, spobjname, tobjname, $
                family_psfile, png_atlasfile,$
                fchartfile, fchartzoomfile, $
                rawjpgfile, rawzoomjpgfile, $
                fchartpsfile,$
                tspecfile, tspecpsfile, $
                thtmlfile, /multi
          ENDIF

          snefchart_makeimages, struct, photoid, markStruct, $
            fchartfile, rawjpgfile, $
            fchartzoomfile, rawzoomjpgfile, $
            family_psfile, png_atlasfile, $
            fchartpsfile, overwrite=overwrite, $
            extra_markstruct=extra_markstruct

          add_arrval, png_atlasFile, atlasFiles
          add_arrval, fchartFile, fchartFiles
          add_arrval, fchartZoomFile, fchartZoomFiles
          add_arrval, fchartPsFile, fchartPsFiles
          add_arrval, tobjname, objNames
          
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Table of photometric info for candidate
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          hdr0=headfits(tsObjFile)
          hdr1=headfits(tsObjFile,ext=1)
          
          mjd_r    = sxpar(hdr1, 'mjd_r')
          phot_ver = sxpar(hdr0, 'PHOT_VER')
          daycnv,(mjd_r)+(2400000.5d), yr, mn, day, hr
          date_string = $
            strn(day,len=2,padchar='0') + '-'+$
            months[mn-1]+'-'+strn(yr)
          
          rmagstr = ntostr(struct[photoid].fibercounts[2])
          rerrstr = ntostr(struct[photoid].fibercountserr[2])
          
          table_entry = $
            '    '+$
            '<td>'+tobjname+'</td>'+$
            '<td>'+date_string+'</td>'+$
            '<td>'+rmagstr+'</td>'+$
            '<td>'+rerrstr+'</td>'+$
            '<td>'+phot_ver+'</td>'
          
          add_arrval, table_entry, table_entries
          
          IF nfound GT 1 THEN BEGIN
              
              printf, mlun, '  <tr>'
              printf, mlun, $
                table_entry
              printf,mlun,'  </tr>'
              
          ENDIF 

      ENDIF 

  ENDFOR 

  ;; same stupid bug with asTrans
  nFound = n_elements(atlasFiles)

  IF nfound GT 1 THEN BEGIN 

      printf,mlun,'</table>'
      printf,mlun,'<p>'

      ttitle = 'Multiple Observations of SNe Candidate '+spObjName

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; The atlas images
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ;; Page with all images
      openw, alun, html_multirun_AtlasFile, /get_lun
      snefchart_print_header, alun, ttitle
      printf, alun, '<h1>'+ttitle + ': Atlas images</h1>'
      FOR i=0L, nfound-1 DO BEGIN 
          printf,alun,'<p><img src="'+atlasFiles[i]+'"><br>'
      ENDFOR 
      snefchart_print_footer, alun
      free_lun, alun

      printf,mlun,'<h2>Atlas Images for each run</h2>'
      printf,mlun,$
        'All images on <a href="'+html_multirun_AtlasFile+'">one page</a>'
      printf,mlun,'<p>'
      printf,mlun,'Individually:<br>'

      FOR i=0L, nfound-1 DO BEGIN 
          printf,mlun,'<a href="'+atlasFiles[i]+'">'+objNames[i]+'</a><br>'
      ENDFOR 
      printf,mlun,'<p>'

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; The finding charts
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ;; Page with all images
      openw, flun, html_multirun_FchartFile, /get_lun
      snefchart_print_header, flun, ttitle
      printf, alun, '<h1>'+ttitle + ': Finding Charts</h1>'
      FOR i=0L, nfound-1 DO BEGIN 
          printf,flun,'<p><img src="'+fchartFiles[i]+'"><br>'
      ENDFOR 
      snefchart_print_footer, flun
      free_lun, flun

      printf,mlun,'<h2>Finding Charts for each run</h2>'
      printf,mlun,'  If matched by ra/dec, dotted circle is nearest object, corresponding to the quoted magnitude<br>'
      printf,mlun,$
        'All images on <a href="'+html_multirun_FchartFile+'">one page</a>'
      printf,mlun,'<p>'
      printf,mlun,'Individually:<br>'

      FOR i=0L, nfound-1 DO BEGIN 
          printf,mlun,'<a href="'+fchartFiles[i]+'">'+objNames[i]+'</a><br>'
      ENDFOR 
      printf,mlun,'<p>'


      ;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Zoomed finding charts
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ;; Page with all images
      openw, fzlun, html_multirun_FchartZoomFile, /get_lun
      snefchart_print_header, fzlun, ttitle
      printf, alun, '<h1>'+ttitle + ': Zoomed Finding Charts</h1>'
      FOR i=0L, nfound-1 DO BEGIN 
          printf,fzlun,'<p><img src="'+fchartZoomFiles[i]+'"><br>'
      ENDFOR 
      snefchart_print_footer, fzlun
      free_lun, fzlun

      printf,mlun,'<h2>Zoomed in finding charts.  Circle is fiber size</h2>'
      printf,mlun,'  If matched by ra/dec, dotted circle is nearest object, corresponding to the quoted magnitude<br>'
      printf,mlun,$
        'All images on <a href="'+html_multirun_FchartZoomFile+'">one page</a>'
      printf,mlun,'<p>'
      printf,mlun,'Individually:<br>'
      FOR i=0L, nfound-1 DO BEGIN 
          printf,mlun,'<a href="'+fchartZoomFiles[i]+'">'+objNames[i]+'</a><br>'
      ENDFOR 
      printf,mlun,'<p>'



      ;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Greyscale
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;

      printf,mlun,'<h2>Greyscale postscript</h2>'
      FOR i=0L, nfound-1 DO BEGIN 
          printf,mlun,'<a href="'+fchartPsFiles[i]+'">'+objNames[i]+'</a><br>'
      ENDFOR 


      snefchart_print_footer, mlun
      free_lun,mlun
  ENDIF 

  ;; If multiple observations, make the multirun page

  IF nsamerun EQ 1 THEN table_entry = table_entries[wsamerun] $
  ELSE table_entry = table_entries[0]

  ;; Write out the html page for the first object found
  snefchart_makehtml, spobjname, $
    htmlfile, html_multirun_file, $
    specpsfile, specfile, $
    fchartpsfile, $
    fchartfile, fchartzoomfile, $
    png_atlasfile, table_entry, $
    std_table_entries=std_table_entries


END 





PRO snefchart, html_overwrite=html_overwrite, overwrite=overwrite, retry_fchart=retry_fchart, snefile=snefile

  COMMON snefchart_block, wr, radarcmin1, radarcmin2, fibrad


  ;; Often we didn't have the runs to make the finding chart
  ;; In that case, no html file is written, so that the web
  ;; page won't try to show a link to that file
  ;; However, we don't want to check every time, so something
  ;; must be created so that snefchart won't try to look for runs
  ;; again **unless requested to do so**
  ;; So we will always output the spectrum, and only try to make
  ;; the finding chart when that file exists if /retry_fchart is
  ;; set or we are overwriting all

  run_status = sdss_runstatus()
  rs = {tsobj_exist: 'y', fpatlas_exist: 'y'}
  wr = sdss_flag_select(run_status.flags, 'runstatus', rs)

  run_status = sdss_runstatus()
  if wr[0] eq -1 then begin
      print,'  * No good runs found for making  fcharts'
      return
  endif else begin 
      ;; Get only latest reruns
      rmd = rem_dup(run_status[wr].run)
      urun = run_status[wr[rmd]].run
      nrun = n_elements(urun)
      FOR i=0L, nrun-1 DO BEGIN 

          wt=where(run_status[wr].run EQ urun[i])

          maxrerun = max(run_status[wr[wt]].rerun, wt2)

          wt = wr[wt[wt2]]

          add_arrval, wt, wr2

      endfor 
      wr = wr2

  endelse 

  make_tsflag_struct, ts
  ts.spectrophoto_std = 'Y'
  
  ;; All in arcminutes
  radarcmin1 = 2
  radarcmin2 = 1
  fibrad = 7.5/2.0*0.4/60.0

  ;; read in coordinates for supernova candidates and make a finding chart and
  ;; atlas images

  sne_dir = sdssidl_config('spec_dir') + '1d_23/sne/'
  outdir = esheldon_config("www_dir") + 'SN/snefchart/'
;  outdir = '~/tmp/snefchart/'
  cd,outdir
  sne_dir = '/net/cheops4/home/garyk/sne/newcode/'
  coord_file = sne_dir + 'allCoords.dat'

  ;; read coordinates here
;  readcol, coord_file, $
;    fermifile, $
;    cra, cdec, $
;    cplate, cfiber, cmjd, $
;    crun, crerun, ccamcol, cfield, cid, $
;    format='A,D,D,L,L,L,L,L,L,L,L'

;  readcol, coord_file, cra, cdec, crun, crerun, ccamcol, cfield, cid, $
;    format='D,D,L,L,L,L,L'

  read_snecand, snestruct, snefile=snefile
  
  nobj = n_elements(snestruct)
  objnames = sdss_objname(snestruct.run,$
                          snestruct.rerun,$
                          snestruct.camcol,$
                          snestruct.field,$
                          snestruct.id)

  spobjnames = strarr(nobj)

  FOR i=0L, nobj-1 DO BEGIN 
      plate = snestruct[i].plate
      fiber = snestruct[i].fiber
      mjd   = snestruct[i].mjd

      spobjnames[i] = $
        strn(mjd)+'-'+$
        strn(plate,len=4,padchar='0')+'-'+$
        strn(fiber,len=3,padchar='0')
  ENDFOR 

  tot_objnames = spobjnames + '-run'+objnames

  rmd = rem_dup(tot_objnames)
  nobj = n_elements(rmd)

  print

  colprint,tot_objnames[rmd]
  print,'Nobj = '+ntostr(nobj)

  plateold = -1L
  FOR ii=0L, nobj-1 DO BEGIN 
;  FOR ii=1L, 1 DO BEGIN 
      i = rmd[ii]

      plate = snestruct[i].plate
      fiber = snestruct[i].fiber
      mjd   = snestruct[i].mjd

      spobjname = $
        strn(mjd)+'-'+$
        strn(plate,len=4,padchar='0')+'-'+$
        strn(fiber,len=3,padchar='0')
      
      snefchart_names, spobjname, objnames[i], $
        family_psfile, png_atlasfile,$
        fchartfile, fchartzoomfile, $
        rawjpgfile, rawzoomjpgfile, $
        fchartpsfile,$
        specfile, specpsfile, $
        htmlfile, html_multirun_file
      
      print
      print,'-----------------------------------------'
      print,'  Cand: '+spobjname
      print,'-----------------------------------------'

      ;; If specfile is there, we have already looked at this object.
      ;; By default, we will skip.
      ;; 
      ;; Other possibilities: 
      ;;    1) try to remake the finding chart if its not ther
      ;;    2) remake the html, even if its there.
      ;;    3) redo everything

      IF ( (NOT fexist(specfile)) OR $
           keyword_set(overwrite) OR $
           keyword_set(html_overwrite) OR $
           keyword_set(retry_fchart)) THEN BEGIN 

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Plot the spectrum, both to the PS and Z devices
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          
          ;; Read the spectrum
          read_spec1d, plate, fiber, spec, mjd=mjd

          begplot,name=specpsfile,/color
          plot_spec1d, spec, nsmooth=10
          endplot 
              
          setupplot, 'Z'
          device, set_resolution=[640,512]
          simpctable, rmap, gmap, bmap
          plot_spec1d, spec, nsmooth=10
          print
          print,'Writing spec png file: ',specfile
          write_png, specfile, tvrd(), rmap, gmap, bmap
          setupplot,'X'

          ;; Still only retrying if it *doesnt' exist* unless /overwrite
          IF ( (NOT fexist(fchartfile)) OR $
               keyword_set(overwrite) OR $
               keyword_set(html_overwrite)) THEN BEGIN 

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; Nearby standard stars
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; Standard Stars in this plate
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              IF plate NE plateold THEN BEGIN 
                  read_spdiag1d, plate, spdiag, mjd=mjd
                  
                  tsflag_select, spdiag, ts, wstd
                  
                  IF wstd[0] NE -1 THEN BEGIN
 
                      nstd = n_elements(wstd)

                      tstd_table_entries = strarr(nstd)
                      std_table_entries = strarr(nstd)
                      FOR jj=0L, nstd-1 DO BEGIN 
                          jstd = wstd[jj]

                          std_objname = $
                            strn(mjd)+'-'+$
                            strn(plate,len=4,padchar='0')+'-'+$
                            strn(spdiag[jstd].fiberid,len=3,padchar='0')

                          ;; will need to augment with distance from
                          ;; object
                          radecstr, $
                            spdiag[jstd].ra, spdiag[jstd].dec,$
                            rastr, decstr
                          tstd_table_entries[jj] = $
                            '    '+$
                            '<td>'+std_objname+'</td>'+$
                            '<td>'+rastr+'</td>'+$
                            '<td>'+decstr+'</td>'+$
                            '<td>'+ntostr(spdiag[jstd].ugrizfibre[2])+'</td>'
                            
                      ENDFOR 
                  ENDIF ELSE BEGIN 
                      delvarx, std_table_entries
                  ENDELSE 
              ENDIF 

              IF wstd[0] NE -1 THEN BEGIN 
                  mygcirc, $
                    spec.raobj, spec.decobj, $
                    spdiag[wstd].ra, spdiag[wstd].dec, $
                    distance

                  distance = distance*60.
                  s = sort(distance)
                  FOR jj=0L, nstd-1 DO BEGIN 
                      std_table_entries[jj] = $
                        '<tr>'+tstd_table_entries[s[jj]] + $
                        '<td>'+ntostr(distance[s[jj]])+'</td></tr>'
                  ENDFOR 
              ENDIF
 
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; Plot the finding chart
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              
              snefchart_makepages, $
                snestruct[i].ra, snestruct[i].dec, $
                snestruct[i].run, snestruct[i].rerun, $
                snestruct[i].camcol, snestruct[i].field, snestruct[i].id,$
                spobjname, objnames[i], astrans=astrans, $
                std_table_entries=std_table_entries, $
                overwrite=overwrite

          ENDIF ;; checked for fchart file or /overwrite
          
      ENDIF ELSE BEGIN ;; overwrite, retry_fchart, or no specfile
          print
          print,'  * Skipping'
      ENDELSE 

  ENDFOR 


  print
  print,'Opening index.html file'
  openw, ilun, 'index.html', /get_lun

  printf,ilun,'<html>'
  printf,ilun,'<head>'
  printf,ilun,'  <title>SDSS Supernova Candidates</title>'
  printf,ilun,'</head>'
  
  printf,ilun,'<BODY bgcolor="#ffffff" link="#0066ff" vlink="#FF0000" text="#000000">'

  FOR ii=0L, nobj-1 DO BEGIN 

      i = rmd[ii]
      
      plate = snestruct[i].plate
      fiber = snestruct[i].fiber
      mjd   = snestruct[i].mjd

      spobjname = $
        strn(mjd)+'-'+$
        strn(plate,len=4,padchar='0')+'-'+$
        strn(fiber,len=3,padchar='0')

      htmlfile = $
        'sneCand-'+spobjname+'-run'+objnames[i]+'.html'

      IF fexist(htmlfile) THEN begin
          printf,ilun,'  <a href="'+htmlfile+'">'+spobjname+'</a><br>'
      ENDIF 

  ENDFOR 

  free_lun, ilun
  print,'Done'

END 
