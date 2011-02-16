PRO mk_ngc_html, new=new

  dir1 = '/sdss7/home/esheldon/ngc/'
  ngc = mrdfits(dir1+'ngc2000.fit',1)
  
  outdir = '/sdss7/home/esheldon/WWW/ngc/'
;  outdir = dir1
  IF keyword_set(new) THEN BEGIN
      mainhtml = outdir+'ngcindex_new.html'
      dir='/sdss7/home/esheldon/ngc/images_new/'
      subdir_string = 'images_new'
  ENDIF ELSE BEGIN 
      mainhtml = outdir+'ngcindex.html'
      dir='/sdss7/home/esheldon/ngc/images/'
      subdir_string = 'images'
  ENDELSE 

  cd,dir
  fullres_files=findfile('*fullres*')
  fchart_files=findfile('*fchart.jpg')

  nf=n_elements(fullres_files)
  IF n_elements(fchart_files) NE nf THEN message,'# of fcharts must match # of fullres'

  ngcnum = intarr(nf)
  name = strarr(nf)
  label = strarr(nf)
  run = strarr(nf)
  camcol = strarr(nf)
  field = strarr(nf)

  FOR i=0L, nf-1 DO BEGIN 
      
      tmp=str_sep(fullres_files[i], '-')
      
      name[i] = tmp[0]
      run[i] = tmp[2]
      camcol[i] = tmp[3]
      field[i] = tmp[4]
      
      label[i] = ' run '+run[i]+ ' camcol '+camcol[i]+ ' field '+field[i]

      num = fix( (str_sep(name[i], 'NGC'))[1] )
      ngcnum[i] = num

  ENDFOR 

  w0=where(ngcnum LT 999,nw)
  w1=where(ngcnum GE 1000 AND ngcnum LT 1999)
  w2=where(ngcnum GE 2000 AND ngcnum LT 2999)
  w3=where(ngcnum GE 3000 AND ngcnum LT 3999)
  w4=where(ngcnum GE 4000 AND ngcnum LT 4999)
  w5=where(ngcnum GE 5000 AND ngcnum LT 5999)
  w6=where(ngcnum GE 6000 AND ngcnum LT 6999)
  w7=where(ngcnum GE 7000 AND ngcnum LT 7999)

  openw, mainlun, mainhtml, /get_lun
  
  printf, mainlun, '<HTML>'
  printf, mainlun, '<HEAD><TITLE>SDSS Images of NGC Galaxies</TITLE></HEAD>'
  printf, mainlun, '<BODY link="#0066ff" vlink="#FF0000" bgcolor="#000000" text="#ffffff">'
  printf, mainlun, '<H2>SDSS Images of NGC Galaxies</H2>'
  printf, mainlun, '<hr size=1><br>'
  printf, mainlun, 'This page contains color images for a subset of the NGC catalog, generated from SDSS "Atlas" images, using code from the <a href="http://sdss4.physics.lsa.umich.edu:8080/~esheldon/weaklens/sdssidl/umich_idl.html">UMich SDSS IDL Libraries</a>. The atlas images and object catalogs were produced by the SDSS data pipeline PHOTO. The images are 6.7X6.7 arcminutes and are in .jpg format.<br><br>'
  printf, mainlun, '<li>All images have the same range. This allows one to directly compare brightness and color. The drawback is that some very bright galaxies may appear "washed out"</li><br><br>'
  printf, mainlun, "<li>The ra,dec of each NGC galaxy was searched for in the SDSS imaging data on the Cheops cluster at UChicago. Because the data is split into different SDSS runs, and there is positional error in the ra,dec of each galaxy, the search code sometimes picks a run adjacent to the correct run. In that case the NGC galaxy may be missing entirely from the image or it may be only partially visible on the left or right edge of the image. (e.g. NGC3521)</li><br><br>"
  printf, mainlun, "<li>Due to positional errors, the code might get the right run but the NGC galaxy may still be missing or only visible on the upper or lower edge of the image.</li><br><br>"
  printf, mainlun, "<li>Some of reconstructions of very large, complex galaxies are not perfect, due to a number of possibilities. e.g. 1) The process of deblending was imperfect. 2) Objects may cover two different fields, and the light may not be detected as an object in the adjacent field. No atlas images exists for that light and the reconstruction may be incomplete.  3) The maximum allowed size for the atlas images was 1000X1000 pixels; if they atlas image is larger than this limit, the image is clipped.</li><br><br>"
  printf, mainlun, '<hr size=1><br>'
  printf, mainlun, "The galaxies are grouped by NGC number. Because some runs overlap, there may be multiple images of some galaxies. The run/camcol/field containing the ra,dec are given. Enjoy!<br><br>"

  ntot=0L
  FOR i=0, 7 DO BEGIN 
      
      min = i*1000
      max = (i+1)*1000
      w=where(ngcnum GE min AND ngcnum LT max, nw)
      ntot = ntot+nw
      IF nw NE 0 THEN BEGIN 
          minstr = ntostr(min)
          maxstr = ntostr(max)
          IF keyword_set(new) THEN BEGIN 
              subhtml = 'NGC'+minstr+'-'+maxstr+'_new.html'
          ENDIF ELSE BEGIN 
              subhtml = 'NGC'+minstr+'-'+maxstr+'.html'
          ENDELSE 
          printf, mainlun, '<a href="./'+subhtml+'">NGC'+minstr+'-'+maxstr+'</a><br>'

          openw, sublun, outdir+subhtml, /get_lun
          printf, sublun, '<HTML>'
          printf, sublun, '<HEAD><TITLE>SDSS Images of NGC Galaxies</TITLE></HEAD>'
          printf, sublun, '<BODY link="#0066ff" vlink="#FF0000" bgcolor="#000000" text="#ffffff">'
          printf,sublun
          printf,sublun,'<TABLE border celpadding=10>'
          printf,sublun,'<TR><TH> <FONT color="#33FF00">Name</FONT></TH>'+ $
                            '<TH> <FONT color="#33FF00">Run</FONT></TH>'+ $
                            '<TH> <FONT color="#33FF00">Camcol</FONT></TH>'+ $
                            '<TH> <FONT color="#33FF00">Field</FONT></TH>'+ $
                            '<TH> <FONT color="#33FF00">RA(J2000)</FONT></TH>'+ $
                            '<TH> <FONT color="#33FF00">DEC(J2000)</FONT></TH>'+ $
                            '<TH> <FONT color="#33FF00">images</FONT></TH>'+ $
                        '</TR>'
          FOR j=0, nw-1 DO BEGIN 

              ww=where(ngc.number EQ ngcnum[w[j]] AND ngc.catalog EQ 'NGC', nww)
              radecstr,ngc[ww].ra,ngc[ww].dec,rastr,decstr
              fullres_url='./'+subdir_string+'/'+fullres_files[w[j]]
              fchart_url='./'+subdir_string+'/'+fchart_files[w[j]]
              printf, sublun, '<TR><TD nowrap><FONT color="yellow">'+name[w[j]]+'</FONT></TD>'+$
                                  '<TD nowrap><FONT color="#00FFFF">'+run[w[j]]+'</FONT></TD>'+$
                                  '<TD nowrap><FONT color="#00FFFF">'+camcol[w[j]]+'</FONT></TD>'+$
                                  '<TD nowrap><FONT color="#00FFFF">'+field[w[j]]+'</FONT></TD>'+$
                                  '<TD nowrap>'+rastr      +                     '</TD>'+$
                                  '<TD nowrap>'+decstr     +                     '</TD>'+$
                                  '<TD nowrap><a href="'+fchart_url+'">8-bit fchart</a>&nbsp&nbsp<a href="'+fullres_url+'">24-bit jpeg</a></TD>'+$
                              '</TR>'

          ENDFOR 
          printf, sublun, '</TABLE>'
          printf, sublun, '</BODY>'
          printf, sublun, '</HTML>'
          free_lun, sublun
      ENDIF 
  ENDFOR 
  
  printf, mainlun, '</BODY>'
  printf, mainlun, '</HTML>'
  free_lun, mainlun
  
  print
  print,'Processed '+ntostr(ntot)+' galaxies'

return
END 
