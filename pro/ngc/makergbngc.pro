PRO makergbngc, astrans=astrans, addu=addu, overwrite=overwrite, sdss=sdss

  indir='/net/cheops1/data0/esheldon/ngc/'
  file=indir+'ngc2000.fit'
  ngc=mrdfits(file,1)

  outdir = '/net/cheops1/data0/esheldon/ngc/'
  nobj = n_elements(ngc)

  radius = 500
  contrast = 45
  
  w=where(ngc.catalog EQ 'NGC' AND ngc.type EQ 'Gx ',nobj)

  s = sort(ngc[w].number)
  w = w[s]

  maxsize = [1000,1000]

  ;; set to Z-buffer (8-bit) to 
  ;; make finding chart
  device_old = !d.name
  setupplot,'Z'
  device, set_resolution=[640,512]*1.2

  FOR i=0L, nobj-1 DO BEGIN 

      ind=w[i]
      print,'>>> Doing '+ngc[ind].name
      ra = ngc[ind].ra
      dec = ngc[ind].dec
      find_radec, ra, dec, run, camcol, field, /silent, astrans=astrans
      IF run[0] NE -1 THEN BEGIN 
          nrun = n_elements(run)
          FOR j=0L, nrun-1 DO BEGIN 

              rstr=run2string(run[j])
              cstr=ntostr(camcol[j])
              fstr=field2string(field[j])
              addname = rstr+'-'+cstr+'-'+fstr
              pngfchart = outdir + ntostr(ngc[ind].name)+'-'+ntostr(ngc[ind].type)$
                +'-'+addname+'-fchart.png'
              jpgfchart = repstr(pngfchart,'png','jpg')
              jpegfile  = outdir + ntostr(ngc[ind].name)+'-'+ntostr(ngc[ind].type)$
                +'-'+addname+'-fullres.jpg'
              name = ngc[ind].name+'  '+ngc[ind].type+'   '+addname

              skip = 0
              jskip = 0
              pskip = 0
              IF NOT keyword_set(overwrite) THEN BEGIN 
                  IF fexist(jpegfile) THEN BEGIN
                      message,'JPEG File '+jpegfile+' already exists; skipping',/inf
                      delvarx, jpegfile
                      jskip = 1
                  ENDIF 
                  IF fexist(pngfchart) THEN BEGIN 
                      message,'PNG fchart '+pngfchart+' already exists; skipping',/inf
                      delvarx, pngfchart
                      pskip=1
                  ENDIF 
                  IF jskip or pskip THEN skip=1
              ENDIF 
              IF NOT skip THEN BEGIN 
                  print,'--- '+name
                  rgbfchart, ngc[ind].ra, ngc[ind].dec, radius=radius, $
                             jpegfile=jpegfile, pngfchart=pngfchart,$
                             title=name,$
                             xtitle='',runuse=j,nodisplay=nodisplay, $
                             maxsize=maxsize, astrans=astrans,$
                             addu=addu, contrast=contrast, status=status, $
                             sdss=sdss
                  IF fexist(pngfchart) AND (NOT pskip) THEN spawn,'convert -verbose '+pngfchart+' '+jpgfchart
              ENDIF 
              
          ENDFOR 
      ENDIF
      IF i MOD 100 EQ 0 THEN print,'.',format='($,a)'

  ENDFOR 

  setupplot,device_old

  return
END 
