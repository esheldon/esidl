PRO large_map, run1, run2, clr, yoverx, slength=slength, scat=scat, noise=noise,check=check

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: large_map, run1, run2, clr, yoverx, slength=slength, scat=scat, noise=noise, check=check'
      return
  ENDIF 

  rdannis,c,/noplot,/silent

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Some parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF NOT keyword_set(check) THEN check=0
  IF NOT keyword_set(noise) THEN BEGIN
      noise=0
      donoise = 0
  ENDIF ELSE BEGIN
      donoise = 1
  ENDELSE 

  IF n_elements(slength) EQ 0 THEN slength = 120. ; Arcseconds ( was 120)
  rfac = 5.                     ; This is trial. (was 20)
  rmax = rfac*slength
  tol = 30./3600.
  colors = ['u','g','r','i','z']
  outdir = '/sdss4/data1/esheldon/LARGEMAPS/tmp/'
  ncent = 100

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Read source catalog
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  r1str = ntostr(run1)
  r2str = ntostr(run2)
  IF n_elements(scat) EQ 0 THEN BEGIN 
      dir='/sdss4/data1/esheldon/CORRECTED/'
      fend='_srcgal_'+colors[clr]+'_overlap.fit'
      file=dir+'run'+r1str+'_'+r2str+fend
      IF NOT exist(file) THEN BEGIN
          file=dir+'run'+r2str+'_'+r1str+fend
          IF NOT exist(file) THEN BEGIN
              print,'No '+colors[clr]+' band overlap file exists for ', $
                    'runs '+r1str+' and '+r2str
              return
          ENDIF 
      ENDIF 

      scat=mrdfits(file, 1, hdr)
  ENDIF 
  ws = where(scat.uncert LT .64, nscat)
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Define map region
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  maxra = max(scat[ws].ra)
  minra = min(scat[ws].ra)
  maxdec = max(scat[ws].dec)
  mindec = min(scat[ws].dec)

  xsize = maxdec - mindec - 2.*tol - 2.*rmax/3600. ;degrees
  ysize = yoverx*xsize

  
  pos_1 = minra + tol + rmax/3600. + ysize/2.
  
  posy = pos_1 + ysize*findgen(ncent)

  ycheck = maxra - tol - rmax/3600. - ysize/2.

  w = where(posy LT ycheck, ncent)
  IF ncent NE 0 THEN posy = posy[w] ELSE BEGIN 
      print,'What!'
      return
  ENDELSE 

  maxx   = 0. + xsize/2.
  minx   = 0. - xsize/2.
  maxy = max(posy) + ysize/2.
  miny = min(posy) - ysize/2.

  left = maxra - ( max(posy) + ysize/2. + rmax/3600.)

  print,'Map Dec range: ',ntostr(minx),' ',ntostr(maxx)
  print,'Map Ra  range: ',ntostr(miny),' ',ntostr(maxy)
  print,'Using ',ntostr(ncent),' Centers'
  print,'Covering ',ntostr( (maxx-minx)*(maxy-miny) ),' square degrees'
  print,'Leftover = ',ntostr(left),' degrees ',ntostr(left*3600.),' arcsec '

  IF check THEN BEGIN 
      reply = ' '
      print,format='($, "Is this Acceptable(y/n)")'
      read, reply
      IF (reply EQ 'n') OR (reply EQ 'N') THEN return
  ENDIF 
  gridsize = xsize*3600.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make 'lens' catalog, centered on our positions.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  radec = replicate(c[0], ncent)
  radec.dec = fltarr(ncent)     ; Dec centered on zero
  radec.ra = posy
  IF NOT noise THEN BEGIN
      radec.name = 'wmap_'+ntostr(radec.ra,8)+'_'+ntostr(radec.dec,6)
  ENDIF ELSE BEGIN
      radec.name = 'wnoise_'+ntostr(radec.ra,8)+'_'+ntostr(radec.dec,6)
  ENDELSE 
  radec.z = 0.

  dir = '/sdss4/data1/esheldon/LARGEMAPS/'
  addstr = colors[clr]+'_N1.fit'
  IF NOT noise THEN BEGIN 
      kname = dir + 'big_map_kappa_'+addstr
      ename = dir + 'big_map_kerr_'+addstr
      dname = dir + 'big_map_dens_'+addstr
  ENDIF ELSE BEGIN 
      kname = dir + 'big_map_noise_'+addstr
      ename = dir + 'big_map_nerr_'+addstr
  ENDELSE 
  WHILE exist(kname) OR exist(ename) DO BEGIN
      kname = newname(kname)
      ename = newname(ename)
      IF n_elements(dname) NE 0 THEN dname = newname(dname)
  ENDWHILE 
  print
  print,'Kappa file: ',kname

  kaplist = strarr(ncent)
  kerrlist = strarr(ncent)
  denslist = strarr(ncent)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make maps around centers
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  fgalfile = '/sdss4/data1/esheldon/CORRECTED/run752_756_fgal_r.fit'
  fgal = mrdfits(fgalfile, 1)

  FOR i=0, ncent-1 DO BEGIN 

      kappa_map, run1, run2, clr, radec[i], $
                 slength=slength, rfac=rfac, yoverx=yoverx, $
                 gridsize=gridsize, outdir=outdir, scat=scat[ws], $
                 donoise=donoise, $
                 /write, /no_ps, $
                 fgal=fgal, $
                 kapfit=kapfit, kerrfit=kerrfit, densfit=densfit
      kaplist[i] = kapfit
      kerrlist[i] = kerrfit
      denslist[i] = densfit
  ENDFOR 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Combine files (faster to do small ones)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  k = mrdfits(kaplist[0], /silent)
  s=size(k)
  sx = s[1]
  sy = s[2]
  ktot = fltarr(sx, ncent*sy)
  etot = ktot
  dtot = ktot

  FOR i=0, ncent-1 DO BEGIN 
      
      k=mrdfits(kaplist[i], /silent)
      e=mrdfits(kerrlist[i], /silent)
      d=mrdfits(denslist[i], /silent)

      ktot[*, i*sy:(i+1)*sy - 1 ] = k
      etot[*, i*sy:(i+1)*sy - 1 ] = e
      dtot[*, i*sy:(i+1)*sy - 1 ] = d
      
  ENDFOR 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output combined files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  mwrfits, ktot, kname, /create
  mwrfits, etot, ename, /create
  IF n_elements(dname) NE 0 THEN mwrfits, dtot, dname, /create
  return 
END 
      
