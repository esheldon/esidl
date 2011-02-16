PRO test_kappa, sigma, slength, rmax, detect, detfrac, s2n, $
                write=write, zlens=zlens, radial=radial


  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: test_kappa, sigma, slength, rmax, detect, detfrac, s2n, write=write, zlens=zlens, radial=radial'
      return
  ENDIF


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Some parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nfiles=130                    ;Number of files to read.

  sstr = ntostr(long(slength))
  rstr = ntostr(long(rmax))
  IF NOT keyword_set(write) THEN write=0
  IF NOT keyword_set(radial) THEN radial=0
  IF n_elements(zlens) EQ 0 THEN zlens=.15

  thresh = [3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7.0]
  nthresh = n_elements(thresh)

  dist = angdist_lambda(zlens)*1000. ;kpc

  nsig = n_elements(sigma)
  a1 = ntostr(long(sigma))
  

  IF radial THEN BEGIN
      nstr = 'rad' 
      estr = 'rerr'
  ENDIF ELSE BEGIN
      nstr='kappa'
      estr = 'kerr'
  ENDELSE 
  dir = '/sdss4/data1/esheldon/CLUSTER/SIM/'
  file=dir+'sim_sig'+a1[0]+'_S'+sstr+'_R'+rstr+'_Z0.15_'+nstr+'_N1.fit'
  t=mrdfits(file, /silent)
  sz = size(t)
  nx = sz[1]
  ny = sz[2]
  cen = [ (float(nx)-1)/2., (float(ny)-1)/2. ]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Define arrays
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  index = lindgen(nx*ny)        ;x and y positions
  x = index MOD nx
  y = index/nx
  
  detect = lonarr( nsig, nthresh )
  posx = fltarr(nfiles) & posx[*] = -1
  posy = posx
  s2ntmp  = posx
  s2n = fltarr(nsig)
  FOR isig=0L, nsig-1 DO BEGIN 
      FOR i=1L, nfiles DO BEGIN 
          
          a2 = ntostr(i)
          kfile=dir+'sim_sig'+a1[isig]+'_S'+sstr+'_R'+rstr+$
                '_Z0.15_'+nstr+'_N'+a2+'.fit'
          nfile=dir+'sim_sig'+a1[isig]+'_S'+sstr+'_R'+rstr+$
                '_Z0.15_'+estr+'_N'+a2+'.fit'
          k=mrdfits(kfile, /silent)
          n=mrdfits(nfile, /silent)

          ksig = k/n

          FOR j=0, nthresh-1 DO BEGIN 
              w=where(ksig GE thresh[j], nw)
              IF nw NE 0 THEN detect[isig, j] = detect[isig, j]+1
          ENDFOR 

          s2ntmp[i-1] = max(ksig)

      ENDFOR 
      s2n[isig] = median(s2ntmp)

  ENDFOR 

  detfrac = detect/float(nfiles)

  IF write THEN BEGIN
      outdir = '/sdss4/data1/esheldon/CLUSTER/'
      name = outdir+nstr+'comp_S'+sstr+'_R'+rstr+'_N1.ps'
      WHILE exist(name) DO BEGIN
          name = newname(name)
      ENDWHILE
      makeps,name,/noland
  ENDIF 

  xtitle = 'S/N threshold'
  ytitle = 'Completeness'

  IF nsig EQ 1 THEN title =nstr+'   Sigma = '+a1+' N = '+ntostr(nfiles)+$
    '  S = '+sstr+' R = '+rstr $
  ELSE title = nstr+'   N = '+ntostr(nfiles)+'  S = '+sstr+'  R = '+rstr

  minx = round( min(thresh)-1 )
  maxx = round( max(thresh)+1 )
  miny = -.1
  maxy = 1.
  plot,[0],[0], xrange = [minx,maxx], xstyle=1, $
    yrange = [miny,maxy], ystyle=1, $
    xtitle=xtitle, ytitle=ytitle, title=title

  psym=1
  FOR isig = 0, nsig-1 DO BEGIN 
      IF psym EQ 3 THEN psym=4
      IF psym EQ 8 THEN psym=1
      IF n_elements(psymkeep) EQ 0 THEN psymkeep=psym $
      ELSE psymkeep = [psymkeep, psym]
      oplot, thresh, detfrac[isig,*], psym=0
      oplot, thresh, detfrac[isig,*], psym=psym
      psym=psym+1
  ENDFOR 
  lmess = replicate('sigma=',nsig) + a1
  legend,lmess,psym=psymkeep,position=[6,.8]

  IF write THEN ep

return
END 
