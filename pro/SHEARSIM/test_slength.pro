PRO test_slength, rfac, write=write, write_kappa, radial=radial

  IF n_params() EQ 0 THEN BEGIN 
      print,'-Syntax: test_slength, rfac, write=write'
      return
  ENDIF 
  
  IF NOT keyword_set(write) THEN write=0
  IF NOT keyword_set(radial) THEN radial=0

  IF NOT radial THEN nstr='kappa' ELSE nstr = 'rad'

  sigma=[0,900,1000,1100,1200]
  slength=[120,140,160,180,200, 240]
  rmax = rfac*slength

  thresh = [3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7.0]
  nthresh = n_elements(thresh)

  nsig = n_elements(sigma)
  nslength = n_elements(slength)

  intfrac = fltarr(nthresh, nslength)

  FOR isl=0, nslength-1 DO BEGIN 

      test_kappa, sigma, slength[isl], rmax[isl], det, detfrac, $
        write=write_kappa, radial=radial
      
      ; integrate over cluster sizes (not including zero)
      FOR ith = 0, nthresh-1 DO BEGIN 
          intfunc, detfrac[1:nsig-1, ith], sigma[1:nsig-1], intdetfrac
          intfrac[ith, isl] = intdetfrac[nsig-2]
      ENDFOR 

  ENDFOR 
 
  rfstr = ntostr(long(rfac))
  dir = '/sdss4/data1/esheldon/CLUSTER/SIM/'
  IF write THEN BEGIN
      outdir = '/sdss4/data1/esheldon/CLUSTER/'
      name = outdir+nstr+'comp_slength_Rfac'+rfstr+'_N1.ps'
      WHILE exist(name) DO BEGIN
          name = newname(name)
      ENDWHILE
      makeps,name,/noland
  ENDIF 
   
  maxx = 1.1*max(slength)
  minx = .9*min(slength)
  maxy = 1.0
  miny = -.1

  plot,[minx],[miny],yrange=[miny,maxy],ystyle=1,xrange=[minx,maxx],xstyle=1, $
    xtitle='Slength',$
    ytitle='Integrated Completeness',$
    title=nstr+'    R = '+rfstr+'*Slength'

  psym=1
  maxint = sigma[nsig-1] - sigma[1]
  FOR ith=0, nthresh-1 DO BEGIN 
      IF psym EQ 3 THEN psym=4
      IF psym EQ 8 THEN psym=1
      IF n_elements(psymkeep) EQ 0 THEN psymkeep=psym $
      ELSE psymkeep = [psymkeep, psym]
      oplot, slength, intfrac[ith, *]/maxint, psym=0
      oplot, slength, intfrac[ith, *]/maxint, psym=psym
      psym=psym+1
  ENDFOR 
  lmess = replicate('thresh=',nthresh) + ntostr(thresh,4)
  pos = [1.1*slength[nslength-2],.95]
  legend,lmess,psym=psymkeep,/right,/top

  IF write THEN ep

  return
END 
