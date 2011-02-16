PRO kappa_noise, nrand, slength, rfac, gridsize, neach, stepfac, noise=noise, write=write, scat=scat, noout=noout

  IF n_params() LT 6 THEN BEGIN 
      print,'-Syntax: kappa_noise, nrand, slength, rfac, gridsize, neach, stepfac, noise=noise, write=write, scat=scat, noout=noout'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Some parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  rmax = slength*rfac
  rdannis, clust, /silent, /noplot

  run1 = 752 & r1str = ntostr(run1)
  run2 = 756 & r2str = ntostr(run2)
  colors = ['u','g','r','i','z']
  clr  = 2

  IF NOT keyword_set(noout) THEN noout=0

  IF n_elements(scat) EQ 0 THEN BEGIN
      indir = '/sdss4/data1/esheldon/CORRECTED/'
      nnm = indir + 'run'+r1str+'_'+r2str+'_srcgal_'+colors[clr]+'_overlap.fit'
      scat = mrdfits(nnm, 1)
  ENDIF 
  time = systime(1)
  COMMON seed, seed

  print,'Doing ',ntostr(nrand),' sets of ',ntostr(neach),' random fields'
  sstr = ntostr(long(slength))
  rstr = ntostr(long(rmax))
  stepstr = ntostr(long(stepfac))

  outdir = '/sdss4/data1/esheldon/GRIDSHEAR/KAPPA/'
  outfile = outdir+'noise_S'+sstr+'_R'+rstr+'_ST'+stepstr+'_N1.dat'
  WHILE exist(outfile) DO BEGIN
      outfile = newname(outfile)
  ENDWHILE 

  print,'Noise file: ',outfile
  
  IF NOT noout THEN BEGIN 
      openw, lun1, outfile, /get_lun
      !textunit = lun1
  ENDIF ELSE print,'  (If you dont set /noout)'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Do the random fields
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nb=0
  FOR jj = 0, nrand-1 DO BEGIN 
      print,'Random #: ',ntostr(jj+1),'/',ntostr(long(nrand))
      check = 1

      ;; Make sure grid is accomodated
      WHILE check DO BEGIN 

          maxd = max(scat.dec)
          mind = min(scat.dec)
          maxr = max(scat.ra)
          minr = min(scat.ra)

          tmp = replicate(clust[0], neach)

          max=maxd-rmax/3600. - 200./3600.
          min=mind+rmax/3600. + 200./3600.
          tmp.dec=arrscl(randomu(seed,neach), min, max, arrmin=0.,arrmax=1.)

          max=maxr-rmax/3600. - 200./3600.
          min=minr+rmax/3600. + 200./3600.
          tmp.ra=arrscl(randomu(seed,neach), min, max, arrmin=0.,arrmax=1.)

          result = where( tag_names(tmp) EQ 'NAME', nw)
          IF nw NE 0 THEN BEGIN
              tmp.name='_rand'
              tmp.z   = 0.
          ENDIF 

          kappa_map, run1, run2, clr, tmp, $
            gridsize=gridsize, slength=slength, rfac=rfac, stepfac=stepfac, $
            scat=scat, $
            write=write, $
            /noprompt, verbose=0, /donoise, $
            noise=noise, radnoise=radnoise, $
            ngood=ngood,$
            check=check
      ENDWHILE 
      fac = sqrt( float(ngood)/neach )
      nb = nb + (1 - fac)
      IF n_elements(outn) EQ 0 THEN BEGIN
          outn=[noise * fac] 
          outr=[radnoise * fac] 
      ENDIF ELSE BEGIN
          outn=[outn,noise*fac]
          outr=[outr,radnoise*fac]
      ENDELSE 
  ENDFOR 
  message = '      noise        radnoise       '
  medn = median(outn)
  medr = median(outr)
  meann = mean_check(outn)
  meanr = mean_check(outr)
  noise = meann*sqrt(neach)
  nbstr = ntostr(long(nb))
  rstr = ntostr(long(nrand))

  IF NOT noout THEN BEGIN 
      printf,lun1,' nstack   : ',neach
      printf,lun1,' Fixed : ',nbstr,'/',rstr,' For accomodation'
      printf,lun1,' grid size: ',gridsize
      printf,lun1,' smoothing length: ',slength
      printf,lun1,' rmax     : ',rmax
      printf,lun1
      printf,lun1,' mean noise: ',meann
      printf,lun1,' sdev      : ',sdev(outn)
      printf,lun1,' mean rad noise: ',meanr
      printf,lun1,' median noise: ',medn
      printf,lun1,' median rad noise: ',medr
      printf,lun1,' mean(noise)*sqrt(n) = ',ntostr(noise)
      printf,lun1
      forprint, outn, outr, TEXT=5, /silent
      close,lun1
      free_lun,lun1
  ENDIF 
  ptime,systime(1)-time

  return
END 
