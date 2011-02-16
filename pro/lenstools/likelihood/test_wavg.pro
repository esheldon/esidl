PRO test_wavg, Nrand, allsame=allsame

  
  zL = 0.1

  sfac = 1.e4

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; generate some source redshifts
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(allsame) THEN BEGIN 
      zs = replicate(0.3, nrand)
  ENDIF ELSE BEGIN 
      ;; draw from some distribution
      zfile = !sdss_shapecorr_dir+'sigmacrit/nzstruct_stripe10.fit'
      nzstruct = mrdfits(zfile,1,/silent)
      genrand, nzstruct.rnz, nzstruct.z, nrand, zs

      plothist,zs,bin=0.01,/norm
      ;;key=get_kbrd(1)
  ENDELSE 
  siginv = sigmacritinv(zl, zs)*sfac

  ;; assign errors

  tsigzs = 0.1
  sigzs = replicate(tsigzs, nrand)*(zs/mean(zs))

  ;; draw from these gaussians
  zsest = fltarr(nrand)
  rand = randomu(seed, nrand,/normal)
  FOR i=0L, nrand-1 DO BEGIN 

      zsest[i] = rand[i]*sigzs[i] + zs[i]

  ENDFOR 

  ;; estimated sigmacritinv based on redshift estimate
  siginvest = sigmacritinv(zl, zsest)*sfac

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; build up redshift histogram based on the estimated errors
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ngrid = 1000
  zsgrid = arrscl(findgen(ngrid), -0.5, 2.0)
  zshist = fltarr(ngrid)

  FOR i=0L, nrand-1 DO BEGIN 
      
      zshist[*] = zshist[*] + gaussprob(zsgrid,zsest[i],sigzs[i])

  ENDFOR 
  zshist = zshist/nrand

  delta = zsgrid[1]-zsgrid[0]
  plot, zsgrid, zshist, color=!red

  wait,2

  genrand, zshist, zsgrid, nrand, randz
  randsiginv = sigmacritinv(zl,randz)*sfac

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; make some plots
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  sisig = sdev(siginvest)
  ninsig = 30
  bin = 2.*sisig/ninsig
  maxplot = nrand

  plothist, siginvest, sigcxhist, sigcyhist, bin=bin
  IF keyword_set(allsame) THEN BEGIN
      oplot,[siginv[0], siginv[0]], [0, maxplot], color=!green, thick=5
  ENDIF ELSE BEGIN 
      plothist, siginv, bin=bin, /overplot, color=!green
  ENDELSE 
  
  plothist, randsiginv, bin=bin, /overplot, color=!magenta

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; now look at mean sigmacrit, from the mean estimates as well
  ;; as the build-up redshift distribution
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  mean_siginv = mean(siginv)
  med_siginvest = median(siginvest)
  mean_siginvest = mean(siginvest)

  IF NOT keyword_set(allsame) THEN BEGIN 
      oplot,[mean_siginv, mean_siginv], [0, maxplot], color=!green, thick=5
  ENDIF 
  oplot, [mean_siginvest, mean_siginvest], [0,maxplot], color=!red
  ;;oplot, [med_siginvest, med_siginvest], [0,maxplot], color=!blue

  mean_randsiginv = mean(randsiginv)
  oplot, [mean_randsiginv, mean_randsiginv], [0,maxplot], color=!magenta

  print,mean(siginv),mean_siginvest, mean_randsiginv

END 
