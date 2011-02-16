
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    PLOTCLRZSHEAR
;       
; PURPOSE:
;    make shear analysis plots for lenses with redshifts in 3 colors.
;    Sorry, this is just a big hodgepodge of stuff I just kept
;    adding without commenting.
;
; CALLING SEQUENCE:
;    
;
; INPUTS: 
;    
;
; OPTIONAL INPUTS:
;    
;
; KEYWORD PARAMETERS:
;    
;       
; OUTPUTS: 
;    
;
; OPTIONAL OUTPUTS:
;    
;
; CALLED ROUTINES:
;    
; 
; PROCEDURE: 
;    
;	
;
; REVISION HISTORY:
;    
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO plot_zshear_old_getwgm, combstruct

  ;; will use r0,err from here as guess
  denscont2wgm, combstruct, wgm, wgmerr, /silent
  denscont2wgm, combstruct, wgm_rebin, wgmerr_rebin, /rebin, /silent

  calc_r0_gamma_wgm, combstruct, $
                     gamma, r0, gammalow, gammahigh, r0low, r0high, $
                     /dolegend

  nr = n_elements(combstruct.meanr)
  nrebin = n_elements(combstruct.meanr_rebin)

  val = 'fltarr('+ntostr(nr)+')'
  val_rebin = 'fltarr('+ntostr(nrebin)+')'

  names = ['wgm', 'wgmerr', 'wgm_rebin', 'wgmerr_rebin', $
           'r0', $
           'r0low', 'r0high', $
           'gamma', $
           'gammalow', 'gammahigh']
  vals = [val, val, val_rebin, val_rebin, $
          '0.0', $
          'fltarr(3)', 'fltarr(3)', $
          '0.0', $
          'fltarr(3)', 'fltarr(3)']

  add_tags, temporary(combstruct), names, vals, combstruct
  combstruct.wgm = wgm
  combstruct.wgmerr = wgmerr
  combstruct.wgm_rebin = wgm_rebin
  combstruct.wgmerr_rebin = wgmerr_rebin
  combstruct.r0 = r0
  combstruct.r0low = r0low
  combstruct.r0high = r0high
  combstruct.gamma = gamma
  combstruct.gammalow = gammalow
  combstruct.gammahigh = gammahigh

END 

PRO plot_zshear_old_meanlum_gmr, lensum, combstruct, lumclr, $
                             tmeanlum, tmeanlumerr, meanlum, meanlumerr, $
                             tmeangmr, tmeangmrerr, meangmr, meangmrerr, $
                             lumfrac, lumfracerr, gmrfrac, gmrfracerr

  ;; absmag cuts
  read_cuts,cuts
  lowcut = cuts.main_minmag_threebin[lumclr,2]
  highcut = cuts.main_maxmag_threebin[lumclr,0]

  psym=8

  nbin = n_elements(lensum[0].rsum)

  meanlum = fltarr(nbin)
  meanlumerr = fltarr(nbin)
  meangmr = fltarr(nbin)
  meangmrerr = fltarr(nbin)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; lum
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Calculating mean lum'

  FOR i=0L, nbin-1 DO BEGIN 

      ww = where(lensum.wsum[i] NE 0.0 AND $
                 lensum.absmag[lumclr] GE lowcut AND $
                 lensum.absmag[lumclr] LE highcut)

      weights = 1./lensum[ww].sigmaerr[i]^2
      wsum = total(weights)
      meanlum[i] = total(weights*lensum[ww].lum[lumclr])/wsum
      lumerr2 = total( weights^2*(lensum[ww].lum[lumclr] - meanlum[i])^2 )
      meanlumerr[i] = sqrt(lumerr2)/wsum
      
  ENDFOR 

  ww = where(lensum.weight NE 0.0 AND $
             lensum.absmag[lumclr] GE lowcut AND $
             lensum.absmag[lumclr] LE highcut)
             
  weights = lensum[ww].weight
  wsum = total(weights)
  tmeanlum = total(weights*lensum[ww].lum[lumclr])/wsum
  tlumerr2 = total( weights^2*(lensum[ww].lum[lumclr] - tmeanlum)^2 )
  tmeanlumerr = sqrt(tlumerr2)/wsum
  
  ;; plot lum

  !p.multi = [0,0,2]
  
  xr = [0.8*min(combstruct.meanr), 1.2*max(combstruct.rmax_act)]
  xr[0] = xr[0] < 10.
  
  yt = 'Mean Lum '+!colors[lumclr]+' [10!U10!N L'+sunsymbol()+']'
  ploterror, combstruct.meanr, meanlum, meanlumerr, psym=psym, $
             xtitle=!kpcxtitle,xrange=xr, /xsty, $
             ytitle=yt, /ynozero, /xlog,xtickf='loglabels'

  oploterror, [combstruct.rmax_act[nbin-1]], [tmeanlum], [tmeanlumerr], $
              psym=psym, color=!blue, errcolor=!blue
  
  oplot, [1, 100000], [tmeanlum, tmeanlum]
  
  lumratio = meanlum/tmeanlum
  lumratioerr = abs(lumratio)*sqrt( $
                                    (meanlumerr/meanlum)^2 + $
                                    (tmeanlumerr/tmeanlum)^2 $
                                  )
  lumfrac = lumratio - 1.
  lumfracerr = lumratioerr
  
  ;;yt = '(L(R) - L)/L'
  yt = 'Fractional Diff'
  ploterror, combstruct.meanr, lumfrac, lumfracerr, psym=psym, $
             xtitle=!kpcxtitle,xrange=xr, /xsty, $
             ytitle=yt, /ynozero, /xlog,xtickf='loglabels'
  oplot, [1, 100000], [0,0]
  
  !p.multi=0


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; g-r
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Calculating mean g-r'

  maxmag = -17.0
  minmag = -23.0

  FOR i=0L, nbin-1 DO BEGIN 

      ww = where(lensum.wsum[i] NE 0.0 AND $
                 lensum.abscounts[1] GT 0.0 AND $
                 lensum.abscounts[2] GT 0.0 AND $
                 lensum.abscounts[1] NE 1000. AND $
                 lensum.abscounts[2] NE 1000.)

      weights = 1./lensum[ww].sigmaerr[i]^2
      wsum = total(weights)
      
      gmr = lensum[ww].abscounts[1] - lensum[ww].abscounts[2]
      meangmr[i] = total(weights*gmr)/wsum
      gmrerr2 = total( weights^2*(gmr - meangmr[i])^2 )
      meangmrerr[i] = sqrt(gmrerr2)/wsum
      
  ENDFOR 

  ww = where(lensum.weight NE 0.0 AND $
             lensum.abscounts[1] GT 0.0 AND $
             lensum.abscounts[2] GT 0.0 AND $
             lensum.abscounts[1] NE 1000. AND $
             lensum.abscounts[2] NE 1000.)
;stop
  weights = lensum[ww].weight
  wsum = total(weights)

  gmr = lensum[ww].abscounts[1] - lensum[ww].abscounts[2]
  tmeangmr = total(weights*gmr)/wsum
  tgmrerr2 = total( weights^2*(gmr - tmeangmr)^2 )
  tmeangmrerr = sqrt(tgmrerr2)/wsum

  ;; plot g-r

  !p.multi = [0,0,2]
  
  xr = [0.8*min(combstruct.meanr), 1.2*max(combstruct.rmax_act)]
  xr[0] = xr[0] < 10.
  
  gmrstr = 'g'+!csym.minus+'r'
  yt = 'Mean '+gmrstr
  ploterror, combstruct.meanr, meangmr, meangmrerr, psym=psym, $
             xtitle=!kpcxtitle,xrange=xr, /xsty, $
             ytitle=yt, /ynozero, /xlog,xtickf='loglabels'

  oploterror, [combstruct.rmax_act[nbin-1]], [tmeangmr], [tmeangmrerr], $
              psym=psym, color=!blue, errcolor=!blue
  
  oplot, [1, 100000], [tmeangmr, tmeangmr]
  
  gmrratio = meangmr/tmeangmr
  gmrratioerr = abs(gmrratio)*sqrt( $
                                    (meangmrerr/meangmr)^2 + $
                                    (tmeangmrerr/tmeangmr)^2 $
                                  )
  gmrfrac = gmrratio - 1.
  gmrfracerr = gmrratioerr

  yt = 'Fractional Diff'
;  yt = '( '+gmrstr+'(R)'+!csym.minus+' '+gmrstr+' )/'+gmrstr
  ploterror, combstruct.meanr, gmrfrac, gmrfracerr, psym=psym, $
             xtitle=!kpcxtitle,xrange=xr, /xsty, $
             ytitle=yt, /ynozero, /xlog,xtickf='loglabels'
  oplot, [1, 100000], [0,0]
  
  !p.multi=0

END 

PRO plot_zshear_old_fitpower, combstruct, bestpow, bestnorm, $
                          powlow, powhigh, normlow, normhigh, $
                          sigallow_low, sigallow_high, combsigdiff_fit, $
                          dojack=dojack, wuse=wuse

  leg_charsize=1.0

  IF keyword_set(dojack) THEN BEGIN 
      errsend = combstruct.covariance
  ENDIF ELSE BEGIN 
      errsend = combstruct.sigmaerr
  ENDELSE 

  print
  print,'Initial power law fitting'
  fitpower, combstruct.meanr[wuse]/1000., $
            combstruct.sigma[wuse], combstruct.sigmaerr[wuse], [5., -.8], $
            tyfit, Aout, Asig

  IF keyword_set(dojack) THEN rangefac = 10.0 ELSE rangefac=5.0
  Asig[1] = Asig[1]*rangefac > 0.3 < 2.0
  Asig[0] = Asig[0]*rangefac > 2.0 < 30.0

  spowrange = Aout[1] + [-Asig[1], +Asig[1]]
  normrange = Aout[0] + [-Asig[0], Asig[0]]

  print
  print,'Power law fits'

  nnorm = 400
  npow = 400

  IF keyword_set(dojack) THEN BEGIN 
      errsend = combstruct.covariance
  ENDIF ELSE BEGIN 
      errsend = combstruct.sigmaerr
  ENDELSE 

  pow_chisq_conf_gen, combstruct.meanr/1000., combstruct.sigma, errsend, $
                      spowrange, normrange, npow, nnorm, $
                      chisq_surf, $
                      bestpow, bestnorm, $
                      powlow, powhigh, $
                      normlow, normhigh, $
                      xtitle='Index', $
                      ytitle='Norm [h M'+sunsymbol()+$
                      ' pc!U'+!csym.minus+'2!N]', $
                      names = ['Index', 'Norm'], $
                      xstyle=1, ystyle=1, charsize=1.5, $
                      aspect=1, /center, $
                      yfit=combsigdiff_fit, $
                      minchisq=minchisq, degfree=degfree, $
                      yallow_low = sigallow_low, $
                      yallow_high=sigallow_high,wuse=wuse,/dolegend, $
                      nkeep = [5, 5]

  clevel = 0

  nhighdiff = normhigh[clevel]-bestnorm & nlowdiff = bestnorm - normlow[clevel]
;  mess = 'Norm = '+ntostr(bestnorm,5,/round)+$
;    '!S!U'+!csym.plus+ntostr(nhighdiff,5,/round)+'!R!D'+!csym.minus+ntostr(nlowdiff,5,/round)
;  mess = [mess, '']
  phighdiff = powhigh[clevel]-bestpow & plowdiff = bestpow - powlow[clevel]
;  mess = [mess, $
;          'Index = '+ntostr(bestpow, 5,/round) + $
;          '!S!U'+!csym.plus+ntostr(phighdiff,5,/round)+'!R!D'+!csym.minus+ntostr(plowdiff,5,/round)]
;  legend, mess, /right, charsize=leg_charsize, box=0
  
;  mess2 = !csym.chi+'!U2!N/'+!csym.nu+' = '+ntostr(minchisq,4,/round)+'/'+ntostr(degfree)+$
;    ' = '+ntostr(minchisq/degfree,4,/round)
;  legend, mess2,charsize=leg_charsize, box=0
  
  print,'Best fit Power Index: '+ntostr(bestpow)+'+'+ntostr(phighdiff)+'-'+ntostr(plowdiff)
  print,'Best fit Power Norm: '+ntostr(bestnorm)+'+'+ntostr(nhighdiff)+'-'+ntostr(nlowdiff)
  
END 

PRO plot_zshear_old_jackknife, lensum, rlensum, combstruct, rcombstruct

  ;; want at least 
  Nsub = long( 500.*n_elements(lensum)/127000. ) > 50

  ;;Nsub = 500L

  print
  print,'Jackknifing   Nsub: '+ntostr(Nsub)
  
  print,'Doing lenses'
  shearjackknife_bystripe, lensum, rlensum, $
                           sigma, sigmaerr, covariance, $
                           rsigma, rsigmaerr, rcovariance, $
                           Nsub=Nsub

;wuse = (lindgen(18))[6:17]
;shearjackknife_bystripe,lensum,rlensum,ts,tse,tcov,nsub=300,wuse=wuse
;stop
  print,'Doing Ortho'
  shearjackknife_bystripe, lensum, rlensum, $
                           orthosig, orthosigerr, orthocov, $
                           ttrsigma, ttrsigmaerr, ttrcovariance, $
                           Nsub=Nsub, /ortho

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Subtract the random points: must include
  ;; their errors
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  sigma = sigma - rsigma
  sigmaerr = sqrt(sigmaerr^2 + rsigmaerr^2)

  ;;combstruct.sigma = sigma
  combstruct.sigmaerr = sigmaerr
  combstruct.orthosigerr = orthosigerr
  rcombstruct.sigmaerr = rsigmaerr

  calc_correlation_matrix, covariance, corr, status
  IF status NE 0 THEN message,'Doh!'
  calc_correlation_matrix, rcovariance, rcorr, rstatus
  IF rstatus NE 0 THEN message,'Doh!'

  ;; increase errors on covariance matrix accordingly
  calc_covfromcorr_matrix, corr, sigmaerr, covariance
  
  add_tag, temporary(combstruct), 'covariance', covariance, $
           combstruct
  add_tag, temporary(combstruct), 'corr', corr, $
           combstruct
  add_tag, temporary(combstruct), 'orthocov', orthocov, $
           combstruct
  add_tag, temporary(rcombstruct), 'covariance', rcovariance, $
           rcombstruct
  add_tag, temporary(rcombstruct), 'corr', rcorr, $
           rcombstruct

END 

PRO plot_zshear_old_rebin, comb, rcomb, combrebin, rcombrebin

  ;; combined files
  zshear_rebin, comb, meanr, sigma, sigmaerr, wsum, orthosig, orthosigerr
  zshear_rebin, rcomb, rmeanr, rsigma, rsigmaerr, rwsum, rorthosig, rorthosigerr

  combrebin = create_struct('meanr',       meanr, $
                            'sigma',       sigma, $
                            'sigmaerr',    sigmaerr, $
                            'wsum',        wsum, $
                            'orthosig',    orthosig, $
                            'orthosigerr', orthosigerr)

  rcombrebin = create_struct('meanr',       rmeanr, $
                             'sigma',       rsigma, $
                             'sigmaerr',    rsigmaerr, $
                             'wsum',        rwsum, $
                             'orthosig',    rorthosig, $
                             'orthosigerr', rorthosigerr)

  add_tag, temporary(comb), 'meanr_rebin', meanr, comb
  add_tag, temporary(comb), 'sigma_rebin', sigma, comb
  add_tag, temporary(comb), 'sigmaerr_rebin', sigmaerr, comb
  add_tag, temporary(comb), 'wsum_rebin', wsum, comb
  add_tag, temporary(comb), 'orthosig_rebin', orthosig, comb
  add_tag, temporary(comb), 'orthosigerr_rebin', orthosigerr, comb

  add_tag, temporary(rcomb), 'meanr_rebin', rmeanr, rcomb
  add_tag, temporary(rcomb), 'sigma_rebin', rsigma, rcomb
  add_tag, temporary(rcomb), 'sigmaerr_rebin', rsigmaerr, rcomb
  add_tag, temporary(rcomb), 'wsum_rebin', rwsum, rcomb
  add_tag, temporary(rcomb), 'orthosig_rebin', rorthosig, rcomb
  add_tag, temporary(rcomb), 'orthosigerr_rebin', rorthosigerr, rcomb

END 


PRO plot_zshear_old, files, rfiles, paramname, combname, rcombname, $
                 lensumfiles=lensumfiles, rlensumfiles=rlensumfiles, $
                 zfiles=zfiles, $
                 struct=struct, rstruct=rstruct, $
                 combstruct=combstruct, $
                 ptitle=ptitle,$
                 sigvrange=sigvrange, cutrange=cutrange, $
                 normrange=normrange, powrange=powrange, $
                 munit=munit, dojack=dojack, lumclr_in=lumclr_in, $
                 comoving=comoving, logbin=logbin

  IF n_params() LT 5 THEN BEGIN 
      print,'-Syntax: plot_zshear_old, files, rfiles, paramname, combname, rcombname, $'
      print,'        lensumfiles=lensumfiles, rlensumfiles=rlensumfiles, $'
      print,'        zfiles=zfiles, $'
      print,'        struct=struct, rstruct=rstruct, $'
      print,'        combstruct=combstruct, $'
      print,'        ptitle=ptitle,$'
      print,'        sigvrange=sigvrange, cutrange=cutrange, $'
      print,'        normrange=normrange, powrange=powrange, $'
      print,'        munit=munit, dojack=dojack, lumclr_in=lumclr_in, $'
      print,'        comoving=comoving'
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;
  ;; Some parameters
  ;;;;;;;;;;;;;;;;;;;
  
  IF n_elements(lumclr_in) NE 0 THEN lumclr=lumclr_in

  IF keyword_set(logbin) THEN BEGIN 
      xlog=1
      ylog=1
      xtickf='loglabels'
      ytickf='loglabels'
  ENDIF 

  setup_mystuff

  leg_charsize=1.0

  !p.noerase=0

  ;; Lets plot some things in color
  simpctable

  xmarg = !x.margin
  !x.margin = [10,10]

  fac = 1000.                   ;convert to Mpc to avoid underflow
  ;; fit guesses for Clusters
  siguess = [40., -.7]
  shguess = [.005, -.7]

  ;; these are cluster defaults
  IF n_elements(sigvrange) EQ 0 THEN sigvrange = [100., 800.]     ;km/s
  IF n_elements(cutrange) EQ 0 THEN cutrange  = [100., 3000.]
  IF n_elements(normrange) EQ 0 THEN normrange = [0., 100.]        ;Msolar/pc^2
  IF n_elements(powrange) EQ 0 THEN powrange = [-.2, 1.5]
  IF n_elements(nresample) EQ 0 THEN nresample = 10000L
  IF n_elements(munit) EQ 0 THEN munit=1.e14

  psym=8
  nclr = 3

  xtitle=!kpcxtitle

  nstring_a = 7
  nstring_b = 5

  ;;;;;;;;;;;;;;;;;;;;;
  ;; get redshifts
  ;;;;;;;;;;;;;;;;;;;;;
  
  nz = n_elements(zfiles)

  IF nz NE 0 THEN BEGIN 

      FOR i=0, nz-1 DO BEGIN 
              
          sz = mrdfits(zfiles[i], 1, /silent)
          
          IF i EQ 0 THEN BEGIN 
              zstruct = sz
          ENDIF ELSE BEGIN
              concat_structs, temporary(zstruct), sz, zstruct
          ENDELSE 
          
      ENDFOR 
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; correct data, read in the corrected files
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  correct_shear_files, files, rfiles, corrfiles

  ;;;;;;;;;;;;;;;;;;;;;
  ;; Get combined data
  ;;;;;;;;;;;;;;;;;;;;;

  combine_zshear, corrfiles, combstruct, comoving=comoving, $
                  already_comoving=already_comoving
  combine_zshear, rfiles, rcombstruct, comoving=comoving, $
                  already_comoving=already_comoving

  IF n_elements(ptitle) EQ 0 THEN ptitle=''
  message1 = 'Nlenses: '+ntostr(long(combstruct.nlenses))


  IF (n_elements(lensumfiles) NE 0) AND (n_elements(rlensumfiles) NE 0) THEN BEGIN 

      IF keyword_set(dojack) THEN BEGIN 
          columns=['zindex', 'ceta','clambda','wsum','sigma']
          combine_zlensum_files, lensumfiles, lensum
          combine_zlensum_files, rlensumfiles, rlensum, columns=columns

          plot_zshear_old_jackknife, lensum, rlensum, combstruct, rcombstruct

          ;; should we convert to comoving?
          IF (keyword_set(comoving) AND $
              NOT keyword_set(already_comoving) ) THEN BEGIN 

              zmean = combstruct.zmean
              combstruct.sigmaerr = combstruct.sigmaerr/(1+zmean)^2
              combstruct.covariance = combstruct.covariance/(1+zmean)^4

              combstruct.orthosigerr = combstruct.orthosigerr/(1+zmean)^2
              combstruct.orthocov = combstruct.orthocov/(1+zmean)^4

              rzmean = rcombstruct.zmean
              rcombstruct.sigmaerr = rcombstruct.sigmaerr/(1+rzmean)^2

          ENDIF 

          ;; free this memory, it may be huge
          delvarx,rlensum

      ENDIF ELSE message,'You sent lensumfiles but did not set /dojack'
      
  ENDIF

  print,'Rebinning lensums and combined structs'
  plot_zshear_old_rebin, combstruct, rcombstruct, combrebin, rcombrebin

                                ;Find out yrange
  yrange=prange(combstruct.sigma, combstruct.orthosig,$
                combstruct.sigmaerr,combstruct.orthosigerr,/slack)
  
  ;; redshift distribution
  IF n_elements(zstruct) NE 0 THEN BEGIN
      erase 
      zrange = [0.0, 0.4]
      !p.multi=[0,0,2]
      plothist, zstruct.z, bin=.01, line=0, $
        xtitle='redshift',ytitle='N(z)', title=ptitle;, xrange=zrange
      !p.multi=0
      zmess = 'Nlenses: '+ntostr(long(combstruct.nlenses))
      legend, zmess, charsize=leg_charsize,/right
      IF !d.name eq 'X' THEN key=get_kbrd(1)
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;
  ;; Shear plots
  ;;;;;;;;;;;;;;;;;;;
  
  nbin = n_elements(combstruct.meanr)

  erase 


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; was wuse sent?
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(wuse) EQ 0 THEN BEGIN 
      wuse = lindgen(nbin)
  ENDIF ELSE BEGIN 
      nbin = n_elements(wuse)
  ENDELSE 
  
  IF n_elements(wuserebin) EQ 0 THEN BEGIN 
      wuserebin = lindgen(n_elements(combrebin.meanr))
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; power law fits
  ;;;;;;;;;;;;;;;;;;;;;;;;;
      
  xt='Cutoff Radius (h!U'+!csym.minus+'1!N kpc)'
  sigvsym=!csym.sigma+'!DV!N'
  yt=sigvsym+' (km/s)'
  IF nbin GT 1 THEN BEGIN 

      plot_zshear_old_fitpower, combstruct, bestpow, bestnorm, $
                            powlow, powhigh, normlow, normhigh, $
                            sigallow_low, sigallow_high, combsigdiff_fit, $
                            dojack=dojack, wuse=wuse

      IF !d.name eq 'X' THEN key=get_kbrd(1)

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; <Sigma>-sigma plots: lower panel
      ;; is orthosig
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;      w = where(combstruct.meanr GE 200,nw)
;      IF nw NE 0 THEN BEGIN 
;          obotyrange = $
;            prange(combstruct.orthosig[w],combstruct.orthosigerr[w],$
;                   /symm,/slack)
;      ENDIF 

;      omax = max(abs(combstruct.orthosig))
      to = prange(combstruct.orthosig,combstruct.orthosigerr,/symm,/slack)
      omax = 2.0*sdev(combstruct.orthosig) < abs(to[0])
      obotyrange = [-omax,omax]
      plotboth_density_contrast, combstruct, logbin=logbin, $
                                 xfit=combstruct.meanr[wuse], $
                                 yfit=combsigdiff_fit,$
                                 botyrange=obotyrange

;      IF !d.name eq 'X' THEN key=get_kbrd(1)

;      plot_zshear_old_fitpower, combrebin, bestp, bestn, $
;                            plow, phigh, nlow, nhigh, $
;                            siglow, sighigh, $
;                            wuse=wuserebin
      
      IF !d.name eq 'X' THEN key=get_kbrd(1)

      plotboth_density_contrast, combrebin, logbin=logbin, $
                                 xfit=combstruct.meanr, yfit=combsigdiff_fit,$
                                 botyrange=obotyrange

      IF !d.name eq 'X' THEN key=get_kbrd(1)

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; truncated isothermal fits
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      bestcut=0.0
      cutlow=0.0
      cuthigh=0.0

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; isothermal fits
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      sismodel, sigvrange, modelx, model, sigv
      ;;print,'------****** ',min(sigv),max(sigv)
      print
      print,'Pure SIS fits'
      chisq_conf_1par, combstruct.meanr,combstruct.sigma,combstruct.sigmaerr,$
                       modelx, model, sigv, $
                       chisq_surf,bestsig,siglow,sighigh,aspect=1,/center, $
                       xtitle=yt, $
                       ytitle=!csym.delta_cap+!csym.chi+'!U2!N/'+!csym.nu, $
                       charsize=1.5
      
      clevel = 0
        
      highdiff = sighigh[clevel]-bestsig & lowdiff = bestsig - siglow[clevel]
      mess = sigvsym+' = '+ntostr(bestsig,6,/round)+$
        '!S!U'+!csym.plus+ntostr(highdiff,5,/round)+'!R!D'+!csym.minus+ntostr(lowdiff,5,/round)
      legend, mess, charsize=leg_charsize
      
      
      IF !d.name eq 'X' THEN key=get_kbrd(1)
      print
      print,'plotting data and best fitting pure sis'
      
      plotmod = sigmasis( bestsig, combstruct.meanr )
      plotboth_density_contrast, combstruct, logbin=logbin, $
                                 xfit=combstruct.meanr,$
                                 yfit=plotmod
      
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; fit ortho to zero
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
      print
      print,'Fitting ortho to zero'
      IF keyword_set(dojack) THEN BEGIN 
          diff = combstruct.orthosig
          cinv = invert(combstruct.orthocov)
          chisq_ortho = (transpose(diff)#(reform(cinv##diff)))[0]
          chisq_ortho_per = chisq_ortho/nbin
      ENDIF ELSE BEGIN 
          chisq_ortho = total( (combstruct.orthosig/combstruct.orthosigerr)^2 )
          chisq_ortho_per = chisq_ortho/nbin
      ENDELSE 
      print,'Fit to zero: '+ntostr(chisq_ortho)+'/'+ntostr(nbin)+$
            ' = '+ntostr(chisq_ortho_per)
      
  ENDIF 

  IF !d.name eq 'X' THEN key=get_kbrd(1)
  
  tsigma_sis_fit, combstruct, sismass, sismasserr, sissigma2,sissigmaerr2

  ;; mean luminosity

  ;; we want to do luminosity thing, even if lumclr not sent, as long
  ;; as lensum is there

  IF n_elements(lumclr) EQ 0 THEN lumclr = 2

  IF n_elements(lensum) NE 0 THEN BEGIN
      IF tag_exist(lensum[0],'lum') THEN BEGIN 

          plot_zshear_old_meanlum_gmr, lensum, combstruct, lumclr, $
                                   tmeanlum, tmeanlumerr, meanlum, meanlumerr, $
                                   tmeangmr, tmeangmrerr, meangmr, meangmrerr, $
                                   lumfrac, lumfracerr, gmrfrac, gmrfracerr
          
          IF !d.name eq 'X' THEN key=get_kbrd(1)
      ENDIF 
  ENDIF 
  
  !p.multi = [0,0,2]

  xr = [0.8*min(combstruct.meanr), 1.2*max(combstruct.rmax_act)]
  plot, combstruct.meanr, combstruct.nlbin, xtitle=xtitle, $
        ytitle='Nlenses', xlog=xlog, xtickf=xtickf, xrange=xr, xsty=1+2,xticklen=0.04
  oplot, combstruct.meanr, combstruct.nlbin, psym=psym

  maxl = max(combstruct.nlbin)
  plot, combstruct.meanr, combstruct.nlbin/maxl, xtitle=xtitle, $
        ytitle='Fraction of Lenses', xlog=xlog, xtickf=xtickf, xrange=xr, xsty=1+2,xticklen=0.04,yrange=[0,1.2]
  oplot, combstruct.meanr, combstruct.nlbin/maxl, psym=psym

  !p.multi=0
  IF !d.name eq 'X' THEN key=get_kbrd(1)

  IF nbin GT 1 THEN BEGIN 

      ;; Plot mass
      IF (munit EQ 1.e12) THEN mstr = '12' ELSE mstr = '14'

      slack = 0.5
      minfac=1.-slack
      maxfac=1.+slack
      xrange=[minfac*min(combstruct.meanr), maxfac*max(combstruct.meanr)]
      yrange = prange(sismass/munit, sismasserr/munit,slack=slack)
      yrange[0] = yrange[0] > 0.01
      ytit='M(R)  (10!U'+mstr+'!N h!U'+!csym.minus+'1!N M'+sunsymbol()+')'
      aploterror, !gratio, $
        combstruct.rmax_act, $
        sismass/munit, sismasserr/munit, psym=psym, $
        ytitle=ytit, xtitle=xtitle, /center,$
        xlog=xlog, ylog=ylog,$
        yst=1+2,xst=1+2, yrange=yrange, xrange=xrange,$
        xtickf=xtickf,ytickf=xtickf,xticklen=0.04,yticklen=0.04
      oplot,[0,10000],[0,0]
      
      ;; overplot best-fit SIS
;      yfit = 7.31e14*(bestsig/1000.)^2*(combstruct.rmax_act/1000.)/munit
;      oplot, combstruct.rmax_act, yfit
      
  ENDIF 

  IF !d.name eq 'X' THEN key=get_kbrd(1)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Output fits: move below
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  print,'--------------------------------'
  print,'Output file: ',paramname
  print,'--------------------------------'

  IF nbin EQ 1 THEN BEGIN 
      paramstruct = create_struct('sismass', sismass, $
                                  'sismasserr', sismasserr)
  ENDIF ELSE BEGIN 
      paramstruct = create_struct('norm', bestnorm, $
                                  'normlow', normlow, $
                                  'normhigh', normhigh, $
                                  'power', bestpow, $
                                  'powlow', powlow, $
                                  'powhigh', powhigh, $
                                  'yfit', combsigdiff_fit, $
                                  'yallow_low', sigallow_low,$
                                  'yallow_high', sigallow_high, $
                                  'sigmav', bestsig, $
                                  'sighigh', sighigh, $
                                  'siglow', siglow, $
                                  'cutoff', bestcut, $
                                  'cuthigh', cuthigh, $
                                  'cutlow', cutlow, $
                                  'sismass', sismass, $
                                  'sismasserr', sismasserr, $
                                  'sissigma2',sissigma2,$
                                  'sissigmaerr2',sissigmaerr2)

      IF n_elements(meanlum) NE 0 THEN BEGIN 
          paramstruct = create_struct(paramstruct, $
                                      'meanlum', meanlum, $
                                      'meanlumerr', meanlumerr, $
                                      'tmeanlum', tmeanlum, $
                                      'tmeanlumerr', tmeanlumerr, $
                                      'lumfrac', lumfrac, $
                                      'lumfracerr', lumfracerr, $
                                      'meangmr', meangmr, $
                                      'meangmrerr', meangmrerr, $
                                      'tmeangmr', tmeangmr, $
                                      'tmeangmrerr', tmeangmrerr, $
                                      'gmrfrac', gmrfrac, $
                                      'gmrfracerr', gmrfracerr)
      ENDIF 
  ENDELSE 

  combine_structs, temporary(combstruct), paramstruct, combstruct

  IF nbin GT 1 THEN BEGIN 
      paramstruct = create_struct(paramstruct, $
                                  'rmax_act', combstruct.rmax_act, $
                                  'tsigma', combstruct.tsigma, $
                                  'tsigmaerr', combstruct.tsigmaerr)
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; calculate wgm, add wgm and r0 to combstruct
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  plot_zshear_old_getwgm, combstruct
;  minfac=0.5
;  maxfac=1.5
  
;  wgmmeanr = combstruct.meanr/1000.
;  wgmxr=[minfac*min(wgmmeanr), maxfac*max(wgmmeanr)]
;  IF wgmxr[0] LT 0.1 THEN BEGIN 
;      IF wgmxr[0] GT 0.01 THEN wgmxr[0] = 0.01
;  ENDIF 
;  wgmyr=prange(combstruct.wgm, combstruct.wgmerr,slack=0.5) > 1

;  aploterror, 1.0, wgmmeanr, combstruct.wgm, combstruct.wgmerr, $
;              psym=8, xrange=wgmxr, /xsty, yrange=wgmyr, /ysty, $
;              xtitle=!mpcxtitle, ytitle=!wgmytitle, /xlog, xtickf='loglabels', /ylog,ytickf='loglabels', $
;              xticklen=0.04, yticklen=0.04
;  wgmfac = combstruct.wgm[0]/combstruct.sigma[0]
;  oplot, wgmmeanr, combstruct.yfit*wgmfac
;  legend, !csym.omega_cap+'!Dm!N = '+ntostr(!omegam,4,/round)+!csym.plusminus+ntostr(!omegamerr,4,/round), /bottom, box=0,charsize=leg_charsize
  
;  IF !d.name eq 'X' THEN key=get_kbrd(1)


  print,'Output param file: ',paramname
  mwrfits, paramstruct, paramname, /create
  print,'Output combined file: ',combname
  mwrfits, combstruct, combname, /create
  print,'Output rand comb file: ',rcombname
  mwrfits, rcombstruct, rcombname, /create

  ;;;;;;;;;;;;;;;;;;;;;;;
  ;; Random point plots
  ;;;;;;;;;;;;;;;;;;;;;;;
      
  print
  print,'Fitting random to zero'
  ;; first fit individual bandpasses.
  
  IF keyword_set(dojack) THEN BEGIN 
      diff = rcombstruct.sigma
      cinv = invert(rcombstruct.covariance)
      chisq_rand = (transpose(diff)#(reform(cinv##diff)))[0]
      chisq_rand_per = chisq_rand/nbin
  ENDIF ELSE BEGIN 
      chisq_rand = total( (rcombstruct.sigma/rcombstruct.sigmaerr)^2 )
      chisq_rand_per = chisq_rand/nbin
  ENDELSE 
  print,'Fit to zero: '+ntostr(chisq_rand)+'/'+ntostr(nbin)+$
        ' = '+ntostr(chisq_rand_per)

  wtt=where(rcombstruct.meanr GE 40)
  tr = prange(rcombstruct.sigma[wtt],rcombstruct.sigmaerr[wtt],/symm,/slack)
  rndmax = 2.0*sdev(rcombstruct.sigma[wtt]) < abs(tr[0])
  rbotyrange = [-rndmax,rndmax]

  aploterror, !gratio, $
              rcombstruct.meanr, rcombstruct.sigma, rcombstruct.sigmaerr, $
              psym=psym, xtitle=xtitle, ytitle=!randdeltaytitle, $
              title=ptitle, ystyle=1, xlog=xlog, xtickf=xtickf, /center,$
              yrange=rbotyrange,xticklen=0.04
  oplot,[1,10000],[0,0]

  rmess = !csym.chi+'!U2!N/'+!csym.nu+' = '+ntostr(chisq_rand,4,/round)+'/'+ntostr(nbin)+$
    ' = '+ntostr(chisq_rand_per,4,/round)
  legend, rmess, /right, charsize=leg_charsize,box=0


  IF !d.name eq 'X' THEN key=get_kbrd(1)

  ;; xyouts here only calibrated for ps

;  tr = prange(rcombstruct.sigma,rcombstruct.sigmaerr,/symm,/slack)

;  rndmax = 2.0*sdev(rcombstruct.sigma) < abs(tr[0])
;  rbotyrange = [-rndmax,rndmax]
  plotboth_density_contrast, combstruct, rcombstruct, logbin=logbin, $
                             xfit=combstruct.meanr, yfit=combsigdiff_fit,$
                             botyrange=rbotyrange

  plotboth_density_contrast, combrebin, rcombrebin, logbin=logbin, $
                             xfit=combstruct.meanr, yfit=combsigdiff_fit,$
                             botyrange=rbotyrange

  plotboth_density_contrast, combstruct, rcombstruct, logbin=logbin, $
                             xoplot=rcombstruct.meanr,$
                             yoplot=rcombstruct.sigma,$
                             oploterr=rcombstruct.sigmaerr,$
                             oplotsym=4,$
                             oplotcolor=!blue,$
                             botyrange=rbotyrange

  IF !d.name eq 'X' THEN key=get_kbrd(1)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Fractional Excess plots
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  minfac=0.7
  maxfac=1.3
  xrange=[minfac*min(combstruct.meanr), maxfac*max(combstruct.meanr)]/1000.


  mincorr = 1.e-4
  yrange=prange(combstruct.clust_corr-1.0,$
                combstruct.clust_corr_err,/slack)
  IF keyword_set(logbin) THEN yrange = yrange > mincorr

  ticklen=0.03
  aploterror, 1, combstruct.meanr/1000., combstruct.clust_corr-1, $
    combstruct.clust_corr_err, $
    yrange=yrange, xrange=xrange,psym=8,$
    ystyle=1+2,xstyle=1+2,$
    ytitle='C(R) '+!csym.minus+' 1', xtitle=!mpcxtitle,$
    xticklen=ticklen, yticklen=ticklen, xlog=xlog, xtickf=xtickf, $
    ylog=ylog,ytickf=ytickf, /center

  print

  IF keyword_set(dojack) THEN BEGIN 
      loadct,0
      cxtit = 'Radial Bin'
      cytit = cxtit
      IF !d.name EQ 'X' THEN key=get_kbrd(1)
      tvim2, combstruct.corr, max_color=!grey100, /scale, $
             xtit=cxtit,ytit=cytit
      IF !d.name EQ 'X' THEN key=get_kbrd(1)
      tvim2, rcombstruct.corr, max_color=!grey100, /scale, $
             xtit=cxtit,ytit=cytit
  ENDIF 

  !x.margin = xmarg 
  
return
END 
