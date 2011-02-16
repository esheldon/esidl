PRO fit_wtheta_tot_lumw,wclr,wtlen,wtrlen, sigonly=sigonly, noweight=noweight, makeplot=makeplot

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: fit_wtheta_tot_lumw,wclr,wtlen,wtrlen, makeplot=makeplot'
      return
  ENDIF 

  setup_mystuff

  indir='/sdss5/data0/lensout/stripe10/'
  outdir = '/sdss5/data0/lensout/wtheta_conv/'

  CASE 1 OF
      keyword_set(sigonly): wstr='sigonly_'
      keyword_set(noweight): wstr='noweight_'
      ELSE: wstr=''
  ENDCASE 

  fend='N2.fit'

  ;; output .fit file
  outfile = outdir+'fit_wtheta_tot_'+wstr+!colors[wclr]+'w_lumw.fit'
  outpsfile = outdir+'fit_wtheta_tot_'+wstr+!colors[wclr]+'w_lumw.ps'
  print
  print,'Output fits file: ',outfile
  print,'Output psfile: ',outpsfile
  print

  IF keyword_set(makeplot) THEN begplot,name=outpsfile

  minpow = -1.1 & maxpow = -0.3
  minnorm = 0.2 & maxnorm = 5.0

  IF wclr EQ 0 THEN BEGIN 
      minpow=-1.1 & maxpow = .3
      minnorm = 0.0 & maxnorm = 5.0
  ENDIF 

  npow = 400L
  nnorm = 400L
  powvals = arrscl( findgen(npow), minpow, maxpow )
  normvals = arrscl( findgen(nnorm), minnorm, maxnorm )

  delvarx, combbeta, combnorm, combpower, combnormlow, combnormhigh, combpowlow, combpowhigh

  setupplot, dtype

  stripestr = ['stripe10']
  nstripe = n_elements(stripestr) - 1
  type='wtheta'

  tfile=indir+'wthetalumw_'+wstr+'stripe10_'+!colors[wclr]+'w_'+fend
  trfile=indir+'wthetarandlumw_'+wstr+'stripe10_'+!colors[wclr]+'w_'+fend
  wtcomb=mrdfits(tfile,1)
  wtrcomb=mrdfits(trfile,1)
  wtsum=mrdfits(indir+'wthetalumw_'+wstr+'stripe10_sum_'+!colors[wclr]+'w_'+fend,1)
  wtrsum=mrdfits(indir+'wthetarandlumw_'+wstr+'stripe10_sum_'+!colors[wclr]+'w_'+fend,1)
  IF n_elements(wtlen) EQ 0 THEN BEGIN 
      wtlen=mrdfits(indir+'wthetalumw_'+wstr+'stripe10_lensum_'+!colors[wclr]+'w_'+fend,1)
      wtrlen=mrdfits(indir+'wthetarandlumw_'+wstr+'stripe10_lensum_'+!colors[wclr]+'w_'+fend,1)
  ENDIF 
  nrsum=n_elements(wtlen[0].rsum)
  wterr=fltarr(nrsum)
  wtrerr=wterr

  w=where(wtcomb.meanr NE 0.,nbin)
  
  nxx = 1000
  xx = arrscl( findgen(nxx), min(wtcomb.meanr[w])/1000., $
               max(wtcomb.meanr[w])/1000. )

  numbeta=n_elements(wtcomb.beta)
  FOR bn=0L,numbeta-1  DO BEGIN 

      diffdensecomb = (wtcomb.density[w,bn]-wtrcomb.density[w,bn])*1000.^2
      
      print
      print,'beta = ',wtcomb.beta[bn]
      title = !tsym.beta+' = '+ntostr(wtcomb.beta[bn],5)

      FOR i=0L, nrsum-1 DO BEGIN 
          
          wterr[i] = sqrt( total( (wtlen.npsum[i,bn]-$
                                   wtlen.wsum[i]*wtcomb.npair[i,bn])^2 ))/wtsum.wsum[i]
          
          wtrerr[i] =sqrt( total( (wtrlen.npsum[i,bn]-$
                                   wtrlen.wsum[i]*wtrcomb.npair[i,bn])^2 ))/wtrsum.wsum[i]
          
      ENDFOR 
      errorcomb = sqrt(wterr^2 + wtrerr^2)/wtcomb.area*1000.^2
;      aploterror,!gratio,wtcomb.meanr[w],diffdensecomb,errorcomb[w],psym=4,$
;        ytitle='Overdensity (# Mpc!U-2!N)',xtitle='Radius (kpc)',$
;        title=title
      
      cumul = fltarr(nbin)
      FOR i=0L, nbin-1 DO BEGIN 
          IF i EQ 0 THEN BEGIN
              cumul[i] = (wtcomb.npair[i,bn] - wtrcomb.npair[i,bn])
          ENDIF ELSE BEGIN 
              cumul[i] = cumul[i-1] + (wtcomb.npair[i,bn] - wtrcomb.npair[i,bn])
          ENDELSE 
      ENDFOR 
      
;      oplot,wtcomb.meanr[w], cumul
;      legend,['Overdensity','Cumulative'],psym=[4,0],/center,/top,charsize=1.0,$
;        thick=[!p.thick,!p.thick]
;      if dtype eq 'X' then key=get_kbrd(1)


      pow_chisq_conf, wtcomb.meanr[w[1:nbin-1]]/1000., $
        diffdensecomb[w[1:nbin-1]], $
        errorcomb[w[1:nbin-1]], $
        powvals, normvals, $
        chisq_surfcomb, $
        pmincomb, nmincomb, $
        powlowcomb, powhighcomb, normlowcomb, normhighcomb, $
        powallow=powallowcomb,normallow=normallowcomb,$
        title=title
      print,'norm comb = ',nmincomb
      print,'pow comb = ',pmincomb
      
      yycomb = nmincomb*(xx)^pmincomb
      if dtype eq 'X' then key=get_kbrd(1)
      
      aploterror,!gratio,wtcomb.meanr[w],diffdensecomb[w],errorcomb[w],psym=4,$
        ytitle='Overdensity (# Mpc!U'+!tsym.minus+'2!N)',$
        xtitle='Projected Radius (h!U'+!tsym.minus+'1!N kpc)',$
        title=title
      oplot,xx*1000.,yycomb
      
      oplot,wtcomb.rmax_act[w],cumul,line=2
      legend,['Overdensity','Cumulative'],line=[-1,2],/center,/top,charsize=1.0,$
        thick=[!p.thick,!p.thick],box=0
      legend,['Overdensity','Cumulative'],psym=[4,3],/center,/top,charsize=1.0,$
        thick=[!p.thick,!p.thick],box=0

      confind=0
      add_arrval, wtcomb.beta[bn], combbeta
      add_arrval, nmincomb, combnorm
      add_arrval, pmincomb, combpower
      add_arrval, normlowcomb[confind], combnormlow
      add_arrval, normhighcomb[confind], combnormhigh
      add_arrval, powlowcomb[confind], combpowlow
      add_arrval, powhighcomb[confind], combpowhigh

      if dtype eq 'X' then BEGIN 
          key=get_kbrd(1)
          IF key EQ 'q' THEN return
      ENDIF 
  ENDFOR 

  if dtype eq 'X' then key=get_kbrd(1)

  erase & multiplot,[1,2]
  plot,combbeta,combnorm,yrange=[0,4],ytitle='Norm'
  oplot,combbeta,combnormlow,line=2
  oplot,combbeta,combnormhigh,line=2

  multiplot

  plot,combbeta, combpower,yrange=[-1.0,-0.5],ytitle='Power',xtitle=!tsym.beta
  oplot,combbeta,combpowlow,line=2
  oplot,combbeta,combpowhigh,line=2

  multiplot,/reset

  IF keyword_set(makeplot) THEN endplot

  ;; make fits file
  print
  print,'Output fits file: ',outfile
  print
  
  st=create_struct('beta',combbeta,$
                   'norm',combnorm,$
                   'power',combpower,$
                   'normlow',combnormlow,$
                   'normhigh',combnormhigh,$
                   'powlow',combpowlow,$
                   'powhigh',combpowhigh)

  mwrfits, st, outfile, /create


END 
