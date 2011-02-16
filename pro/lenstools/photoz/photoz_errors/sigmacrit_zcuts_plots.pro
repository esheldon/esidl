PRO sigmacrit_zcuts_boot, mdenscont, nlens, Nresamp, weights, datamean, dataerr

  Nvar=1L
  resample_data, mdenscont, nlens, Nvar, Nresamp, resamp, $
                 weights=weights

  bootstrap, resamp, Nresamp, Nvar, datamean, dataerr

END 

PRO sigmacrit_zcuts_plots, file, doplot=doplot, refix=refix, unfix=unfix, nosum=nosum

  ;; use /nosum if mdenscontsum was not calculated
  ;; e.g. 
  ;; dir = '/net/cheops1/data0/esheldon/test_sigmacriterr/sim/'
  ;; file = dir + 'test_sigmacrit_N3.fit'

  t=mrdfits(file, 1)
  IF keyword_set(unfix) THEN t.mdenscont = t.mdenscont*t.fracgood
  IF keyword_set(refix) THEN BEGIN 
      rfstr = 'refix'
      t.mdenscont = t.mdenscont*t.fracgood/t.fracgoodtrue
  ENDIF ELSE rfstr = ''
  IF n_elements(t.indices) NE n_elements(t.mdenscont) THEN message,'Adjust for indices'

  maxhist = 100

  fitfile = repstr(file, 'test', 'fits'+rfstr+'_test')
  IF keyword_set(doplot) THEN BEGIN 
      psfile = repstr(fitfile, '.fit', '.ps')
      begplot, name=psfile, /color
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; source redshifts used
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; for case where we used all sources, this is
  ;; the redshift distribution
  deltazs = t.zs[1]-t.zs[0]

  ;; we made a redshift cut at nsig*sigzs
  deltazsi = t.zsi[1]-t.zsi[0]

  ;; plot them together
  ytitle = 'P(z!Ds!N)'
  xtitle = 'z!Ds!N'
  aplot, !gratio, t.zsi, t.pofzitot, $
        xtitle=xtitle, ytitle=ytitle
  oplot, t.zsi, t.pofzitottrue, color=!blue
  oplot, t.zs, t.pofzs, color=!red
  legend,['All', $
          'Zcut est: '+ntostr(t.nsig,3)+'*'+!csym.sigma+'(z!Ds!N)', $
          'Zcut true'],$
         /right, line=[0,0,0], color=[!red,!p.color,!blue],$
         thick=[!p.thick,!p.thick,!p.thick],$
         charsize=1

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; When we make zcuts, we are biased towards objects that have
  ;; scattered toward the high side of their distribution. Plot the
  ;; estimated versus true for the case where we cut on redshift
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF !d.name EQ 'X' THEN key=get_kbrd(1)
  
  title = 'Cut based on estimated redshift > '+ntostr(t.nsig,3)+'*'+$
    !csym.sigma+'(z!Ds!N)'
  xtitle = 'z!Ds!N True'
  ytitle = 'z!Ds!N estimated'
;  plot, t.zstrue,t.zsest, psym=3, $
;        xtitle=xtitle,ytitle=ytitle,title=title,/iso
;  oplot, [0, max(t.zstrue)],[0, max(t.zstrue)], color=!red

  ploth, t.zstrue, t.zsest, xtitle=xtitle,ytitle=ytitle,title=title, $
         xrange = [0, max(t.zstrue)], yrange=[0, max(t.zstrue)],/silent
  oplot, [0, max(t.zstrue)],[0, max(t.zstrue)], color=!red

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Plot estimated density contrast for each lens, using the two 
  ;; methods. Spread estimated multiple ways.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF !d.name EQ 'X' THEN key=get_kbrd(1) ELSE !p.multi=[0,0,2]

  bin=0.05

  ;; for bootstrap
  Nresamp = 1000L

  ;; redshift cut
  w=where(t.nsources NE 0)
  ;;sigma = sqrt(t.nsources[w])
  sigma = 1./sqrt(t.nsources[w])
  weights = 1./sigma^2
  ;; sdev is error here
  wmom, t.mdenscont[w], sigma, mean_mdenscont, sig_mdenscont, err_mdenscont, $
        /calcerr

  ;; all sources
  mdenscont2 = t.mdenscontsum/t.nsourcesum
  w2=where(t.nsourcesum NE 0, nw2)

  psig = 4.
  IF nw2 NE 0 THEN BEGIN 
      sigma2 = 1./sqrt(t.nsourcesum[w2])
      weights2 = 1./sigma2^2
      wmom, mdenscont2[w2], sigma2, mean_mdenscont2, sig_mdenscont2, $
            err_mdenscont2, $
            /calcerr

      minxhist = mean_mdenscont2-psig*sig_mdenscont2
      maxxhist = mean_mdenscont2+psig*sig_mdenscont2
  ENDIF ELSE BEGIN 
      minxhist = mean_mdenscont-psig*sig_mdenscont
      maxxhist = mean_mdenscont+psig*sig_mdenscont
  ENDELSE 

 
  xtitle = !csym.delta_cap+!csym.sigma_cap
  ytitle = 'N'
  plothist, t.mdenscont, xhist, yhist, bin=bin, /norm, charsize=1,$
            xtitle=xtitle,ytitle=ytitle, $
            min=minxhist, $
            max=maxxhist

  ;; oplot actual mean
  oplot, [t.meandenscont, t.meandenscont], [0, maxhist], $
         color=!blue,thick=!p.thick*4
  oplot, [mean_mdenscont, mean_mdenscont], [0, maxhist], $
         thick=!p.thick;*2

  IF nw2 NE 0 THEN BEGIN 
      plothist, mdenscont2, xhist2,yhist2,bin=bin, /overplot, /norm, color=!red

      oplot, [mean_mdenscont2, mean_mdenscont2], [0, maxhist], $
             color=!red,thick=!p.thick ;*2
  ENDIF 

  ;; what decimal to round?
  rndp = 5
  ;; how many string positions to keep?
  IF mean_mdenscont GE 10 THEN BEGIN 
      strp = 8
  ENDIF ELSE strp = 7
  IF nw2 NE 0 THEN BEGIN 
      IF mean_mdenscont2 GE 10 THEN BEGIN 
          strp2 = 8
      ENDIF ELSE strp2 = 7
  ENDIF 
  ;; for sig
  rnds = 2
  strs = 4
  ;; string positions to keep for error
  stre = 7

  IF !d.name EQ 'X' THEN lchar = 1.0 ELSE lchar=0.7
  IF nw2 NE 0 THEN BEGIN 
      legend, ['Photoz: '+!csym.sigma+' = '+ntostr(rnd(sig_mdenscont, rnds),strs)+' mean = '+ntostr(rnd(mean_mdenscont, rndp),strp )+!csym.plusminus+ntostr(rnd(err_mdenscont, rndp),stre), $
               'All:    '+!csym.sigma+' = '+ntostr(rnd(sig_mdenscont2,rnds),strs)+' mean = '+ntostr(rnd(mean_mdenscont2,rndp),strp2)+!csym.plusminus+ntostr(rnd(err_mdenscont2,rndp),stre), $
               'True: '+ntostr(rnd(t.meandenscont,rndp),strp)], $
              line=[0,0,0], color=[!p.color, !red, !blue], $
              thick=[!p.thick,!p.thick,!p.thick],$
              /clear,charsize=lchar
  ENDIF ELSE BEGIN 
      legend, ['mean = '+ntostr(rnd(mean_mdenscont, rndp),strp )+!csym.plusminus+ntostr(rnd(err_mdenscont, rndp),stre), $
               'True: '+ntostr(rnd(t.meandenscont,rndp),strp)], $
              line=[0,0], color=[!p.color, !blue], $
              thick=[!p.thick,!p.thick],$
              /clear,charsize=lchar


      legend,!csym.sigma+' = '+ntostr(rnd(sig_mdenscont, rnds),strs),$
             /right,/clear,charsize=1
  ENDELSE 
  ;; Fit a gaussian to the histrogram
  ;; objects not weighted in making the histogram
  
  ;; zcut
  w=where(yhist NE 0,nw)
  tyfit = gaussfit(xhist[w], yhist[w], A, nterms=3)
  vals = (xhist[w] - mean_mdenscont)^2/2./sig_mdenscont^2 < 9.0
  tyfit2 = A[0]*exp(-vals)
  oplot, xhist[w], tyfit2, line=2
  ;;oplot, xhist[w], tyfit, line=2

  print
  print,'Straight: sigma = '+ntostr(sig_mdenscont)+' mean = '+$
        ntostr(mean_mdenscont)+!plusminus+ntostr(err_mdenscont)
  nsigoff = abs( (mean_mdenscont-t.meandenscont)/err_mdenscont )
  print,'Sigma2 = ',A[2]
  print,'Nsig off: '+ntostr(nsigoff)

  ;; all
  IF nw2 NE 0 THEN BEGIN 
      w2=where(yhist2 NE 0,nw2)
      tyfit = gaussfit(xhist2[w2], yhist2[w2], A2, nterms=3)
      vals = (xhist2[w2] - mean_mdenscont2)^2/2./sig_mdenscont2^2 < 9.0
      tyfit = A2[0]*exp(-vals)
      oplot, xhist2[w2], tyfit, color=!red

      print,'Straight: sigma = '+ntostr(sig_mdenscont2)+' mean = '+$
            ntostr(mean_mdenscont2)+!plusminus+ntostr(err_mdenscont2)
      nsigoff = abs( (mean_mdenscont2-t.meandenscont)/err_mdenscont2 )
      print,'Nsig off: '+ntostr(nsigoff)
      print,'Sigma2 = ',A2[2]
  ENDIF 
  print

  t=create_struct('zL', t.zL[0], $
                  'nlens', t.nlens,$
                  'mean_mdensconttrue', t.meandenscont, $
                  'mean_mdenscont', mean_mdenscont, $
                  'sig_mdenscont', sig_mdenscont, $
                  'err_mdenscont', err_mdenscont, $
                  'denscont', t.meandenscont)
  IF nw2 NE 0 THEN BEGIN 

      t=create_struct(t, $
                      'mean_mdenscont2', mean_mdenscont2, $
                      'sig_mdenscont2', sig_mdenscont2, $
                      'err_mdenscont2', err_mdenscont2)
                      
  ENDIF 

  print,'Output fits file: ',fitfile
  mwrfits, t, fitfile, /create

  !p.multi=0

  IF keyword_set(doplot) THEN endplot

END 
