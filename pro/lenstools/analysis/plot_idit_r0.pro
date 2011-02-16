PRO plot_idit_r0, nbin, new=new

  CASE nbin OF 
      2: BEGIN 
          lumstr = 'twobin'
      END 
      3: BEGIN 
          lumstr = 'threebin'
      END 
      4: BEGIN 
          lumstr = 'fourbin'
      END 
      ELSE: message,'What!'
  ENDCASE 


  lumclr=2

  combdir = $
    '/net/cheops2/home/esheldon/lensout/combstripe/comb/sublum/'+$
    !colors[lumclr]+'/'


  f1 = combdir + 'lum1'+lumstr+'_zgal_gal_stripe10_11_12_35_36_37_ri_jack_comb_N4.fit'
  f2 = combdir + 'lum2'+lumstr+'_zgal_gal_stripe10_11_12_35_36_37_ri_jack_comb_N4.fit'
  f3 = combdir + 'lum3'+lumstr+'_zgal_gal_stripe10_11_12_35_36_37_ri_jack_comb_N4.fit'
  f4 = combdir + 'lum4'+lumstr+'_zgal_gal_stripe10_11_12_35_36_37_ri_jack_comb_N4.fit'

  t1=mrdfits(f1,1) 
  add_arrval, t1.r0, r0
  add_arrval, t1.r0low, r0low
  add_arrval, t1.r0high, r0high
  add_arrval, t1.tmeanlum, meanlum

  t2=mrdfits(f2,1)
  add_arrval, t2.r0, r0
  add_arrval, t2.r0low, r0low
  add_arrval, t2.r0high, r0high
  add_arrval, t2.tmeanlum, meanlum

  IF nbin GT 2 THEN BEGIN
      t3=mrdfits(f3,1)
      add_arrval, t3.r0, r0
      add_arrval, t3.r0low, r0low
      add_arrval, t3.r0high, r0high
      add_arrval, t3.tmeanlum, meanlum
  ENDIF 
  IF nbin GT 3 THEN BEGIN 
      t4=mrdfits(f4,1)
      add_arrval, t4.r0, r0
      add_arrval, t4.r0low, r0low
      add_arrval, t4.r0high, r0high
      add_arrval, t4.tmeanlum, meanlum
  ENDIF 

  ;; new stuff
  IF keyword_set(new) THEN BEGIN 

      iabsmaglow = [-18., -19., -20., -21., -22.]
      iabsmaghigh = [-19., -20., -21., -22., -23.]
      iabsmag = (iabsmaghigh+iabsmaglow)/2.

      ir0 = [3.76, 5.14, 5.70, 6.32, 8.74]
      ir0err = [0.28, 0.39, 0.41, 0.19, 0.43]

  ENDIF ELSE BEGIN 
      iabsmaglow = [-20., -21.5, -23.]
      iabsmaghigh = [-18.5, -20.0, -21.5]
      iabsmag=(iabsmaghigh+iabsmaglow)/2.
      
      ir0 = [4.72,6.28,7.42]
      ir0err = [0.44,0.77,0.33]
  ENDELSE 

  xrange = [-17, -24]
  xtitle = 'M!Dr!N - 5 log!D10!Nh'
  ytitle = 'r!D0!N [h!U'+!csym.minus+'1!N Mpc]'
  plotderror,iabsmag,ir0,iabsmaglow,iabsmaghigh,ir0-ir0err,ir0+ir0err,psym=8,$
             xrange=xrange, /xsty, yrange=[1,10],xtitle=xtitle,ytitle=ytitle

  sun=[6.38,5.06,4.64,4.53,4.52]
  absmag = sun[lumclr] - 2.5*alog10(meanlum*1.e10)

  colprint,absmag,r0

  IF !d.name EQ 'X' THEN cc = !green ELSE cc = !blue
  oplotderror, absmag, r0, iabsmaglow, iabsmaghigh, r0low,r0high, psym=8,$
               color=cc,errc=cc

END 
