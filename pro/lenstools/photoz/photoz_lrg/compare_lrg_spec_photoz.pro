PRO compare_lrg_spec_photoz, lrg, south

  dir = sdssidl_config('shapecorr_dir') + 'combined/'
  IF n_elements(lrg) EQ 0 THEN BEGIN 
;      lrg=mrdfits(dir+'stripe82_lrg_vagc_photoz.fit',1)
      lrg=mrdfits(dir+'stripe09_10_11_12_13_14_15_16_26_27_28_29_30_31_32_33_34_35_36_37_42_43_44_76_82_86_lrg_vagc_photoz.fit',1)
  ENDIF 
  IF n_elements(south) EQ 0 THEN BEGIN 
      south=mrdfits(dir+'stripe82_southern_vagc_photoz.fit',1)
  ENDIF 

  photoz_dist_struct, lrg, lrg_z, lrg_pofz, wlrg
  photoz_dist_struct, south, south_z, south_pofz, wsouth

  binner_bynum, lrg[wlrg].z, lrg[wlrg].photoz_z, 2000, lxo, lyo, lsig, $
    yin_err = lrg[wlrg].photoz_zerr
  binner_bynum, south[wsouth].z, south[wsouth].photoz_z, 500, sxo, syo, ssig,$
    yin_err = south[wsouth].photoz_zerr

  rng = [0,0.65]
  tit = 'LRG'
  xtit='Spectro Z'
  ytit='Photo Z'
  botytit = 'Photo Z/Spectro Z'
;  plottwo,$
;    lrg[wlrg].z,lrg[wlrg].photoz_z,$
;    lrg[wlrg].z,lrg[wlrg].photoz_z/lrg[wlrg].z,$
;    xrange=rng,topyrange=rng,$
;    psym=3,/center,$
;    /ystyle,/xstyle,$
;    xtitle=xtit,$
;    topytitle=ytit,$
;    botytitle=botytit, $
;    title=tit,$
;    botyrange=[0.5,1.5],topaspect=1,frac=0.75, $
;    xoplot=[0,1],yoplot=[0,1],oplotcolor=!blue
;  oplot,[0,1],[1,1], color=!blue

 aplot, 1, lrg[wlrg].z,lrg[wlrg].photoz_z, $
    xrange=rng, yrange=rng, $
    psym=3, $
    /xstyle, /ystyle, $
    xtitle=xtit, ytitle=ytit, title=tit
  oplot, [0,1],[0,1],color=!blue
  oploterror, lxo, lyo, lsig,color=!red,errc=!red


  key=prompt_kbrd()

  IF !d.name EQ 'PS' THEN loadct,0
  ploth, lrg[wlrg].z,lrg[wlrg].photoz_z, $
    xrange=rng,yrange=rng, /xstyle, /ystyle, $
    xtitle=xtit,ytitle=ytit, title=tit

;  binner, lrg[wlrg].z, lrg[wlrg].photoz_z, 0.025, xo, yo, sig
  simpctable
  oploterror, lxo, lyo, lsig, color=!red,errc=!red
  oplot,[0,1],[0,1], color=!blue
  key=prompt_kbrd()

  tit = 'Stripe 82 LRG Special'
;  plottwo,$
;    south[wsouth].z,south[wsouth].photoz_z,$
;    south[wsouth].z,south[wsouth].photoz_z/south[wsouth].z,$
;    xrange=rng,topyrange=rng,$
;    psym=3,/center,$
;    /ystyle,/xstyle,$
;    xtitle=xtit,$
;    topytitle=ytit,$
;    botytitle=botytit,$
;    title=tit,$
;    botyrange=[0.5,1.5],topaspect=1,frac=0.75, $
;    xoplot=[0,1],yoplot=[0,1],oplotcolor=!blue
;  oplot,[0,1],[1,1], color=!blue
  
  aplot, 1, south[wsouth].z,south[wsouth].photoz_z, $
    xrange=rng, yrange=rng, $
    psym=3, $
    /xstyle, /ystyle, $
    xtitle=xtit, ytitle=ytit, title=tit
  oplot, [0,1],[0,1],color=!blue
  oploterror, sxo, syo, ssig,color=!red,errc=!red



  key=prompt_kbrd()

  IF !d.name EQ 'PS' THEN loadct,0
  ploth, south[wsouth].z,south[wsouth].photoz_z, $
    xrange=rng,yrange=rng, /xstyle, /ystyle, $
    xtitle=xtit,ytitle=ytit, title=tit

;  binner, south[wsouth].z, south[wsouth].photoz_z, 0.025, xo, yo, sig
  simpctable
  oploterror, sxo, syo, ssig,color=!red,errc=!red
  oplot,[0,1],[0,1], color=!blue
  key=prompt_kbrd()

  !p.multi=[0,0,2]

  plothist, lrg[wlrg].z, bin=0.01, /norm
  oplot, lrg_z, lrg_pofz
  legend,'lrg',/right,box=0

  plothist, south[wsouth].z, bin=0.01, /norm
  oplot, south_z, south_pofz
  legend,'Stripe 82 lrg special',/right,box=0
  !p.multi=0

END

