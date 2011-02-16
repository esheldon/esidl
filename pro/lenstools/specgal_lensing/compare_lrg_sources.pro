
PRO compare_lrg_sources, type, all, lrgs

  ;; Compare the results for lrg_sources and all

  setup_mystuff

  CASE type OF
      'eclass1_': yrange = [0,50]
      'match_': yrange = [0,40]
      ELSE: message,'Unknown type: ',type
  ENDCASE 

  indir= '/net/cheops1/data3/lensout/combstripe/comb/'
  all = $
    mrdfits(indir+$
            type + 'zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_comb_N2.fit',1)
  lrgs = $
    mrdfits(indir+$
            type + 'zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_lrg_recorr_h_comb_N1.fit',1)

  nbin = n_elements(all.meanr)
;  allfrac = all.nlbin/all.nlenses
;  lrgsfrac = lrgs.nlbin/lrgs.nlenses
;  minfrac = 0.95
;  w = where(allfrac GE minfrac AND $
;            lrgsfrac GE minfrac, nw)
  w=lindgen(nbin)

  !p.multi=[0,0,2]

  yrange = [0,50]
  xrange = [0, 1100]
  xtitle = !kpcxtitle2
  ytitle = !deltaytitle
  IF !d.name EQ 'PS' THEN lclr = !blue ELSE lclr = !green
  lline = 2

  chisq = total( (all.sigma-lrgs.sigma)^2/(all.sigmaerr^2 + lrgs.sigmaerr^2))
  print,chisq

  aploterror, !gratio, all.meanr[w], all.sigma[w], all.sigmaerr[w], $
    xrange=xrange, yrange=yrange, xtitle=xtitle, ytitle=ytitle, /center,$
    psym=8
  oplot, all.meanr[w], all.sigma[w]
  oploterror, lrgs.meanr[w], lrgs.sigma[w], lrgs.sigmaerr[w], $
    psym=4, color=lclr, errc=lclr
  oplot, lrgs.meanr[w], lrgs.sigma[w], line=lline, color=lclr

  oplot, [0,10000], [0,0]
  legend,$
    ['All sources', 'LRG sources'],$
    line=[0, lline],$
    color=[!p.color,lclr], /right,$
    thick=[!p.thick,!p.thick], box=0

  ratio = lrgs.sigma[w]/all.sigma[w]
  ratioerr = ratio*sqrt( (lrgs.sigmaerr[w]/lrgs.sigma[w])^2 + $
                         (all.sigmaerr[w]/all.sigma[w])^2 )
  
  rytitle = 'LRG/ALL'
  aploterror, !gratio, all.meanr[w], ratio, ratioerr, $
    xrange=xrange, xtitle=xtitle, ytitle=rytitle, /center,psym=8, $
    yrange=[0,2]
  oplot, [0,10000],[1,1]

  key = prompt_kbrd()

  ;; Cumulative plots

;  wmom, all.sigma[w], all.sigmaerr[w], all_tsigma, blah, all_tsigmaerr
;  wmom, lrgs.sigma[w], lrgs.sigmaerr[w], lrgs_tsigma, blah, lrgs_tsigmaerr

;  print,lrgs_tsigma/all_tsigma
;  print,lrgs.tsigma/all.tsigma

  ytitle = 'Cumulative '+ytitle
  aploterror, !gratio, all.rmax_act[w], all.tsigma[w], all.tsigmaerr[w], $
    xrange=xrange, yrange=yrange, xtitle=xtitle, ytitle=ytitle, /center,$
    psym=8
  oplot, all.rmax_act[w], all.tsigma[w]
  oploterror, lrgs.rmax_act[w], lrgs.tsigma[w], lrgs.tsigmaerr[w], $
    psym=4, color=lclr, errc=lclr
  oplot, lrgs.rmax_act[w], lrgs.tsigma[w], line=lline, color=lclr

  oplot, [0,10000], [0,0]
  legend,$
    ['All sources', 'LRG sources'],$
    line=[0, lline],$
    color=[!p.color,lclr], /right,$
    thick=[!p.thick,!p.thick],box=0

  ratio = lrgs.tsigma[w]/all.tsigma[w]
  ratioerr = ratio*sqrt( (lrgs.tsigmaerr[w]/lrgs.tsigma[w])^2 + $
                         (all.tsigmaerr[w]/all.tsigma[w])^2 )

  rytitle = 'LRG/ALL'
  aploterror, !gratio, all.rmax_act[w], ratio, ratioerr, $
    xrange=xrange, xtitle=xtitle, ytitle=rytitle, /center,psym=8, $
    yrange=[0,2]
  oplot, [0,10000],[1,1]

  !p.multi=0

END 

