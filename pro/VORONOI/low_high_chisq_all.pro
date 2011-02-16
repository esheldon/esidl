PRO low_high_chisq_all, nops=nops,paper=paper, help=help

  IF keyword_set(help) THEN BEGIN 
      print,'-Syntax: low_high_chisq, nops=nops,paper=paper'
      print
      return
  ENDIF 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; This plots the chisq surface for all bands for
  ;; High, Low, and All density regions
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  gg = 3.50959599
  rr = 2.235299689
  ii = 2.549197961
  gr = 1.586466533 
  gi = 1.6983151 
  ri = 1.512359722

  errfac=1.32

  outdir='/sdss5/data0/lensout/wtheta_conv/'
  typestr = ['high','low','tot']
  types = [2, 0, 1]

  ntype = n_elements(typestr)

  IF keyword_set(paper) THEN BEGIN 
      outfile1 = outdir + 'conv_isotrunc_fits_all_paper.eps'
      outfile2 = outdir + 'conv_denscont_all_paper.eps'
  ENDIF ELSE BEGIN 
      outfile1 = outdir + 'conv_isotrunc_fits_all.eps'
      outfile2 = outdir + 'conv_denscont_all.eps'
  ENDELSE 
  IF NOT keyword_set(nops) THEN nops = 0
  IF NOT nops THEN begplot, name=outfile1,/color,/encap

  ;!p.charsize=1.2
  charsize=1.3
  colors=['u','g','r','i','z']
  colorsp=colors+!tsym.prime
  pos = [75., 105.]

  clr = [1,2,3]
  nclr=n_elements(clr)

  ;; read models
  type=0
  read_wthetaconv_func, type, cutoff, sigma, hmodelx, hmodel,$
    hneigh, hcentral, hcentral_dens, hmodel_ratio, hmodel_ratio_cum
  type=1
  read_wthetaconv_func, type, cutoff, sigma, lmodelx, lmodel,$
    lneigh, lcentral, lcentral_dens, lmodel_ratio, lmodel_ratio_cum
  type=2
  read_wthetaconv_func, type, cutoff, sigma, tmodelx, tmodel,$
    tneigh, tcentral, tcentral_dens, tmodel_ratio, tmodel_ratio_cum

  yrange=[min(sigma),max(sigma)]
  xrange=[min(cutoff),max(cutoff)]

  datadir = '/sdss5/data0/lensout/stripe10/'

  erase & multiplot, [1, ntype]

  xt='Cutoff Radius (h!U'+!tsym.minus+'1!N kpc)'
  yt=!tsym.sigma+'!DV!N (km/s)'
  tt=['High Density Regions','Low Density Regions','All Galaxies']
  IF NOT keyword_set(paper) THEN tt=tt+'   95% conf. levels.'

;IF type EQ 2 THEN model=model*2.88/0.92

  stripestr = 'stripe36_stripe37_stripe42_stripe43_stripe82'
  ext = 'N6.fit'

  
  FOR i=0L, ntype-1 DO BEGIN

      type = types[i]

      xtitle=''
      ytitle=''
      title=''
      IF i EQ ntype-1 THEN BEGIN
          xtitle=xt
          ytitle=yt
          title=''
      ENDIF ELSE IF i EQ 0 THEN BEGIN
          ;title=tt[type]
          title=''
      ENDIF 

      ;; read data

      addstr = colorsp[1]+'+'+colorsp[2]+'+'+colorsp[3]
      
      infileg=datadir + typestr[type]+'_zgal_gal_'+stripestr+'_g_corr_'+ext
      corrdata = mrdfits(infileg,1)
      w=where(corrdata.meanr NE 0,nw)
      dataxg = corrdata.meanr[w]
      datag = corrdata.sigma[w]
      dataerrg = corrdata.sigmaerr[w]
      
      infiler=datadir + typestr[type]+'_zgal_gal_'+stripestr+'_r_corr_'+ext
      corrdata=mrdfits(infiler,1,/silent)
      w=where(corrdata.meanr NE 0,nw)
      dataxr = corrdata.meanr[w]
      datar = corrdata.sigma[w]
      dataerrr = corrdata.sigmaerr[w]
      
      infilei=datadir + typestr[type]+'_zgal_gal_'+stripestr+'_i_corr_'+ext
      corrdata=mrdfits(infilei,1,/silent)
      w=where(corrdata.meanr NE 0,nw)
      dataxi = corrdata.meanr[w]
      datai = corrdata.sigma[w]
      dataerri = corrdata.sigmaerr[w]
      
      infilecomb=datadir + typestr[type]+'_zgal_gal_'+stripestr+'_comb_corr_'+ext
      corrdata=mrdfits(infilecomb,1,/silent)
      w=where(corrdata.meanr NE 0,nw)
      combdatax = corrdata.meanr[w]
      combdata = corrdata.sigma[w]
      combdataerr = corrdata.sigmaerr[w]

      CASE type OF
          0: BEGIN & modelx=hmodelx & model=hmodel & END
          1: BEGIN & modelx=lmodelx & model=lmodel & END
          2: BEGIN & modelx=tmodelx & model=tmodel & END
          ELSE:
      ENDCASE 

      chisq_conf_3band, dataxg, datag, dataerrg, $
        dataxr, datar, dataerrr, $
        dataxi, datai, dataerri, $
        gg, rr, ii, gr, gi, ri, $
        modelx, model, cutoff, sigma, $
        chisq_surf, min1, min2, low1, high1, low2, high2, $
        xtitle=xtitle, title=title,$
        yrange=yrange,xrange=xrange, xstyle=1,$
        /noplotmin,charsize=charsize,$
        xtitle=xtitle,ytitle=ytitle,title=title,$
        wparam1=wparam1,wparam2=wparam2 ;,yrange=[50,180],ystyle=1
            
      IF NOT keyword_set(paper) THEN BEGIN 
          message = strarr(3)
          message[0] = tt[type]
          message[1] = !tsym.sigma+'!DV!N '+'['+$
            ntostr(long(rnd(low2[1])))+', '+$
            ntostr(long(rnd(min2)))+', '+$
            ntostr( long(rnd(high2[1])))+']'
          message[2] = 'Cutoff '+'['+$
            ntostr(long(rnd(low1[1])))+', '+$
            ntostr(long(rnd(min1)))+', '+$
            ntostr(long(rnd(high1[1])))+']'
          legend, message,/right,box=0,charsize=1.2
      ENDIF ELSE BEGIN 
          message = tt[type]
          legend, message,/right,box=0
      ENDELSE 

      print,typestr[type]+ ' cutoff'
      print,min1
      forprint,low1,high1
      print,typestr[type]+ ' sigma'
      print,min2
      forprint,low2,high2

      CASE i OF
          0: BEGIN
              totdatax=combdatax
              totdata=combdata
              totdataerr=combdataerr
              totwcut=wparam1
              totwsig=wparam2
          END 
          1: BEGIN
              highdatax=combdatax
              highdata=combdata
              highdataerr=combdataerr
              highwcut=wparam1
              highwsig=wparam2
          END 
          2: BEGIN
              lowdatax=combdatax
              lowdata=combdata
              lowdataerr=combdataerr
              lowwcut=wparam1
              lowwsig=wparam2
          END 
          ELSE:
      ENDCASE 
      IF i NE ntype-1 THEN multiplot
  ENDFOR 
  multiplot, /reset

  IF NOT nops THEN endplot

  IF NOT nops THEN BEGIN
      begplot, name=outfile2,/color,/encap
  ENDIF

  IF nops THEN key=get_kbrd(1)

  ytitle = '!S'+!tsym.sigma_cap+'!R!A'+!tsym.minus+'!N ('+!tsym.ltequal+'R) '+!tsym.minus+' !S'+$
    !tsym.sigma_cap+'!R!A'+!tsym.minus+$
    '!N (R) (h M'+sunsymbol()+' pc!U'+!tsym.minus+'2!N)!X'
  xtitle='Projected Radius (h!U'+!tsym.minus+'1!N kpc)'
;  yrange1=prange(gdata,rdata,gdataerr,rdataerr)
;  yrange2=prange(gdata,idata, gdataerr,idataerr)
;  yrange=yrange1
;  yrange[0] = min([yrange1[0], yrange2[0]])
;  yrange[1] = max([yrange2[1], yrange2[1]])

  
  yrange = [-9.5,40]
  pos=[300,37]

;  erase & multiplot, [1,3]
  ;; all
  type=2
  erase & multiplot, [1,3]
  w=where(tmodelx LE max(totdatax) AND tmodelx GE min(totdatax))
  tmpmod=reform(tmodel[totwcut,totwsig,*])
  tmpcent = reform(tcentral[totwcut,totwsig,*])
  ploterror,totdatax, totdata, totdataerr*errfac, psym=1, yrange=yrange,$
    ystyle=1,charsize=charsize
  oplot,tmodelx[w],tmpmod[w]
  oplot,tmodelx[w],tmpcent[w],line=2
  oplot,[0,10000],[0,0]
  legend,tt[type],/right,box=0

  ;; lable curves
  legend,['Best Fit','Central Galaxy'],line=[0,2],charsize=1.0,$
    pos=pos,box=0,thick=[!p.thick,!p.thick]

  ;; High
  multiplot
  type=0
  tmpmod=reform(hmodel[highwcut,highwsig,*])
  tmpcent = reform(hcentral[highwcut,highwsig,*])
  ploterror,highdatax, highdata, highdataerr*errfac, psym=1, yrange=yrange,$
    ystyle=1,charsize=charsize
  oplot,hmodelx[w],tmpmod[w]
  oplot,hmodelx[w],tmpcent[w],line=2
  oplot,[0,10000],[0,0]
  legend,tt[type],/right,box=0

  ;; low
  ;; force to use tot cutoff?

  cdiff = abs(cutoff - 230.0)
  lowwcut= ( where( min(cdiff) EQ cdiff ) )[0]

  minc=min(chisq_surf[lowwcut, *])
  lowwsig=( where( chisq_surf[lowwcut,*] EQ minc ) )[0]

  type=1
  multiplot
  tmpmod=reform(lmodel[lowwcut,lowwsig,*])
  tmpcent = reform(lcentral[lowwcut,lowwsig,*])
  ploterror,lowdatax, lowdata, lowdataerr*errfac, psym=1, yrange=yrange,$
    ystyle=1,charsize=charsize,xtitle=xtitle,ytitle=ytitle
  oplot,lmodelx[w],tmpmod[w]
  oplot,lmodelx[w],tmpcent[w],line=2
  oplot,[0,10000],[0,0]
  legend,tt[type],/right,box=0

  multiplot,/reset

  IF NOT keyword_set(nops) THEN endplot


  !p.charsize=0

return
END 
