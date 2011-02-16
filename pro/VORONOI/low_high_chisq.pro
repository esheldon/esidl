PRO low_high_chisq, type, nops=nops,paper=paper, logbin=logbin

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: low_high_chisq, type, nops=nops,paper=paper'
      print
      print,' type=0 (high dens) type=1 (low dens) type=2 (tot)'
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
  ;typestr = ['high','low', 'tot']
  typestr = ['high','low','main']

  IF keyword_set(paper) THEN BEGIN 
      outfile1 = outdir + 'conv_isotrunc_fits_'+typestr[type]+'_paper.eps'
      outfile2 = outdir + 'conv_denscont_'+typestr[type]+'_paper.eps'
  ENDIF ELSE BEGIN 
      outfile1 = outdir + 'conv_isotrunc_fits_'+typestr[type]+'.eps'
      outfile2 = outdir + 'conv_denscont_'+typestr[type]+'.eps'
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


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; All regions
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  datadir = '/sdss5/data0/lensout/stripe10/'

  nclr=nclr+1                   ;extra for combined
  erase & multiplot, [1, nclr]

  xt='Cutoff Radius (h!U'+!tsym.minus+'1!N kpc)'
  yt=!tsym.sigma+'!DV!N (km/s)'
  tt=['High Density Regions','Low Density Regions','All Galaxies']

  clevel=0
  clevstr = ['   68% conf. levels.', '   95% conf. levels.', '   99% conf. levels.']
  IF NOT keyword_set(paper) THEN tt=tt+clevstr[clevel]

  ;; read model
  read_wthetaconv_func, type, cutoff, sigma, modelx, model,$
    neigh, central, central_dens, model_ratio, model_ratio_cum

;IF type EQ 2 THEN model=model*2.88/0.92

  yrange=[min(sigma),max(sigma)]
  xrange=[min(cutoff),max(cutoff)]

  ;;stripestr = 'stripe36_stripe37_stripe42_stripe43_stripe82'
  stripestr = 'stripe10'
  ext = 'N6.fit'

  FOR i=0L, nclr-1 DO BEGIN

      xtitle=''
      ytitle=''
      title=''
      IF i EQ nclr-1 THEN BEGIN
          xtitle=xt
          ytitle=yt
          title=''
      ENDIF ELSE IF i EQ 0 THEN BEGIN
          title=tt[type]
      ENDIF 

      ;; read data
      IF i EQ nclr-1 THEN BEGIN 

          addstr = colorsp[1]+'+'+colorsp[2]+'+'+colorsp[3]

          infileg=datadir + typestr[type]+'_zgal_gal_'+stripestr+'_g_corr_'+ext
          gcorrdata = mrdfits(infileg,1)
          w=where(gcorrdata.meanr NE 0,nw)
          dataxg = gcorrdata.meanr[w]
          datag = gcorrdata.sigma[w]
          dataerrg = gcorrdata.sigmaerr[w]

          infiler=datadir + typestr[type]+'_zgal_gal_'+stripestr+'_r_corr_'+ext
          rcorrdata=mrdfits(infiler,1,/silent)
          w=where(rcorrdata.meanr NE 0,nw)
          dataxr = rcorrdata.meanr[w]
          datar = rcorrdata.sigma[w]
          dataerrr = rcorrdata.sigmaerr[w]

          infilei=datadir + typestr[type]+'_zgal_gal_'+stripestr+'_i_corr_'+ext
          icorrdata=mrdfits(infilei,1,/silent)
          w=where(icorrdata.meanr NE 0,nw)
          dataxi = icorrdata.meanr[w]
          datai = icorrdata.sigma[w]
          dataerri = icorrdata.sigmaerr[w]

          infilecomb=datadir + typestr[type]+'_zgal_gal_'+stripestr+'_comb_corr_'+ext
          ccorrdata=mrdfits(infilecomb,1,/silent)
          w=where(ccorrdata.meanr NE 0,nw)
          combdatax = ccorrdata.meanr[w]
          combdata = ccorrdata.sigma[w]
          combdataerr = ccorrdata.sigmaerr[w]

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
            wparam1=wparam1,wparam2=wparam2;,yrange=[50,180],ystyle=1
      ENDIF ELSE BEGIN
          addstr = colorsp[clr[i]]
          infile = datadir + typestr[type]+'_zgal_gal_'+stripestr+'_'+!colors[clr[i]]+'_corr_'+ext
          corrdata = mrdfits(infile,1,/silent)
          
          w=where(corrdata.meanr NE 0,nw)
          datax = corrdata.meanr[w]
          data = corrdata.sigma[w]
          dataerr = corrdata.sigmaerr[w]

          chisq_conf, datax, data, dataerr, modelx, model, cutoff, sigma, $
            chisq_surf, min1, min2, low1, high1, low2, high2, $
            xtitle=xtitle, title=title,/noplotmin,$
            yrange=yrange,xrange=xrange, xstyle=1,charsize=charsize,$
            xtitle=xtitle,ytitle=ytitle,title=title,$
            wparam1=wparam1,wparam2=wparam2;,yrange=[50,180],ystyle=1
      ENDELSE 


;      nyticks=7
;      ytickn=[' ','120','140','160','180','200',' ']
;      axis, yaxis=0,yticks=nyticks-1,ytickn=ytickn, ytitle=ytitle
;      axis, yaxis=1, yticks=nyticks-1, ytickn=[replicate(' ',nyticks)]
      
      IF NOT keyword_set(paper) THEN BEGIN 
          message = strarr(2)
          message[0] = addstr+' '+!tsym.sigma+'!DV!N '+'['+$
            ntostr(long(rnd(low2[clevel])))+', '+$
            ntostr(long(rnd(min2)))+', '+$
            ntostr( long(rnd(high2[clevel])))+']'
          message[1] = 'Cutoff '+'['+$
            ntostr(long(rnd(low1[clevel])))+', '+$
            ntostr(long(rnd(min1)))+', '+$
            ntostr(long(rnd(high1[clevel])))+']'
          legend, message,/right,box=0,charsize=1.2
      ENDIF ELSE BEGIN 
          message = addstr
          legend, message,/right,box=0
      ENDELSE 


      ;xyouts, pos[0], pos[1], 'Cutoff '+'['$
      ;  +ntostr(low1[1])+', '+ntostr(high1[1])+']', charsize=0.8
      ;xyouts, pos[0], pos[1]+7., 'Vel. '+'['$
      ;  +ntostr(low2[1])+', '+ntostr(high2[1])+']', charsize=0.8
      ;xyouts, pos[0], pos[1]+14., '95% conf. levels.',charsize=0.8

      print,'cutoff'
      print,min1
      forprint,low1,high1
      print,'sigma'
      print,min2
      forprint,low2,high2

      CASE i OF
          0: BEGIN
              gdatax=datax
              gdata=data
              gdataerr=dataerr
              gwcut=wparam1
              gwsig=wparam2
          END 
          1: BEGIN
              rdatax=datax
              rdata=data
              rdataerr=dataerr
              rwcut=wparam1
              rwsig=wparam2
          END 
          2: BEGIN
              idatax=datax
              idata=data
              idataerr=dataerr
              iwcut=wparam1
              iwsig=wparam2
          END 
          3: BEGIN 
              combwcut=wparam1
              combwsig=wparam2
          END 
          ELSE:
      ENDCASE 
      IF i NE nclr-1 THEN multiplot
  ENDFOR 
  multiplot, /reset

  IF NOT nops THEN endplot

  IF NOT nops THEN BEGIN
      begplot, name=outfile2,/color,/encap
  ENDIF

  IF nops THEN key=get_kbrd(1)

  ytitle = '!S'+!tsym.sigma_cap+'!R!A'+!tsym.minus+'!N ('+!tsym.ltequal+'R) '+!tsym.minus+' !S'+$
    !tsym.sigma_cap+'!R!A'+!tsym.minus+$
    '!N (R) (h M!DSun!N pc!U'+!tsym.minus+'2!N)!X'
  xtitle='Projected Radius (h!U'+!tsym.minus+'1!N kpc)'

;  yrange1=prange(gdata,rdata,gdataerr,rdataerr)
;  yrange2=prange(gdata,idata, gdataerr,idataerr)
;  yrange=yrange1
;  yrange[0] = min([yrange1[0], yrange2[0]])
;  yrange[1] = max([yrange2[1], yrange2[1]])
;  yrange = [-9.5,40]

  yrange=[0.1,80]

  pos=[400,37]

;  erase & multiplot, [1,3]
  IF keyword_set(logbin) THEN square=1 ELSE square=0
  erase & multiplot, [2,2], square=square
  w=where(modelx LE max(gdatax) AND modelx GE min(gdatax))
  tmpmod=reform(model[gwcut,gwsig,*])
  tmpcent = reform(central[gwcut,gwsig,*])
  ;ploterror,gdatax, gdata, gdataerr, psym=1, yrange=yrange,$
  ;  ystyle=1,charsize=charsize,title=tt[type]
  plot_density_contrast, gcorrdata, logbin=logbin,/noxtitle,wuse=w,yrange=yrange
  oplot,modelx[w],tmpmod[w]
  oplot,modelx[w],tmpcent[w],line=2
  oplot,[0,10000],[0,0]
  legend,'g'+!tsym.prime,/right,box=0

  ;; lable curves
  legend,['Best Fit','Central Galaxy'],line=[0,2],charsize=1.0,$
    pos=pos,box=0,thick=[!p.thick,!p.thick]

  multiplot
  tmpmod=reform(model[rwcut,rwsig,*])
  tmpcent = reform(central[rwcut,rwsig,*])
;  ploterror,rdatax, rdata, rdataerr, psym=1, yrange=yrange,$
;    ystyle=1,charsize=charsize
  plot_density_contrast, rcorrdata, logbin=logbin,/noxtitle,/noytitle,wuse=w,yrange=yrange
  oplot,modelx[w],tmpmod[w]
  oplot,modelx[w],tmpcent[w],line=2
  oplot,[0,10000],[0,0]
  legend,'r'+!tsym.prime,/right,box=0

  multiplot
  tmpmod=reform(model[iwcut,iwsig,*])
  tmpcent = reform(central[iwcut,iwsig,*])
;  ploterror,idatax, idata, idataerr, psym=1, yrange=yrange,$
;    ystyle=1,charsize=charsize;,xtitle=xtitle,ytitle=ytitle
  plot_density_contrast, icorrdata, logbin=logbin,wuse=w,yrange=yrange
  oplot,modelx[w],tmpmod[w]
  oplot,modelx[w],tmpcent[w],line=2
  oplot,[0,10000],[0,0]
  legend,'i'+!tsym.prime,/right,box=0

  multiplot
  tmpmod=reform(model[combwcut,combwsig,*])
  tmpcent = reform(central[combwcut,combwsig,*])
;  ploterror,combdatax, combdata, combdataerr*errfac, psym=1, yrange=yrange,$
;    ytitle=ytitle,xtitle=xtitle,ystyle=1,charsize=charsize
  plot_density_contrast, ccorrdata, logbin=logbin,wuse=w,/noytitle,yrange=yrange
  oplot,modelx[w],tmpmod[w]
  oplot,modelx[w],tmpcent[w],line=2
  oplot,[0,10000],[0,0]
  legend,addstr,/right,box=0

  multiplot,/reset

  IF NOT keyword_set(nops) THEN endplot


  !p.charsize=0

return
END 
