PRO beta_sigmav_chisq, wclr, nops=nops,paper=paper, sigonly=sigonly, noweight=noweight

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: beta_sigmav_chisq, lum_color, nops=nops,paper=paper'
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

  outdir='/sdss5/data0/lensout/wtheta_conv/'
  type=0
  typestr = ['lum','lum1', 'lum2','lum3','lum4']

  colors=['u','g','r','i','z']
  CASE 1 OF
      keyword_set(sigonly): wstr='_sigonly'
      keyword_set(noweight): wstr='_noweight'
      ELSE: wstr=''
  ENDCASE 
  IF keyword_set(paper) THEN BEGIN 
      outfile1 = outdir + 'conv_beta_fits_'+colors[wclr]+'w_'+typestr[type]+wstr+'_paper.eps'
      outfile2 = outdir + 'conv_betadenscont_'+colors[wclr]+'w_'+typestr[type]+wstr+'_paper.eps'
      noplotmin=1
  ENDIF ELSE BEGIN 
      outfile1 = outdir + 'conv_beta_fits_'+colors[wclr]+'w_'+typestr[type]+wstr+'.eps'
      outfile2 = outdir + 'conv_betadenscont_'+colors[wclr]+'w_'+typestr[type]+wstr+'.eps'
  ENDELSE 
  IF NOT keyword_set(nops) THEN nops = 0
  IF NOT nops THEN begplot, name=outfile1,/color;,/encap

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

  datadir = '/sdss5/data0/lensout/stripe10/sublum/'+colors[wclr]+'/'

  nclr=nclr+1                   ;extra for combined
  erase & multiplot, [1, nclr]

;  xt='Cutoff Radius (h!U'+!tsym.minus+'1!N kpc)'
  xt = !tsym.beta
  yt=!tsym.sigma+'!D*!N (km/s)'
;  tt=['High Density Regions','Low Density Regions','All Galaxies']
  tt=replicate(colorsp[wclr]+' lum',3)
  IF NOT keyword_set(paper) THEN tt=tt+'   95% conf. levels.'

  ;; read model; they are same for lensing in different bandpasses,
  ;; use red (2)
  read_wthetaconv_func_lumw, 2, wclr, beta, sigma, modelx, model,$
    neigh, central, central_dens, model_ratio, model_ratio_cum, sigonly=sigonly, noweight=noweight

;IF type EQ 2 THEN model=model*2.88/0.92

  yrange=[min(sigma),max(sigma)]
  xrange=[min(beta),max(beta)]

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
          addstr = 'comb  '
          addstr = colorsp[1]+'+'+colorsp[2]+'+'+colorsp[3]

          infileg=datadir + typestr[type]+'_zgal_gal'+wstr+'_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_g_corr_N1.fit'
          corrdata = mrdfits(infileg,1)
          w=where(corrdata.meanr NE 0,nw)
          dataxg = corrdata.meanr[w]
          datag = corrdata.sigma[w]
          dataerrg = corrdata.sigmaerr[w]

          infiler=datadir + typestr[type]+'_zgal_gal'+wstr+'_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_r_corr_N1.fit'
          corrdata=mrdfits(infiler,1,/silent)
          w=where(corrdata.meanr NE 0,nw)
          dataxr = corrdata.meanr[w]
          datar = corrdata.sigma[w]
          dataerrr = corrdata.sigmaerr[w]

          infilei=datadir + typestr[type]+'_zgal_gal'+wstr+'_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_i_corr_N1.fit'
          corrdata=mrdfits(infilei,1,/silent)
          w=where(corrdata.meanr NE 0,nw)
          dataxi = corrdata.meanr[w]
          datai = corrdata.sigma[w]
          dataerri = corrdata.sigmaerr[w]

          chisq_conf_3band, dataxg, datag, dataerrg, $
            dataxr, datar, dataerrr, $
            dataxi, datai, dataerri, $
            gg, rr, ii, gr, gi, ri, $
            modelx, model, beta, sigma, $
            chisq_surf, min1, min2, low1, high1, low2, high2, $
            xtitle=xtitle, title=title,$
            yrange=yrange,xrange=xrange, xstyle=1,$
            noplotmin=noplotmin,charsize=charsize,$
            xtitle=xtitle,ytitle=ytitle,title=title,$
            wparam1=wparam1,wparam2=wparam2;,yrange=[50,180],ystyle=1
      ENDIF ELSE BEGIN
          addstr = colorsp[clr[i]]
          infile = datadir + typestr[type]+'_zgal_gal'+wstr+'_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_'+colors[clr[i]]+'_corr_N1.fit'
          corrdata = mrdfits(infile,1,/silent)
          
          w=where(corrdata.meanr NE 0,nw)
          datax = corrdata.meanr[w]
          data = corrdata.sigma[w]
          dataerr = corrdata.sigmaerr[w]

          chisq_conf, datax, data, dataerr, modelx, model, beta, sigma, $
            chisq_surf, min1, min2, low1, high1, low2, high2, $
            xtitle=xtitle, title=title,noplotmin=noplotmin,$
            yrange=yrange,xrange=xrange, xstyle=1,charsize=charsize,$
            xtitle=xtitle,ytitle=ytitle,title=title,$
            wparam1=wparam1,wparam2=wparam2;,yrange=[50,180],ystyle=1
      ENDELSE 

;(
;      nyticks=7
;      ytickn=[' ','120','140','160','180','200',' ']
;      axis, yaxis=0,yticks=nyticks-1,ytickn=ytickn, ytitle=ytitle
;      axis, yaxis=1, yticks=nyticks-1, ytickn=[replicate(' ',nyticks)]
      
      IF NOT keyword_set(paper) THEN BEGIN 
          message = strarr(2)
          message[0] = addstr+' '+!tsym.sigma+'!D*!N '+'['+$
            ntostr(long(rnd(low2[1])))+', '+$
            ntostr(long(rnd(min2)))+', '+$
            ntostr( long(rnd(high2[1])))+']'
          message[1] = !tsym.beta+' '+'['+$
            ntostr(rnd(low1[1],2),4)+', '+$
            ntostr(rnd(min1,2),4)+', '+$
            ntostr(rnd(high1[1],2),4)+']'
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

      print,'beta'
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

  
  yrange = [-9.5,40]
  pos=[400,37]

  erase & multiplot, [1,3]
  w=where(modelx LE max(gdatax) AND modelx GE min(gdatax))
  tmpmod=reform(model[gwcut,gwsig,*])
  tmpcent = reform(central[gwcut,gwsig,*])
  ploterror,gdatax, gdata, gdataerr, psym=1, yrange=yrange,$
    ystyle=1,charsize=charsize,title=tt[type]
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
  ploterror,rdatax, rdata, rdataerr, psym=1, yrange=yrange,$
    ystyle=1,charsize=charsize
  oplot,modelx[w],tmpmod[w]
  oplot,modelx[w],tmpcent[w],line=2
  oplot,[0,10000],[0,0]
  legend,'r'+!tsym.prime,/right,box=0

  multiplot
  tmpmod=reform(model[iwcut,iwsig,*])
  tmpcent = reform(central[iwcut,iwsig,*])
  ploterror,idatax, idata, idataerr, psym=1, yrange=yrange,$
    ytitle=ytitle,xtitle=xtitle,ystyle=1,charsize=charsize
  oplot,modelx[w],tmpmod[w]
  oplot,modelx[w],tmpcent[w],line=2
  oplot,[0,10000],[0,0]
  legend,'i'+!tsym.prime,/right,box=0

  multiplot,/reset

  IF NOT keyword_set(nops) THEN endplot


  !p.charsize=0

return
END 
