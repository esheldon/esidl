PRO beta_sigmav_chisq_fixsig, wclr, nops=nops,paper=paper, sigonly=sigonly, noweight=noweight

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: beta_sigmav_chisq_fixsig, lum_color, nops=nops,paper=paper'
      print
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; This plots the chisq surface for all bands for
  ;; High, Low, and All density regions
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; fix sigma_star ~ 117 km/s, or wsig = 74

  wsig = 47

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
      outfile1 = outdir + 'conv_beta_fits_'+colors[wclr]+'w_'+typestr[type]+wstr+'_fixsig_paper.eps'
      outfile2 = outdir + 'conv_betadenscont_'+colors[wclr]+'w_'+typestr[type]+wstr+'_fixsig_paper.eps'
      noplotmin=1
  ENDIF ELSE BEGIN 
      outfile1 = outdir + 'conv_beta_fits_'+colors[wclr]+'w_'+typestr[type]+wstr+'_fixsig.eps'
      outfile2 = outdir + 'conv_betadenscont_'+colors[wclr]+'w_'+typestr[type]+wstr+'_fixsig.eps'
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

;  yt=!tsym.chi+'!U2!N/'+!tsym.nu
  yt=!tsym.delta_cap+!tsym.chi+'!U2!N'

;  tt=['High Density Regions','Low Density Regions','All Galaxies']
  tt=replicate(colorsp[wclr]+' lum',3)
  IF NOT keyword_set(paper) THEN tt=tt+'   95% conf. levels.'

  ;; read model; they are same for lensing in different bandpasses,
  ;; use red (2)
  read_wthetaconv_func_lumw, 2, wclr, beta, sigma, modelx, model,$
    neigh, central, central_dens, model_ratio, model_ratio_cum, sigonly=sigonly, noweight=noweight
;(
  model = reform(model[*,wsig,*])

;IF type EQ 2 THEN model=model*2.88/0.92

  duse=[0,1,2,3]
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
;datag=(datar+datai)/2.0
          chisq_conf_1par_3band, dataxg, datag, dataerrg, $
            dataxr, datar, dataerrr, $
            dataxi, datai, dataerri, $
            gg, rr, ii, gr, gi, ri, $
            modelx, model, beta, $
            chisq_surf, min1, low1, high1, $
            xtitle=xtitle, title=title,$
            xrange=xrange, xstyle=1,$
            noplotmin=noplotmin,charsize=charsize,$
            xtitle=xtitle,ytitle=ytitle,title=title,$
            wparam=wparam
      ENDIF ELSE BEGIN
          addstr = colorsp[clr[i]]
          infile = datadir + typestr[type]+'_zgal_gal'+wstr+'_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_'+colors[clr[i]]+'_corr_N1.fit'
;          IF clr[i] EQ 1 THEN infile = datadir + typestr[type]+'_zgal_gal_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_comb_corr_N1.fit'
          corrdata = mrdfits(infile,1,/silent)
          
          w=where(corrdata.meanr NE 0,nw)
          datax = corrdata.meanr[w]
          data = corrdata.sigma[w]
          dataerr = corrdata.sigmaerr[w]

          chisq_conf_1par, datax, data, dataerr, modelx, model, $
            beta, $
            chisq_surf, min1, low1, high1, $
            xtitle=xtitle, title=title,noplotmin=noplotmin,$
            xrange=xrange, xstyle=1,charsize=charsize,$
            xtitle=xtitle,ytitle=ytitle,title=title,$
            wparam=wparam
      ENDELSE 

;(
;      nyticks=7
;      ytickn=[' ','120','140','160','180','200',' ']
;      axis, yaxis=0,yticks=nyticks-1,ytickn=ytickn, ytitle=ytitle
;      axis, yaxis=1, yticks=nyticks-1, ytickn=[replicate(' ',nyticks)]
      
      IF NOT keyword_set(paper) THEN BEGIN 
          message = strarr(1)
          message[0] = addstr+' '+!tsym.beta+' '+'['+$
            ntostr(rnd(low1[1],2),4)+', '+$
            ntostr(rnd(min1,2),4)+', '+$
            ntostr(rnd(high1[1],2),4)+']'
          legend, message,/right,box=0,charsize=1.2
      ENDIF ELSE BEGIN 
          message = addstr
          legend, message,/right,box=0
      ENDELSE 

      print,'beta'
      print,min1
      forprint,low1,high1

      CASE i OF
          0: BEGIN
              gdatax=datax
              gdata=data
              gdataerr=dataerr
              gwcut=wparam
          END 
          1: BEGIN
              rdatax=datax
              rdata=data
              rdataerr=dataerr
              rwcut=wparam
          END 
          2: BEGIN
              idatax=datax
              idata=data
              idataerr=dataerr
              iwcut=wparam
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
  tmpmod=reform(model[gwcut,*])
  tmpcent = reform(central[gwcut,wsig,*])
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
  tmpmod=reform(model[rwcut,*])
  tmpcent = reform(central[rwcut,wsig,*])
  ploterror,rdatax, rdata, rdataerr, psym=1, yrange=yrange,$
    ystyle=1,charsize=charsize
  oplot,modelx[w],tmpmod[w]
  oplot,modelx[w],tmpcent[w],line=2
  oplot,[0,10000],[0,0]
  legend,'r'+!tsym.prime,/right,box=0

  multiplot
  tmpmod=reform(model[iwcut,*])
  tmpcent = reform(central[iwcut,wsig,*])
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
