PRO low_high_chisq_powerlaw, type, nops=nops,paper=paper

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
  typestr = ['high','low', 'main']

  IF keyword_set(paper) THEN BEGIN 
      outfile1 = outdir + 'fits_power_'+typestr[type]+'_paper.eps'
      outfile2 = outdir + 'denscont_power_'+typestr[type]+'_paper.eps'
  ENDIF ELSE BEGIN 
      outfile1 = outdir + 'fits_power_'+typestr[type]+'.eps'
      outfile2 = outdir + 'denscont_power_'+typestr[type]+'.eps'
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

  xt=!tsym.alpha
  yt='A (h M'+sunsymbol()+' pc!U'+!tsym.minus+'2!N)'
;  tt=['High Density Regions','Low Density Regions','All Galaxies']
  tt=['','','']
  IF NOT keyword_set(paper) THEN tt=tt+'   95% conf. levels.'

  npow=400
  nnorm=400
  power=arrscl( findgen(npow), -1.3, -0.3 )
  norm = arrscl( findgen(nnorm), 0, 6 )

;IF type EQ 2 THEN model=model*2.88/0.92

  yrange=[min(norm),max(norm)]
  xrange=[min(power),max(power)]

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

          infileg=datadir + typestr[type]+'_zgal_gal_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_g_corr_N1.fit'
          corrdata = mrdfits(infileg,1)
          w=where(corrdata.meanr NE 0,nw)
          dataxg = corrdata.meanr[w]
          datag = corrdata.sigma[w]
          dataerrg = corrdata.sigmaerr[w]

          infiler=datadir + typestr[type]+'_zgal_gal_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_r_corr_N1.fit'
          corrdata=mrdfits(infiler,1,/silent)
          w=where(corrdata.meanr NE 0,nw)
          dataxr = corrdata.meanr[w]
          datar = corrdata.sigma[w]
          dataerrr = corrdata.sigmaerr[w]

          infilei=datadir + typestr[type]+'_zgal_gal_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_i_corr_N1.fit'
          corrdata=mrdfits(infilei,1,/silent)
          w=where(corrdata.meanr NE 0,nw)
          dataxi = corrdata.meanr[w]
          datai = corrdata.sigma[w]
          dataerri = corrdata.sigmaerr[w]

          infilecomb=datadir + typestr[type]+'_zgal_gal_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_comb_corr_N1.fit'
          corrdata=mrdfits(infilecomb,1,/silent)
          w=where(corrdata.meanr NE 0,nw)
          combdatax = corrdata.meanr[w]
          combdata = corrdata.sigma[w]
          combdataerr = corrdata.sigmaerr[w]


          pow_chisq_conf_3band, dataxg/1000., datag, dataerrg, $
            dataxr, datar, dataerrr, $
            dataxi, datai, dataerri, $
            gg, rr, ii, gr, gi, ri, $
            power, norm, $
            chisq_surf, min1, min2, low1, high1, low2, high2, $
            xtitle=xtitle, title=title,$
            yrange=yrange,xrange=xrange, xstyle=1,$
            charsize=charsize,$
            xtitle=xtitle,ytitle=ytitle,title=title,$
            wparam1=wparam1,wparam2=wparam2;,yrange=[50,180],ystyle=1
      ENDIF ELSE BEGIN
          addstr = colorsp[clr[i]]
          infile = datadir + typestr[type]+'_zgal_gal_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_'+colors[clr[i]]+'_corr_N1.fit'
          corrdata = mrdfits(infile,1,/silent)
          
          w=where(corrdata.meanr NE 0,nw)
          datax = corrdata.meanr[w]
          data = corrdata.sigma[w]
          dataerr = corrdata.sigmaerr[w]

          pow_chisq_conf, datax/1000., data, dataerr, power, norm, $
            chisq_surf, min1, min2, low1, high1, low2, high2, $
            xtitle=xtitle, title=title,$
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
          message[0] = addstr+'A  '+'['+$
            ntostr(low2[1],5)+', '+$
            ntostr(min2,5)+', '+$
            ntostr( high2[1],5)+']'
          message[1] = !tsym.alpha+' '+'['+$
            ntostr(low1[1],6)+', '+$
            ntostr(min1,6)+', '+$
            ntostr(high1[1],6)+']'
          legend, message,/right,box=0,charsize=1.2
      ENDIF ELSE BEGIN 
          message = addstr
          legend, message,/right,box=0
      ENDELSE 


      print,'Power'
      print,min1
      forprint,low1,high1
      print,'Norm'
      print,min2
      forprint,low2,high2

      CASE i OF
          0: BEGIN
              gdatax=datax
              gdata=data
              gdataerr=dataerr
              gpow=min1
              gnorm=min2
          END 
          1: BEGIN
              rdatax=datax
              rdata=data
              rdataerr=dataerr
              rpow=min1
              rnorm=min2
          END 
          2: BEGIN
              idatax=datax
              idata=data
              idataerr=dataerr
              ipow=min1
              inorm=min2
          END 
          3: BEGIN 
              combpow=min1
              combnorm=min2
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
    '!N (R) (h M'+sunsymbol()+' pc!U'+!tsym.minus+'2!N)!X'
  xtitle='Projected Radius (h!U'+!tsym.minus+'1!N kpc)'
  
  yrange = [-9.5,40]
  pos=[400,37]

  xx=arrscl(findgen(1000), min(gdatax)/1000., max(gdatax)/1000.)

;  erase & multiplot, [1,3]
  erase & multiplot, [1,4]

  model = gnorm*xx^gpow
  ploterror,gdatax, gdata, gdataerr, psym=1, yrange=yrange,$
    ystyle=1,charsize=charsize,title=tt[type]
  oplot,xx*1000.,model
  oplot,[0,10000],[0,0]
  legend,'g'+!tsym.prime,/right,box=0

  ;; lable curves
;  legend,['Best Fit','Central Galaxy'],line=[0,2],charsize=1.0,$
;    pos=pos,box=0,thick=[!p.thick,!p.thick]

  multiplot

  model = rnorm*xx^rpow
  ploterror,rdatax, rdata, rdataerr, psym=1, yrange=yrange,$
    ystyle=1,charsize=charsize
  oplot,xx*1000.,model
  oplot,[0,10000],[0,0]
  legend,'r'+!tsym.prime,/right,box=0

  multiplot

  model = inorm*xx^ipow
  ploterror,idatax, idata, idataerr, psym=1, yrange=yrange,$
    ystyle=1,charsize=charsize;,xtitle=xtitle,ytitle=ytitle
   oplot,xx*1000.,model
  oplot,[0,10000],[0,0]
  legend,'i'+!tsym.prime,/right,box=0

  multiplot

  model = combnorm*xx^combpow
  ploterror,combdatax, combdata, combdataerr*errfac, psym=1, yrange=yrange,$
    ytitle=ytitle,xtitle=xtitle,ystyle=1,charsize=charsize
   oplot,xx*1000.,model
  oplot,[0,10000],[0,0]
  legend,addstr,/right,box=0

  multiplot,/reset

  IF NOT keyword_set(nops) THEN endplot


  !p.charsize=0

return
END 
