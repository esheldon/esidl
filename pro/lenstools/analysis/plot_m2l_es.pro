PRO plot_m2l_es_addstruct, file, lensum

  ;; make sure file was found (some stripes don't have any
  ;; matches in certain luminosity bins)
  lensumtmp = mrdfits(file,1,/silent)
  print,file
  IF datatype(lensumtmp) NE 'INT' THEN BEGIN 
      IF n_elements(lensum) EQ 0 THEN lensum = lensumtmp ELSE BEGIN
          concat_dstructs, lensum, lensumtmp, tmp
          lensum=0 & lensumtmp=0 & lensum=tmp
      ENDELSE 
      lensumtmp=0
  ENDIF 

END 

PRO plot_m2l_es, radbin, stripes=stripes, doallow=doallow, check=check, $
                 readfromfile=readfromfile, neighborlum=neighborlum,$
                 nocolor=nocolor, redolum=redolum

  IF n_params() LT 1 THEN BEGIN 
      print,'Syntax - plot_m2l_es, radbin, stripes=stripes,doallow=doallow, check=check, readfromfile=readfromfile, nocolor=nocolor, redolum=redolum'
      return
  END
                                ;Find out device

  fend = 'N1.fit'
  
  errfac = 1.32                 ;Multiply "combined" error bars by this factor

  IF n_elements(stripes) EQ 0 THEN stripes=[10,36,37,42,43,82]
  nstripe = n_elements(stripes)
  stripestr=''
  FOR i=0L, nstripe-1 DO stripestr = stripestr+'stripe'+ntostr(stripes[i])+'_'
  print,'stripestr = ',stripestr

  types=['lum1_','lum2_','lum3_','lum4_']
  names=['lum1 (low)','lum2','lum3','lum4 (high)']
  colors=['u','g','r','i','z']
  nclr=n_elements(colors)

  basedir = '/sdss5/data0/lensout/'
  ;; combined parameters are output to single directory: first stripe
  paramdir = basedir+'stripe'+ntostr(stripes[0])+'/sublum/' + colors +'/'
  ;; lensums are in each individual stripe dir, must be combined
  stdir = basedir + 'stripe'+ntostr(stripes)+'/'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; radial bin max(R)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  rmax_act = [100,180,260,340,420,500,$
              580,660,740,820,900,980]
  rmax_actstr = ['100','180','260','340','420','500',$
                 '580','660','740','820','900','980']
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; output fits file names
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(neighborlum) THEN addstr='_neigh' ELSE addstr=''
  lumbinfile = basedir+'mass2light/lumbin_rad'+rmax_actstr[radbin]+addstr+'.fit'
  typebinfile = basedir+'mass2light/typebin_rad'+rmax_actstr[radbin]+addstr+'.fit'

  IF keyword_set(nocolor) THEN ccstr='_bw' ELSE ccstr=''
  psfile = basedir+'mass2light/'+'plot_m2l_rad'+rmax_actstr[radbin]+addstr+$
    ccstr+'_N1.ps'
  IF keyword_set(check) THEN BEGIN 
      WHILE fexist(psfile) DO psfile=newname(psfile)
  ENDIF 
  begplot,name=psfile,/color
  !p.thick=7


  IF (!d.flags AND 1) EQ 0 THEN doX=1 ELSE doX=0

  nlum=n_elements(types)

  ;; binned by lum differently in 
  ;; each bandpass

  lum=fltarr(nclr,nlum)
  lumerr=fltarr(nclr,nlum)
  lumhigh=lumerr
  lumlow=lumerr
  m2l=fltarr(nclr,nlum)
  m2lerr = m2l
  mass=fltarr(nclr,nlum)
  masserr=mass
  sissigma=mass
  sissigmaerr=mass

  neigh = mrdfits(basedir+'mass2light/sublum_lumdens_fitpar.fit',1)
  IF keyword_set(neighborlum) THEN BEGIN 
      nxx = 1000.
      xxn=arrscl( findgen(nxx), 50./1000., 1000./1000. ) ;Mpc
  ENDIF 

  ;; read in allowed parameters for mass vs each luminosity
  IF keyword_set(doallow) THEN BEGIN 
      allowdir = '/sdss4/data1/esheldon/GAL_GAL/spectra/'
      uallow = mrdfits(allowdir + 'massvs_ulum_allowedpar.fit',1)
      gallow = mrdfits(allowdir + 'massvs_glum_allowedpar.fit',1)
      rallow = mrdfits(allowdir + 'massvs_rlum_allowedpar.fit',1)
      iallow = mrdfits(allowdir + 'massvs_ilum_allowedpar.fit',1)
      zallow = mrdfits(allowdir + 'massvs_zlum_allowedpar.fit',1)
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; loop over luminosity bins
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(readfromfile) THEN BEGIN 

      lumbin = mrdfits(lumbinfile, 1)

      mass = lumbin.mass
      masserr = lumbin.masserr
      lum = lumbin.lum
      lumerr = lumbin.lumerr
      m2l = lumbin.m2l
      m2lerr = lumbin.m2lerr
      sissigma = lumbin.sissigma
      sissigmaerr = lumbin.sissigmaerr

      
  ENDIF ELSE BEGIN 
      ngals=0L
      FOR i=0L,nlum-1 DO BEGIN
          print

          CASE i OF
              0:BEGIN & lumalpha = neigh.lum1alpha & lumnorm = neigh.lum1norm & END
              1:BEGIN & lumalpha = neigh.lum2alpha & lumnorm = neigh.lum2norm & END
              2:BEGIN & lumalpha = neigh.lum3alpha & lumnorm = neigh.lum3norm & END
              3:BEGIN & lumalpha = neigh.lum4alpha & lumnorm = neigh.lum4norm & END
              ELSE:
          ENDCASE 

          ;; loop over colors in this luminosity bin
          FOR clr=0L, nclr-1 DO BEGIN 
              
              file=paramdir[clr]+types[i]+'zgal_gal_'+stripestr+'fitparam_'+fend
              l=0
              l=mrdfits(file,1,/silent)
              mass[clr,i]=l.sismass[radbin]
              sissigma[clr,i] = sqrt( l.sissigma2[radbin] > 0. )
              ;; Multiply by sqrt(3), because each bandpass is highly correlated
              masserr[clr,i]=l.sismasserr[radbin]*sqrt(errfac)
              sissigmaerr[clr,i] = sqrt( (l.sissigmaerr2[radbin] > 0.)*sqrt(errfac) )
              
              
              FOR ist=0L, nstripe-1 DO BEGIN 
                  file=stdir[ist]+'sublum/'+colors[clr]+'/'+$
                    types[i]+'zgal_gal_stripe'+ntostr(stripes[ist])+'_r_lensum_'+fend
                  
                  plot_m2l_es_addstruct, file, lensum
              ENDFOR 
              IF clr EQ 2 THEN ngals=ngals+n_elements(lensum)

              IF keyword_set(redolum) THEN BEGIN 
                  wtheta_absmag_diffz, lensum.z1d, clr, lensum.petrocounts[clr],$
                    lensum.counts_model[1]-lensum.counts_model[2], $
                    absmag, tmplum
                  lensum.lum[clr] = tmplum
              ENDIF 

              tave=0
              terr=0
              lensave,lensum,'lum',tave,terr,element=clr,pairelement=radbin

              IF keyword_set(neighborlum) THEN BEGIN 
                  ww=where(xxn LE rmax_act[radbin]/1000., nww)
                  la=lumalpha[clr]
                  lumdens=lumnorm[clr]*xxn^(-la)
                  neighlumin = lumdens*2./(2.-la)*!pi*xxn^2
                  print,'Cent = ',tave/1.e10,' Neigh = ',neighlumin[nww-1]
                  lumtot = tave + neighlumin*1.e10 ;lum
                  tave = lumtot[nww-1]
              ENDIF 

              lum[clr, i] = tave
              lumerr[clr,i] = terr
              lumlow[clr,i] = min(lensum.lum[clr])
              lumhigh[clr,i] = max(lensum.lum[clr])
                                ;print,'meanlum = ',lum[clr,i],' lumlow = ',lumlow[clr,i],' lumhigh=',lumhigh[clr,i]
              
              m2l[clr,i] = mass[clr,i]/lum[clr,i]
              m2lerr[clr,i] = m2l[clr,i]*sqrt( (masserr[clr,i]/mass[clr,i])^2 + $
                                               (lumerr[clr,i]/lum[clr,i] )^2 )
              ;; delete variable
              delvarx,lensum
          ENDFOR
          
      ENDFOR
      
      mass = mass/1.e12
      masserr = masserr/1.e12
      lum = lum/1.e10
      lumerr = lumerr/1.e10
      lumlow=lumlow/1.e10
      lumhigh=lumhigh/1.e10
  ENDELSE 
  
  print
  ;print,'Ngals = ',ngals

;  pcolor=[!black,!blue,!green,!magenta,!red]
  IF keyword_set(nocolor) THEN BEGIN 
      pcolor = replicate(!black, 5)
      typecolor=replicate(!black,3)
      lines = [0,1,2,3,4]
  ENDIF ELSE BEGIN 
      pcolor=[!blue,!lightblue,!green,!magenta,!red]
      typecolor=[!blue,!black,!red]
      lines = replicate(0,5)
  ENDELSE 
  
;  xtitle='Petrosian Luminosity  (h!U'+!tsym.minus+'2!N 10!U10!N L!DSun!N)'
;  ytitle = 'M(< '+rmax_actstr[radbin]+' kpc) / L (h M!DSun!N/L!DSun!N)'
  xtitle='Petrosian Luminosity  (h!U'+!tsym.minus+'2!N 10!U10!N L'+sunsymbol()+')'
  ytitle = 'M(< '+rmax_actstr[radbin]+' kpc) / L (h M'+sunsymbol()+' /L'+sunsymbol()+')'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; m/l vs l, each overplotted
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  CASE radbin OF
      5: m2lyrange=[0,1000]
      2: m2lyrange=[0,800]
      1: m2lyrange=[0,600]
      0: m2lyrange=[0,300]
      ELSE: m2lyrange=[0,1200]
  ENDCASE 

  erase
  o=sort(lum[3,*])
  lum=lum[*,o]
  yrange=m2lyrange
  xrange=[0,max(lum,wmax)]
  FOR clr=nclr-1,0L,-1 DO BEGIN 
      IF clr EQ nclr-1 THEN BEGIN 
          aploterror,1,lum[clr,*],m2l[clr,*],m2lerr[clr,*],$
            xrange=xrange,yrange=yrange,charsize=1.5,$
            xtitle=xtitle,ytitle=ytitle,ystyle=1,subtitle=subtitle,title=title,$
            line=lines[clr]
      ENDIF 
      
      oploterror,lum[clr,*],m2l[clr,*],m2lerr[clr,*],color=pcolor[clr], $
        errcolor=pcolor[clr],line=lines[clr]
  ENDFOR 
  legend,[!colorsp],/top,/right,color=pcolor,linestyle=lines,charsize=1.5,thick=replicate(!p.thick,5)


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; also do each separate plot
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF doX THEN key=get_kbrd(1)
  
  !p.multi = [0,2,3]
  
  FOR clr=0L,nclr-1 DO BEGIN 
      aploterror,1,lum[clr,*],m2l[clr,*],m2lerr[clr,*], $
        charsize=2,yrange=yrange,$
        xtitle=xtitle,ytitle=ytitle,ystyle=1,$
        subtitle=subtitle,title=title,ticklen=0.03
      oploterror,lum[clr,*],m2l[clr,*],m2lerr[clr,*],$
        color=pcolor[clr],errcolor=pcolor[clr]
      legend, !colorsp[clr],/right,box=0
  ENDFOR 
  !p.multi=0
  

  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;; m vs l, overplotted
  ;;;;;;;;;;;;;;;;;;;;;;;;

  IF doX THEN key=get_kbrd(1)

  xrange=[min(lum),max(lum)]
  yrange=prange(mass,masserr)
;  mytitle='M(< '+rmax_actstr[radbin]+' kpc) (10!U12!N M!DSun!N)'
  mytitle='M(< '+rmax_actstr[radbin]+' kpc) (10!U12!N M'+sunsymbol()+')'
  aplot,1,xrange,yrange,/nodata,xtitle=xtitle,$
    ytitle=mytitle, charsize=1.5,subtitle=subtitle,title=title,line=lines[nclr-1]
  slopes=fltarr(5)
  FOR clr=nclr-1,0L,-1 DO BEGIN 
      oploterror,lum[clr,*],mass[clr,*],lumerr[clr,*],masserr[clr,*],$
        color=pcolor[clr],errcolor=pcolor[clr], line=lines[clr]
;      oplotderror, lum[clr,*], mass[clr,*], lumlow[clr,*], lumhigh[clr,*], $
;        mass[clr,*]-masserr[clr,*], mass[clr,*]+masserr[clr,*], $
;        color=pcolor[clr],errcolor=pcolor[clr] ;,psym=1

  ENDFOR
  
  legend,!colorsp,$
    color=pcolor,charsize=1.5,/top,/left,thick=replicate(!p.thick,5),$
    line=lines
;      legend,['Best Fit slopes',colors+': '+ntostr(slopes)],/bottom,/right
  
  IF doX THEN key=get_kbrd(1)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; now each bandpass separately
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF radbin NE 2 THEN !p.multi=[0,0,2] ELSE chsz=1.5
;  !p.multi=[0,2,0]

  npow = 400L
  nnorm = 400L
  bestnorm = fltarr(5)
  bestpow = fltarr(5)
  nxx = 1000
  xx=findgen(nxx)
  modely_low = xx
  modely_high = xx

  CASE radbin OF
      0: tnormmax=200
      1: tnormmax=300
      2: tnormmax=800
      3: tnormmax=800
      4: tnormmax=800
      5: tnormmax=1000
      9: tnormmax=1000
      ELSE: tnormmax=1000
  ENDCASE 

  FOR clr=0L,nclr-1 DO BEGIN 

      CASE clr OF
          0:BEGIN & normmax=tnormmax & normmin=0. & powmin=-.5 & powmax=1.0 & END 
          1:BEGIN & normmax=tnormmax & normmin=100. & powmin=0. & powmax=2.0 & END
          ELSE: BEGIN & normmax=tnormmax & normmin=0. & powmin=0. & powmax=2.0 & END
      ENDCASE 
      powvals = arrscl( findgen(npow), powmin, powmax )
      normvals = arrscl( findgen(nnorm), normmin, normmax )
;      addtit=' (h M!DSun!N/L!DSun!N)'
      addtit=' (h M'+sunsymbol()+' /L'+sunsymbol()+')'
      pow_chisq_conf,lum[clr,*],mass[clr,*]*100.,masserr[clr,*]*100.,$
        powvals, normvals,chisq_surf, $
        bestp,bestn,powlow,powhigh,normlow,normhigh, $
        xtitle=!tsym.beta,ytitle=!tsym.upsilon_cap+addtit,$
        /center,charsize=chsz, aspect=1,xrange=[0,2]
      bestpow[clr] = bestp
      bestnorm[clr] = bestn
      range2error, powlow,bestp,powhigh,perrlow,perrhigh
      range2error, normlow, bestn, normhigh, nerrlow, nerrhigh

      modelx = arrscl( xx, min(lum[clr,*]), max(lum[clr,*]) )
      IF radbin NE 2 THEN BEGIN 
          message =  ['norm: ',' ','pow: ']
          message[0]=message[0]+ntostr(long(rnd(bestn,-1)))+'!S!U'+$
            '+'+ntostr(long(rnd(nerrhigh[0])))+$
            '!R!D'+!tsym.minus+ntostr(long(rnd(nerrlow[0])))
          message[2]=message[2]+ntostr(rnd(bestp,2),4)+'!S!U'+$
            '+'+ntostr(rnd(perrhigh[0],2),4)+$
            '!R!D'+!tsym.minus+ntostr(rnd(perrlow[0],2),4)
          legend,message,/right,box=0

          aploterror,1,lum[clr,*],mass[clr,*],masserr[clr,*],xtitle=xtitle,$
            ytitle=ytitle, charsize=1.5,subtitle=subtitle,title=title,psym=1,/center
          legend, !colorsp[clr],/left,box=0,charsize=1.5
          oploterror, lum[clr,*],mass[clr,*],lumerr[clr,*],masserr[clr,*],psym=3
;      oplotderror, lum[clr,*],mass[clr,*],lumlow[clr,*],lumhigh[clr,*],$
;        mass[clr,*]-masserr[clr,*], mass[clr,*]+masserr[clr,*], psym=3
          oplot,modelx,(bestnorm[clr]/100.)*(modelx)^bestpow[clr],$
            color=pcolor[clr]
      ENDIF ELSE BEGIN 

          print,'norm: ',ntostr(long(rnd(bestn,-1))),$
            ' + '+ntostr(long(rnd(nerrhigh[0])))+$
            ' - '+ntostr(long(rnd(nerrlow[0])))
          print,'pow: ',ntostr(rnd(bestp,2),4)+$
            ' + '+ntostr(rnd(perrhigh[0],2),4)+$
            ' - '+ntostr(rnd(perrlow[0],2),4)

          IF clr EQ 1 THEN xrg=[0.5, 4.0] ELSE delvarx,xrg

          
;          IF clr EQ 0 THEN BEGIN 
;              pos=[ [.18, .54], [0.48, 0.775]]
;              legend, !colorsp[clr],/right,box=0,charsize=2.0
;          ENDIF ELSE BEGIN 
              pos=[ [.55, .54], [0.85, 0.775]]
              legend, !colorsp[clr],/left,box=0,charsize=2.0
;          ENDELSE 
          ploterror,lum[clr,*],mass[clr,*],masserr[clr,*],$
            subtitle=subtitle,title=title,psym=1,$
            /noerase, pos=pos, charsize=1,xrange=xrg,$
$;            xtitle='L (10!U10!N L!DSun!N)',ytitle='M!D260!N (10!U12!N M!DSun!N)'
            xtitle=xtitle,ytitle=ytitle
          oplot,modelx,(bestnorm[clr]/100.)*(modelx)^bestpow[clr];,$
;            color=pcolor[clr]
      ENDELSE 

      
      IF keyword_set(doallow) THEN BEGIN 
          command = 'allow = '+colors[clr]+'allow'
          print,command
          IF NOT execute(command) THEN message,'Dohh!'
          FOR ix=0L, nxx-1 DO BEGIN 
              mody = allow.normallow*(modelx[ix])^allow.powallow
              modely_low[ix] = min(mody)
              modely_high[ix] = max(mody)
          ENDFOR 
          oplot,modelx,modely_low,line=2
          oplot,modelx,modely_high,line=2
      ENDIF 
      IF doX THEN key=get_kbrd(1)
  ENDFOR 
  !p.multi=0
  

;  colprint,lum[2,*],lumerr[2,*]
;  print
;  print,'Signal to noise u g r i z'
;  colprint,mass[0,*]/masserr[0,*], mass[1,*]/masserr[1,*], mass[2,*]/masserr[2,*], mass[3,*]/masserr[3,*], mass[4,*]/masserr[4,*]
 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; copy into outputs structure
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  lumbin = create_struct('lum', lum, $
                         'lumerr', lumerr, $
                         'mass', mass, $
                         'masserr', masserr, $
                         'sissigma', sissigma, $
                         'sissigmaerr',sissigmaerr,$
                         'm2l', m2l, $
                         'm2lerr', m2lerr)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Now do by type
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  paramdir = basedir+'stripe'+ntostr(stripes[0])+'/'
  types = ['spiral_', 'main_', 'ellip_']
  ntypes = n_elements(types)

  lum = fltarr(nclr, ntypes)
  lumerr = lum
  m2l = lum
  m2lerr = lum
  mass = fltarr(ntypes)
  masserr = mass
  sissigma=mass
  sissigmaerr=mass

  maxlum = ntostr([10.0e10, 10.e10, 15.0e10, 30.0e10, 45.0e10])

  IF keyword_set(readfromfile) THEN BEGIN 

      typebin = mrdfits(typebinfile, 1)
      mass = typebin.mass
      masserr = typebin.masserr
      lum = typebin.lum
      lumerr = typebin.lumerr
      m2l = typebin.m2l
      m2lerr = typebin.m2lerr
      sissigma = typebin.sissigma
      sissigmaerr = typebin.sissigmaerr
  ENDIF ELSE BEGIN 
      FOR i=0L, ntypes-1 DO BEGIN 

          CASE i OF
              0: BEGIN & lumalpha=neigh.spiralalpha & lumnorm=neigh.spiralnorm & END
              1: BEGIN & lumalpha=neigh.lumalpha & lumnorm=neigh.lumnorm & END
              2: BEGIN & lumalpha=neigh.ellipalpha & lumnorm=neigh.ellipnorm & END
              ELSE:
          ENDCASE 
          file=paramdir+types[i]+'zgal_gal_'+stripestr+'fitparam_'+fend
          
          l=0
          l=mrdfits(file, 1, /silent)
          mass[i] = l.sismass[radbin]
          masserr[i] = l.sismasserr[radbin]*sqrt(errfac)
          sissigma[i] = sqrt( l.sissigma2[radbin] > 0. )
          sissigmaerr[i] = sqrt( (l.sissigmaerr2[radbin] > 0.)*sqrt(errfac) )
          
          FOR ist=0L, nstripe-1 DO BEGIN 
              
              file = stdir[ist]+types[i]+'zgal_gal_stripe'+$
                ntostr(stripes[ist])+'_r_lensum_'+fend
              plot_m2l_es_addstruct, file, lensum
          ENDFOR 
          
          FOR clr=0L, nclr-1 DO BEGIN 

              IF keyword_set(redolum) THEN BEGIN 
                  wtheta_absmag_diffz, lensum.z1d, clr, lensum.petrocounts[clr],$
                    lensum.counts_model[1]-lensum.counts_model[2], $
                    absmag, tmplum
                  lensum.lum[clr] = tmplum
              ENDIF 
              w=where(lensum.lum[clr] LT maxlum[clr] AND $
                      lensum.lum[clr] GT 0.0 )

              lensave, lensum[w], 'lum', tave, terr, element=clr,$
                pairelement=radbin

              IF keyword_set(neighborlum) THEN BEGIN 
                  ww=where(xxn LE rmax_act[radbin]/1000., nww)
                  la=lumalpha[clr]
                  lumdens=lumnorm[clr]*xxn^(-la)
                  neighlumin = lumdens*2./(2.-la)*!pi*xxn^2
                  print,types[i]+' Cent = ',tave/1.e10,' Neigh = ',neighlumin[nww-1]
                  lumtot = tave + neighlumin*1.e10 ;lum
                  tave = lumtot[nww-1]
              ENDIF 

              lum[clr,i] = tave
              lumerr[clr,i] = terr
              m2l[clr,i] = mass[i]/lum[clr,i]
              m2lerr[clr,i] = sqrt( (masserr[i]/mass[i])^2 + $
                                    (lumerr[clr,i]/lum[clr,i])^2)*m2l[clr,i]
          ENDFOR 
          ;; delete variable
          delvarx, lensum
      ENDFOR 
      
      mass = mass/1.e12
      masserr = masserr/1.e12
      lum = lum/1.e10
      lumerr = lumerr/1.e10
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; m/l vs l, each type of galaxy overplotted
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;ytitle = 'M(< '+rmax_actstr[radbin]+' kpc) / L'
  symbols = [2,1,6]
  symsize=replicate(1.3,3)
  erase
  o=sort(lum[3,*])
  lum=lum[*,o]
  xrange=[min(lum,wmin),max(lum,wmax)]
  yrange=m2lyrange
  FOR clr=nclr-1,0L,-1 DO BEGIN 
      IF clr EQ nclr-1 THEN BEGIN 
          aploterror,1,lum[clr,*],m2l[clr,*],m2lerr[clr,*],$
            xrange=xrange,yrange=yrange,charsize=1.5,$
            xtitle=xtitle,ytitle=ytitle,ystyle=1,subtitle=subtitle,title=title,$
            line=lines[clr]
      ENDIF 
      oploterror,lum[clr,*],m2l[clr,*],m2lerr[clr,*],color=pcolor[clr], $
        errcolor=pcolor[clr],line=lines[clr]

      myusersym,'fill_triangle'
      oplot,[lum[clr,0]],[m2l[clr,0]],psym=8,color=typecolor[0],symsize=1.3
      oplot,[lum[clr,1]],[m2l[clr,1]],psym=symbols[1],color=typecolor[1],symsize=1.3
      myusersym,'fill_circle'
      oplot,[lum[clr,2]],[m2l[clr,2]],psym=8,color=typecolor[2],symsize=1.3

  ENDFOR 
  lthick = replicate(!p.thick,3)
  legend,[!colorsp],/top,/right,color=pcolor,charsize=1.5,thick=replicate(!p.thick,5),line=lines
  myusersym,'fill_triangle'

  legend,['spiral','all galaxies','ellipticals'],symsize=symsize,$
    psym=[8,3,3],colors=typecolor, thick=lthick, /top, pos=[1.5,775.]
  legend,['spiral','all galaxies','ellipticals'],symsize=symsize,$
    psym=[3,symbols[1],3],colors=typecolor,thick=lthick,/top, $
    pos=[1.5,775.]
  myusersym,'fill_circle'
  legend,['spiral','all galaxies','ellipticals'],symsize=symsize,$
    psym=[3,3,8],colors=typecolor, thick=lthick, /top, pos=[1.5,775.]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Now each individually
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF doX THEN key=get_kbrd(1)
  
  !p.multi = [0,2,3]
  
  FOR clr=0L,nclr-1 DO BEGIN 
      aploterror,1,lumbin.lum[clr,*],lumbin.m2l[clr,*],lumbin.m2lerr[clr,*], $
        charsize=2,yrange=yrange,$
        xtitle=xtitle,ytitle=ytitle,ystyle=1,$
        subtitle=subtitle,title=title,ticklen=0.03
      oploterror,lumbin.lum[clr,*],lumbin.m2l[clr,*],lumbin.m2lerr[clr,*],$
        color=pcolor[clr],errcolor=pcolor[clr]
      myusersym,'fill_triangle'
      oploterror,[lum[clr,0]],[m2l[clr,0]],[m2lerr[clr,0]],$
        color=typecolor[0],psym=8,symsize=1.3
      myusersym,'fill_circle'
      oploterror,[lum[clr,2]],[m2l[clr,2]],[m2lerr[clr,2]],$
        color=typecolor[2],psym=8,symsize=1.3

      IF clr EQ 1 THEN BEGIN 
          myusersym,'fill_triangle'
          legend,['spiral','ellipticals'],symsize=symsize,thick=[!p.thick,!p.thick],$
            psym=[8,3],colors=[typecolor[0],typecolor[2]],/top,right=right,left=left
          myusersym,'fill_circle'
          legend,['spiral','ellipticals'],symsize=symsize,thick=[!p.thick,!p.thick],$
            psym=[3,8],colors=[typecolor[0],typecolor[2]],/top,right=right,left=left
      ENDIF 
      ;;legend, !colorsp[clr],/right,box=0
  ENDFOR 
  !p.multi=0

  xrange=fltarr(2)
  FOR clr=0L,nclr-1 DO BEGIN 

      yrange=prange(m2l[clr,*],m2lerr[clr,*])
      yrange[0] = yrange[0]*0.9
      yrange[1] = yrange[1]*1.1
      xrange[0] = 0.9*min(lum[clr,*])
      xrange[1] = 1.1*max(lum[clr,*])
      xx=arrscl( findgen(100), min(lum[clr,*]),max(lum[clr,*]) )
      aploterror,1,lum[clr,*],m2l[clr,*],m2lerr[clr,*],$
        xrange=xrange,yrange=yrange,charsize=1.5,psym=3,$
        xtitle=xtitle,ytitle=ytitle,ystyle=1,subtitle=subtitle,title=title

      oplot,xx,bestnorm[clr]*xx^(bestpow[clr]-1),color=pcolor[clr]

      myusersym,'fill_triangle'
      oplot,[lum[clr,0]],[m2l[clr,0]],psym=8,color=typecolor[0],symsize=1.3
      oplot,[lum[clr,1]],[m2l[clr,1]],psym=symbols[1],color=typecolor[1],symsize=1.3
      myusersym,'fill_circle'
      oplot,[lum[clr,2]],[m2l[clr,2]],psym=8,color=typecolor[2],symsize=1.3

      myusersym,'fill_triangle'
      IF clr EQ 0 THEN BEGIN
          right=1
          left=0
      ENDIF ELSE BEGIN 
          right=0
          left=1
      ENDELSE 
      legend,['spiral','all galaxies','ellipticals'],symsize=symsize,thick=lthick,$
        psym=[8,3,3],colors=typecolor,/top,right=right,left=left
      legend,['spiral','all galaxies','ellipticals'],symsize=symsize,thick=lthick,$
        psym=[3,symbols[1],3],colors=typecolor,/top,right=right,left=left
      myusersym,'fill_circle'
      legend,['spiral','all galaxies','ellipticals'],symsize=symsize,thick=lthick,$
        psym=[3,3,8],colors=typecolor,/top,right=right,left=left

  ENDFOR 

  typebin = create_struct('lum', lum, $
                          'lumerr', lumerr, $
                          'mass', mass, $
                          'masserr', masserr, $
                          'sissigma', sissigma, $
                          'sissigmaerr',sissigmaerr,$
                          'm2l', m2l, $
                          'm2lerr', m2lerr)

  print
  print,'Sigma_v in each luminosity bin (split differently in each bp)'
  colprint,'  '+ntostr(long(rnd(lumbin.sissigma[0,*])))+'+/-'+ntostr(long(rnd(lumbin.sissigmaerr[0,*]))), $
    '  '+ntostr(long(rnd(lumbin.sissigma[1,*])))+'+/-'+ntostr(long(rnd(lumbin.sissigmaerr[1,*]))), $
    '  '+ntostr(long(rnd(lumbin.sissigma[2,*])))+'+/-'+ntostr(long(rnd(lumbin.sissigmaerr[2,*]))), $
    '  '+ntostr(long(rnd(lumbin.sissigma[3,*])))+'+/-'+ntostr(long(rnd(lumbin.sissigmaerr[3,*]))), $
    '  '+ntostr(long(rnd(lumbin.sissigma[4,*])))+'+/-'+ntostr(long(rnd(lumbin.sissigmaerr[4,*]))),$
    format='(5A10)'
    
  print
  print,'Sigma_v in each of the types (spiral,all,elliptical)'
  colprint,'  '+ntostr(long(rnd(typebin.sissigma)))+'+/-'+ntostr(long(rnd(typebin.sissigmaerr))),$
    format='(A10)'
  print
  print,'M/L for each bandpass'
  colprint,'  '+ntostr(long(rnd(lumbin.m2l[0,*])))+'+/-'+ntostr(long(rnd(lumbin.m2lerr[0,*]))), $
    '  '+ntostr(long(rnd(lumbin.m2l[1,*])))+'+/-'+ntostr(long(rnd(lumbin.m2lerr[1,*]))), $
    '  '+ntostr(long(rnd(lumbin.m2l[2,*])))+'+/-'+ntostr(long(rnd(lumbin.m2lerr[2,*]))), $
    '  '+ntostr(long(rnd(lumbin.m2l[3,*])))+'+/-'+ntostr(long(rnd(lumbin.m2lerr[3,*]))), $
    '  '+ntostr(long(rnd(lumbin.m2l[4,*])))+'+/-'+ntostr(long(rnd(lumbin.m2lerr[4,*]))),$
    format='(5A10)'
;  colprint,lumbin.m2l[0,*],lumbin.m2l[1,*],lumbin.m2l[2,*],$
;    lumbin.m2l[3,*],lumbin.m2l[4,*]
  print

  print,'M/L by type (down) and bandpass (across)'
  colprint,'  '+ntostr(long(rnd(typebin.m2l[0,*])))+'+/-'+ntostr(long(rnd(typebin.m2lerr[0,*]))), $
    '  '+ntostr(long(rnd(typebin.m2l[1,*])))+'+/-'+ntostr(long(rnd(typebin.m2lerr[1,*]))), $
    '  '+ntostr(long(rnd(typebin.m2l[2,*])))+'+/-'+ntostr(long(rnd(typebin.m2lerr[2,*]))), $
    '  '+ntostr(long(rnd(typebin.m2l[3,*])))+'+/-'+ntostr(long(rnd(typebin.m2lerr[3,*]))), $
    '  '+ntostr(long(rnd(typebin.m2l[4,*])))+'+/-'+ntostr(long(rnd(typebin.m2lerr[4,*]))),$
    format='(5A10)'
;  colprint,typebin.m2l[0,*],typebin.m2l[1,*],typebin.m2l[2,*],$
;    typebin.m2l[3,*],typebin.m2l[4,*]
  print

  IF NOT keyword_set(readfromfile) THEN BEGIN 
      print,'Outputting: ',lumbinfile
      mwrfits, lumbin, lumbinfile, /create
      print,'Outputting: ',typebinfile
      mwrfits, typebin, typebinfile, /create
  ENDIF 
  endplot

  return
END
