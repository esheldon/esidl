PRO runomega_printheader, m2l_lun

printf, m2l_lun, '\begin{deluxetable}{cccccccc}'
printf, m2l_lun, '\tabletypesize{\small}'
printf, m2l_lun, '\tablecaption{M/L and $\Omega_M$ fits \label{m2lomegatable} }'
printf, m2l_lun, '\tablewidth{0pt}'
printf, m2l_lun, '\tablehead{'
printf, m2l_lun, '\colhead{Sample \tablenotemark{a}} &'
printf, m2l_lun, '\colhead{Bandpass \tablenotemark{b}} &'
printf, m2l_lun, '\colhead{$\alpha$ \tablenotemark{c}} &'
printf, m2l_lun, '\colhead{$L_0$ \tablenotemark{d}} &'
printf, m2l_lun, '\colhead{$R_s$ \tablenotemark{e}} &'
printf, m2l_lun, '\colhead{$M_0$ \tablenotemark{f}} &'
printf, m2l_lun, '\colhead{$M_0/L_0$ \tablenotemark{g}} &'
printf, m2l_lun, '\colhead{$\Omega_M$ \tablenotemark{h} } \\'
printf, m2l_lun, '&'
printf, m2l_lun, '&'
printf, m2l_lun, '&'
printf, m2l_lun, '\colhead{$10^{10} h^{-2} $L$_{\sun}$} &'
printf, m2l_lun, '\colhead{$h^{-1}$ kpc} &'
printf, m2l_lun, '\colhead{$10^{12} h^{-1} $M$_{\sun}$} &'
printf, m2l_lun, '\colhead{$h$ M$_{\sun}$/L$_{\sun}$} &'
printf, m2l_lun, '}'
printf, m2l_lun, '\startdata'


END 


PRO runomega, lumextno, mextno

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: runomega, lumextno, mextno'
      return
  ENDIF 

  meanlum = [0.784103, 0.889751, 1.51780, 2.07778, 2.58957]

  types = ['all','spiral','ellip','lum1','lum2']
  ntype = n_elements(types)
  clrs = [2,3,4]
  nclr=n_elements(clrs)

  outdir = '/sdss5/data0/lensout/mass2light/'

  m2l = fltarr(nclr, ntype) & m2lel=m2l & m2leh=m2l
  M0 = fltarr(nclr, ntype) & M0el=M0 & M0eh=M0
  L0 = fltarr(nclr, ntype) & L0el=L0 & L0eh=L0
  pow = fltarr(nclr, ntype) & powel=pow & poweh=pow
  omega = m2l & omegael=m2l & omegaeh=m2l

  m2ltex = outdir+'m2ltable.tex'
  omegatex = outdir+'omegatable.tex'

  print,'Tex file: ',m2ltex
  openw, m2l_lun, m2ltex, /get_lun
  
  runomega_printheader, m2l_lun

;  openw, omega_lun, omegatex, /get_lun

  FOR i=0L, nclr-1 DO BEGIN 
      clr=clrs[i]

      calc_m2l_omega_lumdens, clr, $
        tm2l, tm2lel, tm2leh,$
        tomega, tomegael, tomegaeh, $
        tM0, tM0el, tM0eh, tL0, tL0el, tL0eh,$
        tpow,tpowel,tpoweh,$
        /fixlum,/doplot,lumextno=lumextno,mextno=mextno

      m2l[i,*] = tm2l & m2lel[i,*]=tm2lel & m2leh[i,*]=tm2leh
      M0[i,*] = tM0 & M0el[i,*]=tM0el & M0eh[i,*]=tM0eh
      L0[i,*] = tL0 & L0el[i,*]=tL0el & L0eh[i,*]=tL0eh
      pow[i,*] = tpow & powel[i,*]=tpowel & poweh[i,*]=tpoweh
      ;;print,'tpowel = ',tpowel

      omega[i,*] = tomega & omegael[i,*]=tomegael & omegaeh[i,*]=tomegaeh

  ENDFOR 

  ;; convert M0 to 10^13 M_{\sun}
;  M0 = M0/10.0 & M0el = M0el/10.0 & M0eh = M0eh/10.0
;  L0 = L0/10.0 & L0el = L0el/10.0 & L0eh = L0eh/10.0

  ;; output to m2l,omega files
  FOR t=0L, ntype-1 DO BEGIN 
      FOR i=0L, nclr-1 DO BEGIN 

          clr=clrs[i]

          ;; calculate characteristic ratius where
          ;; mass to light reaches fraction f of the
          ;; value at infinity

          ml = meanlum[clr]
          f = 0.5
          Rs = ( f/(1.-f) * ml/L0[i,t] )^(1./(2.-pow[i,t]))*1000.
          Rserr1 = Rs*sqrt( 1./(2.-pow[i,t])^2*(L0el[i,t]/L0[i,t])^2 + $
                            (alog(L0[i,t]))^2/(2.-pow[i,t])^4 * powel[i,t]^2 )
          Rserr2 = Rs*sqrt( 1./(2.-pow[i,t])^2*(L0eh[i,t]/L0[i,t])^2 + $
                            (alog(L0[i,t]))^2/(2.-pow[i,t])^4 * poweh[i,t]^2 )
          Rserr3 = Rs*sqrt( 1./(2.-pow[i,t])^2*(L0el[i,t]/L0[i,t])^2 + $
                            (alog(L0[i,t]))^2/(2.-pow[i,t])^4 * poweh[i,t]^2 )
          Rserr4 = Rs*sqrt( 1./(2.-pow[i,t])^2*(L0eh[i,t]/L0[i,t])^2 + $
                            (alog(L0[i,t]))^2/(2.-pow[i,t])^4 * powel[i,t]^2 )
          Rserr = max([Rserr1,Rserr2,Rserr3,Rserr4])
          IF t EQ 0 THEN BEGIN
              print,'Rs = '+ntostr(Rs)+' +/- '+ ntostr(Rserr)
          ENDIF 

          IF i EQ 0 THEN tstr = types[t] ELSE tstr='   '
          m2lmes=tstr+' & '+!colors[clr]+'$^{\prime}$ & '+$
                 ntostr(rnd(pow[i,t],2), 4)+$
                 '$\pm '+ntostr(rnd(poweh[i,t],2), 4)+'$' + $
                 ' & '+$
                 ntostr(rnd(L0[i,t],1), 4)+$
                 '$^{+'+ntostr(rnd(L0eh[i,t],1), 3)+'}' + $
                  '_{-'+ntostr(rnd(L0el[i,t],1), 3)+'}$' + $
                 ' & '+$
                 ntostr(long(rnd(Rs)))+$
                 '$\pm '+ntostr(long(rnd(Rserr)))+'$' + $
                 ' & '+$
                 ntostr(rnd(M0[i,t],1), 4)+$
                 '$^{+'+ntostr(rnd(M0eh[i,t],1), 3)+'}' + $
                  '_{-'+ntostr(rnd(M0el[i,t],1), 3)+'}$' + $
                 ' & '+$
                 ntostr(long(rnd(m2l[i,t])), 5)+$
                 '$^{+'+ntostr(long(rnd(m2leh[i,t])), 5)+'}' + $
                  '_{-'+ntostr(long(rnd(m2lel[i,t])), 5)+'}$' + $
                 ' & '+$
                 ntostr(rnd(omega[i,t],2),4)+$
                 '$ \pm '+ntostr(rnd(omegaeh[i,t],2), 4)+'$'
          IF NOT ((t EQ ntype-1) AND (i EQ nclr-1) ) THEN m2lmes=m2lmes+' \\'
;                 ntostr(rnd(omega[i,t],2),4)+$
;                 '$^{+'+ntostr(rnd(omegaeh[i,t],2), 4)+'}' + $
;                  '_{-'+ntostr(rnd(omegael[i,t],2), 4)+'}$' + ' \\'

          printf,m2l_lun, m2lmes



      ENDFOR 
  ENDFOR 

  printf, m2l_lun, '\enddata'
  free_lun, m2l_lun

  extstr = 'N'+ntostr(long(lumextno))+'.fit'

  indir = '/sdss5/data0/lensout/mass2light/'
  lall = mrdfits(indir+'lumdense_allplot_'+extstr,1,/silent)
  alpha = abs(lall.pow)
  print,'alpha = ',alpha
  xx = arrscl(findgen(1000), 0.0, 2.)
  xxf = xx^(2. - alpha[2])
  m2lfunc_r = m2l[0, 0]*( xxf/( meanlum[2]/L0[0,0] + xxf ) )
  xxf = xx^(2. - alpha[3])
  m2lfunc_i = m2l[1, 0]*( xxf/( meanlum[3]/L0[1,0] + xxf ) )
  xxf = xx^(2. - alpha[4])
  m2lfunc_z = m2l[2, 0]*( xxf/( meanlum[4]/L0[2,0] + xxf ) )

  outps = indir + 'm2lfunc.ps'
  begplot,name=outps

  xtitle='Aperture Radius (h!U'+!tsym.minus+'1!N kpc)'
  ytitle='M/L (h M'+sunsymbol()+' /L'+sunsymbol()+')'

  xx=xx*1000.
  yrange=[0,300]
  lines = [0, 2, 3]
  aplot,!gratio,xx,m2lfunc_r, lines=lines[0],$
    yrange=yrange,xtitle=xtitle,ytitle=ytitle
  oplot,xx,m2lfunc_i, lines=lines[1]
  oplot,xx,m2lfunc_z, lines=lines[2]

  legend,!colorsp[2:4],lines=lines,thick=replicate(!p.thick,3)

  endplot


END 
