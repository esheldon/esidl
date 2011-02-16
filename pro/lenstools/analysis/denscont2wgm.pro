PRO denscont2wgm, struct, wgm, wgmerr, r0, r0errl, r0errh, r0low, r0high, r0_wgg=r0_wgg, pow_wgg=pow_wgg, doplot=doplot, rebin=rebin, level=level, silent=silent, plot_wggdata=plot_wggdata

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: denscont2wgm, struct, wgm, wgmerr, r0, r0errl, r0errh, r0low, r0high, r0_wgg=r0_wgg, pow_wgg=pow_wgg, doplot=doplot, rebin=rebin, silent=silent, plot_wggdata=plot_wggdata'
      print
      print,'For flux limited sample r0_wgg = 5.77, pow_wgg = 1.8'
      return
  ENDIF 

  IF keyword_set(rebin) THEN BEGIN 
      meanr = struct.meanr_rebin
      sigma = struct.sigma_rebin
      sigmaerr = struct.sigmaerr_rebin
  ENDIF ELSE BEGIN 
      meanr = struct.meanr
      sigma = struct.sigma
      sigmaerr = struct.sigmaerr
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; calculate w_{gm} in units of Mpc (same as Idit)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  ;; convert density contrast to M_{\sun} Mpc^{-2}
  dsigunits = double(1.e12)
  ;; convert rho_crit to M_{\sun} Mpc^{-3}
  rhounits = double(1.e18)

  ;; convert from density contrast to excess density 
  ;; based on the power law fit
  power = abs(struct.power)
  normfac = (2. - power)/power

  rho_crit = double(2.77545e-7) ; Msolar/pc^3
  rho_crit = rho_crit*rhounits    ; Msolar/Mpc^3

  rhobar = rho_crit*!omegam
  rhobarerr = 0.0;rho_crit*!omegamerr

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; calculate excess density and convert to wgm
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  sigma = sigma*dsigunits*normfac
  sigmaerr = sigmaerr*dsigunits*normfac
  norm = struct.norm*dsigunits*normfac

  wgm = sigma/rhobar
  wgmerr = sigmaerr/rhobar

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; calculate r0
  ;; remember, power is negative
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  powlow = abs(struct.powhigh[0])
  powhigh = abs(struct.powlow[0])

  normlow = struct.normlow[0]*dsigunits*normfac
  normhigh = struct.normhigh[0]*dsigunits*normfac

  calc_r0, norm, normlow,normhigh,$
           power, powlow,powhigh, $
           rhobar,rhobar-rhobarerr,rhobar+rhobarerr, $
           r0, r0errl, r0errh, r0low, r0high, level=level, silent=silent

  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; calculate Idit's w_gg
  ;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(r0_wgg) NE 0 AND n_elements(pow_wgg) NE 0 THEN BEGIN 
;      r0_wgg = 5.77              ;Mpc
;      pow_wgg = 1.8
      
      dowgg = 1
      f=gamma(0.5)*gamma(0.5*(pow_wgg-1))/gamma(0.5*pow_wgg)
      
      wggnorm = r0_wgg^pow_wgg*f
      
      wgg = wggnorm*(meanr/1000.)^(1.-pow_wgg)
      wgg_rad = meanr/1000.

  ENDIF ELSE IF keyword_set(plot_wggdata) THEN BEGIN 
      file = '~/wp_full_mjk_new.dat'
      readcol,file,wgg_rad,junk,wgg,junk,wggerr
      wtmp = where(wgg_rad LT 10.0)
      wgg_rad = wgg_rad[wtmp]
      wgg = wgg[wtmp]
      wggerr = wggerr[wtmp]
      dowgg = 1
  ENDIF ELSE dowgg=0

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; the plot
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(doplot) THEN BEGIN 
      yfit = struct.yfit*dsigunits*normfac/rhobar

      xrange = [0.01,20.]

      aplot, 1.0, meanr/1000., yfit,$
             /xlog,/ylog,xrange=xrange,xsty=1,$
             ytit=!wgmytitle,xtit=!mpcxtitle2,yrange=[5,4000],ysty=1,$
             xtickf='loglabels',ytickf='loglabels',xticklen=0.04, yticklen=0.04
      
      IF dowgg THEN BEGIN
          IF !d.name EQ 'X' THEN oclr=!green ELSE oclr=!p.color
          mess = ['Best Fit', $
                  'w!Dgg!N']
          colors = [!p.color, oclr]
          IF keyword_set(plot_wggdata) THEN BEGIN 
              line = 2
              ;;oplot, wgg_rad, wgg, color=oclr, psym=4
              oplot, wgg_rad, wgg, color=oclr, line=line
              
              ;; psym = [8,4]
              lines = [0,line]
              legend, mess, line=lines, colors=colors, $
                      /right, box=0,thick=[!p.thick,!p.thick]
          ENDIF ELSE BEGIN 
              line = 2
              oplot, wgg_rad, wgg, color=oclr, line=line
              
              lines = [0,line]
              legend, mess, line=lines, colors=colors, $
                      /right, box=0,thick=[!p.thick,!p.thick]
          ENDELSE 
      ENDIF 
      rr = meanr/1000.
      oploterror, rr, wgm, wgmerr, psym=8
      
      ;; wierd bug, must replot
      w=(lindgen(18))[12:17]
      oploterror, rr[w],wgm[w],wgmerr[w],psym=8

      IF dowgg AND NOT keyword_set(plot_wggdata) THEN BEGIN 
          IF !d.name EQ 'X' THEN key=get_kbrd(1)
          
          ratio = wgm/wgg
          ratioerr = ratio*(wgmerr/wgm)
          
          ytit = 'w!Dgm!N / w!Dgg!N'
          aploterror, !gratio, meanr/1000., ratio, ratioerr, $
                      psym=8,/xlog,xrange=xrange,/xsty, $
                      ytitle=ytit,xtitle=!mpcxtitle2
          oplot,[1.e-3,100.],[1,1], line=line
      ENDIF 


  ENDIF 


END 
