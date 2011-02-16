PRO generate_meanscinv, stripe, clr, npts, use_lambda=use_lambda, nops=nops, hirata=hirata, hudson=hudson, lrg_sources=lrg_sources, rlrg_sources=rlrg_sources

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: generate_meanscinv, stripe(s), clr(s), npts [, /use_lambda, /hirata, /nops, /hudson, /lrg_sources, /rlrg_sources]'
      print,'Npts=200 is good'
      return
  ENDIF 

  ;; generate array of mean_scritinv as a function of
  ;; lens redshift, mean source redshift and redshift
  ;; error
  ;;
  ;; npts=200 is good

  meanscinv_names, stripe, clr, fitfile, psfile, $
    use_lambda=use_lambda, hirata=hirata, hudson=hudson, $
    lrg_sources=lrg_sources, rlrg_sources=rlrg_sources
  print,fitfile
  print,psfile

  IF NOT keyword_set(nops) THEN begplot,name=psfile,/color
  setup_mystuff

  get_nz_photoz, stripe, clr, struct, hirata=hirata, hudson=hudson, $
    lrg_sources=lrg_sources, rlrg_sources=rlrg_sources
  prior_zs = struct.z
  prior_pofzs = struct.pofz

  ;; slightly broader than cuts
  ;; in our photoz catalogs

  nzl = 100
  nzs = 50
  nzsErr = 50

  zlMin = 0.01
  zlMax = 1.0
  zlStep = (zlMax-zlMin)/(nzl-1)

  zsMin = 0.01
  zsMax = 1.0
  zsStep = (zsMax-zsMin)/(nzs-1)

  zsErrMin = 0.01
  zsErrMax = 0.5
  zsErrStep = (zsErrMax-zsErrMin)/(nzsErr-1)

  zl = arrscl( findgen(nzl), zlmin, zlmax )
  zs = arrscl( findgen(nzs), zsmin, zsmax )
  zserr = arrscl( findgen(nzserr), zserrmin, zserrmax )
  
  IF n_elements(mean_scinv) EQ 0 THEN BEGIN 

      mean_scinv = dblarr(nzl, nzs, nzserr)
      
      tt=systime(1)
      FOR il = 0L, nzl-1 DO BEGIN 
          FOR is = 0L, nzs-1 DO BEGIN 
              FOR ise = 0L, nzserr-1 DO BEGIN 
                  mean_scinv[il, is, ise] = $
                    mean_sigmacritinv(zl[il], zs[is], zserr[ise],$
                                      npts, $
                                      prior_zs=prior_zs,  $
                                      prior_pofzs=prior_pofzs, $
                                      use_lambda=use_lambda)
              ENDFOR 
          ENDFOR 
      ENDFOR 
      ptime,systime(1)-tt
  ENDIF 

  ytit = !csym.sigma_cap+'!S!DCrit!N!U'+!csym.minus+'1!N  [10!U'+$
    !csym.minus+'4!N M'+sunsymbol()+'!U'+!csym.minus+'1!N pc!U2!N]'
  xtit = 'Z!DS!N Error'
  myusersym,'fill_circle'
  symsize=0.7

  IF !d.name EQ 'X' THEN BEGIN 
      R = long( arrscl( findgen(nzserr), 100L, 255L ) )
      G = 0L
      B = 0L
      
      colors = R + 256L*(G+256L*B)
  ENDIF ELSE BEGIN 
      colors = long( arrscl( lindgen(nzserr), !grey0, !grey90 ) )
  ENDELSE 
  fac = 1.e4

  FOR il=0L, nzl-1 DO BEGIN 

      minsc = min(mean_scinv[il,*,*]*fac)
      maxsc = max(mean_scinv[il,*,*]*fac)
      FOR is=0L, nzs-1 DO BEGIN 

          IF is EQ 0 THEN BEGIN  
              aplot, 1,zserr, mean_scinv[il, is, *]*fac, psym=8, $
                     symsize=symsize,$
                     xtit=xtit,ytit=ytit,$
                     yrange=[minsc,maxsc],/center
          ENDIF ELSE BEGIN 
              oplot, zserr, mean_scinv[il, is, *]*fac, psym=8, $
                     symsize=symsize,$
                     color=colors[is]
          ENDELSE 

      ENDFOR 
      legend,['Z!DL!N = '+ntostr(zl[il],6,/round),$
              'Npts = '+ntostr(npts)],/right, charsize=1
      IF !d.name EQ 'X' THEN key=get_kbrd(1)
  ENDFOR 

  IF NOT keyword_set(nops) THEN endplot

  struct = create_struct('npts',      long(npts), $
                         $
                         'zlStep',    float(zlStep), $
                         'zlMin',     float(zlMin),$
                         'zlMax',     float(zlMax), $
                         'zl',        float(zl), $
                         'zli',       findgen(nzl), $
                         $
                         'zsStep',    float(zsStep), $
                         'zsMin',     float(zsMin),$
                         'zsMax',     float(zsMax), $
                         'zs',        float(zs), $
                         'zsi',       findgen(nzs), $
                         $
                         'zsErrStep', float(zsErrStep), $
                         'zsErrMin',  float(zsErrMin),$
                         'zsErrMax',  float(zsErrMax), $
                         'zsErr',     float(zsErr), $
                         'zsErri',    findgen(nzserr), $
                         $      ; dblarr(nzl, nzs, nzserr)
                         'mean_scinv', double(mean_scinv)) 
  
  
  print,'output file: ',fitfile
  mwrfits, struct, fitfile, /create

END 
