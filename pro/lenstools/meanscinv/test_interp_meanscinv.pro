PRO test_interp_meanscinv, stripe, clr, mean_scinv, mean_scinv2

  get_meanscinv, stripe, clr, scinv_struct, /use_lambda

  num = 10000
  npts = 200

  zlmin = 0.01
  zlmax = 1.0                   ;good enough for main (see below for lrg)

  get_nz_photoz, 10, [2,3], struct
  prior_zs = struct.z
  prior_pofzs = struct.pofz

  ;; slightly broader than cuts
  ;; in our photoz catalogs
  zsmin = 0.01
  zsmax = 1.0
  zserrmin = 0.01
  zserrmax = 0.5

  zl = arrscl( randomu(seed,num), zlmin, zlmax, arrmin=0.0, arrmax=1.0)
  zs = arrscl( randomu(seed,num), zsmin, zsmax, arrmin=0.0, arrmax=1.0)
  zserr = arrscl( randomu(seed,num), zserrmin, zserrmax, arrmin=0.0, arrmax=1.0)
      
  tt=systime(1)
  mean_scinv = interp_meanscinv(zl, zs, zserr, scinv_struct)
  print,systime(1)-tt

  mean_scinv2 = dblarr(num)

  tt=systime(1)
  FOR i = 0L, num-1 DO BEGIN 

      mean_scinv2[i] = $
        mean_sigmacritinv(zl[i], zs[i], zserr[i],$
                          npts, $
                          prior_zs=prior_zs,  $
                          prior_pofzs=prior_pofzs)

  ENDFOR 
  print,systime(1)-tt

return

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



END 
