pro esheldon_setup, silent=silent

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; misc. stuff that only need be defined once
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  defsysv, '!COLORS', exist=colors_exist
  IF NOT colors_exist THEN BEGIN 

      defsysv, '!COLORS', ['u','g','r','i','z']

      myusersym, 'fill_circle'

      ;; radians and degrees in double
      defsysv, '!d2r', !dpi/180.0d0 ;Change degrees to radians.
      defsysv, '!r2d', 180.0d0/!dpi ;radians to degrees

      ;; ASCII plusminus
      defsysv, '!plusminus', string(177b)
      ;defsysv, '!plusminus', '+/-' 

      ;; typical shape noise
;      defsysv, '!shapenoise', 0.32

      ;; boundary flags
;      defsysv, '!MINLAMFLAG', 2b^0
;      defsysv, '!MAXLAMFLAG', 2b^1
;      defsysv, '!MINETAFLAG', 2b^2
;      defsysv, '!MAXETAFLAG', 2b^3

;      defsysv, '!FLAGS_MASKED', '1'X
;      defsysv, '!FLAGS_QUAD1_MASKED', '2'X
;      defsysv, '!FLAGS_QUAD2_MASKED', '4'X
;      defsysv, '!FLAGS_QUAD3_MASKED', '8'X
;      defsysv, '!FLAGS_QUAD4_MASKED', '10'X

;      defsysv, '!FLAGS_QUAD1_MASKED_MONTE', '20'X
;      defsysv, '!FLAGS_QUAD2_MASKED_MONTE', '40'X
;      defsysv, '!FLAGS_QUAD3_MASKED_MONTE', '80'X
;      defsysv, '!FLAGS_QUAD4_MASKED_MONTE', '100'X

;      defsysv, '!FLAGS_MASKED_BASIC','200'X
;      defsysv, '!FLAGS_MASKED_SIMPLE','400'X
;      defsysv, '!FLAGS_MASKED_BOUND','1000'X
;      defsysv, '!FLAGS_MASKED_COMBINED','2000'X

      ;; convert rho_crit to M_{\sun} Mpc^{-3}
      rhounits = double(1.e18)
      rhocrit = double(2.77545e-7) ; Msolar/pc^3
      rhocrit = rhocrit*rhounits ; Msolar/Mpc^3

      defsysv, '!rhocrit', rhocrit

      ;; hard cut on smear polarizability and galaxy probability
;      defsysv, '!hardrcut', 0.8
;      defsysv, '!hardprobcut', 0.8

;      defsysv, '!PROB_PSF_DEFAULT', -9999.
;      defsysv, '!PROBFLAG_NOITER_REQUESTED', 2b^0
;      defsysv, '!PROBFLAG_NOITER', 2b^1
;      defsysv, '!PROBFLAG_NOPROB', 2b^2

      ;; sun absmag
      sunmag = [6.38d,5.06d,4.64d,4.53d,4.52d]

      ;; from blanton
      mstar = [-18.34d, -20.04d, -20.83d, -21.26d, -21.55d]
      mstar_err = [0.08d, 0.04d, 0.03d, 0.04d,0.04d]

      lstar = 10d^( (mstar - sunmag)/(-2.5d) )/1.e10

      defsysv, '!SUNMAG', sunmag
      defsysv, '!MSTAR', mstar
      defsysv, '!MSTAR_ERR', mstar_err
      defsysv, '!LSTAR', lstar

  ENDIF 



  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Depends on device: redefine each time
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;defsymbols

  if 0 then begin
      defsysv, '!COLORSP', !COLORS+!csym.prime
	  kpcxtitle='R [h^{-1} kpc]'
      kpcxtitle2=kpcxtitle

      mpcxtitle='R [h^{-1} Mpc]'
	  mpcxtitle2=mpcxtitle

      mpcxtitle3d = 'r [h^{-1} Mpc]'

	  shytitle = '\gamma_T'
	  orthoytitle = 'gamma_{\times}'

      ;sigytitle = '!S'+!csym.sigma_cap+'!R!A'+!csym.minus+'!N ('+!csym.ltequal+'R) '+$
      ;    !csym.minus+' !S'+!csym.sigma_cap+'!R!A'+!csym.minus+$
      ;    '!N (R) [h M'+sunsymbol()+' pc!U'+!csym.minus+'2!N]'
	  deltaytitle = '\Delta\Sigma [h M_{\sun} pc^{-2}'
	  rdeltaytitle = '\Delta\Sigma_{\times} [h M_{\sun} pc^{-2}'
	  randdeltaytitle = '\Delta\Sigma_{rand} [h M_{\sun} pc^{-2}'

	  lumytitle = 'L/Area [10^{10} L_{\sun} Mpc^{-2}'

      wgmytitle = 'w!Dgm!N  [h!U'+!csym.minus+'1!N Mpc]'

      xigmxtitle = 'r [h!U'+!csym.minus+'1!N Mpc]'
      xigmytitle = !csym.xi+'!Dgm!N(r) '+!csym.times+' '+!csym.omega_cap+'!Dm!N /0.27'

	  rhoytitle = '\rho(r) - !S\rho!R!U-!N} [h^{2} 10^{12} M_{\sun} Mpc^{-3}]'
	  rhoytitle2 = '\rho(r) - !S\rho!R!U-!N} [h^{2} 10^{12} M_{\sun} pc^{-3}]'

	  mytitle = 'M(<r) [10^{12} h^{-1} M_{\sun}]'
	  mytitle2 = 'M(<r) [10^{14} h^{-1} M_{\sun}]'

	  m200ytitle = 'M_{200} [10^{14} h^{-1} M_{\sun}]'
      r200ytitle = 'r_{200} [h^{-1} Mpc]'
      mVirYtitle = 'M_{Vir} [10^{14} h^{-1} M_{\sun}]'
      rVirYtitle = 'r_{Vir} [h^{-1} Mpc]'

      defsysv, '!kpcxtitle',       kpcxtitle
      defsysv, '!kpcxtitle2',      kpcxtitle2
      defsysv, '!mpcxtitle',       mpcxtitle
      defsysv, '!mpcxtitle2',      mpcxtitle2
      defsysv, '!mpcxtitle3d',      mpcxtitle3d

      ;defsysv, '!sigytitle',       sigytitle
      defsysv, '!deltaytitle',     deltaytitle
      defsysv, '!rdeltaytitle',    rdeltaytitle
      defsysv, '!randdeltaytitle', randdeltaytitle

      defsysv, '!shytitle',        shytitle
      defsysv, '!orthoytitle',     orthoytitle

      defsysv, '!lumytitle',       lumytitle

      defsysv, '!wgmytitle',       wgmytitle

      defsysv, '!xigmxtitle',      xigmxtitle
      defsysv, '!xigmytitle',      xigmytitle

      defsysv, '!rhoytitle',       rhoytitle
      defsysv, '!rhoytitle2',      rhoytitle2
      defsysv, '!mytitle',         mytitle
      defsysv, '!mytitle2',        mytitle2
      defsysv, '!m200ytitle',      m200ytitle
      defsysv, '!r200ytitle',      r200ytitle
      defsysv, '!mVirytitle',      mVirytitle
      defsysv, '!rVirytitle',      rVirytitle

  endif


  return


END 
