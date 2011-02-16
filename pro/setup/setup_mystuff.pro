PRO setup_mystuff, silent=silent

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; run the sdssidl_setup program
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  sdssidl_setup

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; setup the defined structure
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  defsysv, '!MYIDL_DEF', exist=defexist

  IF NOT defexist THEN BEGIN 
      myidl_def=create_struct('config_file_read',0b,$
                              'mybindir_defined', 0b, $
                              'lensout_dir_defined',0b,$
                              'compea4_sofile_defined', 0b,$
                              'compea4_entry_defined',0b,$
                              'jack_sofile_defined',0b,$
                              'jack_entry_defined',0b,$
                              'gauleg_sofile_defined',0b,$
                              'gauleg_entry_defined',0b, $
                              'www_dir_defined', 0b)


      defsysv, '!MYIDL_DEF', myidl_def
      IF NOT keyword_set(silent) THEN BEGIN 
          print,'Structure !MYIDL_DEF defined.'
      ENDIF 
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Read config file
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF NOT !MYIDL_DEF.CONFIG_FILE_READ THEN BEGIN 

      config_file = getenv('MYIDL_CONFIG')
      IF config_file[0] EQ '' THEN BEGIN 
          config_file = $
            '~esheldon/.idl_config/myidl_setup_'+$
            sdssidl_config('hostname')+'.config'
          print, "Trying "+config_file
      ENDIF
      
      IF NOT fexist(config_file) THEN BEGIN
          message,'Cannot find config file: '+config_file,/inf
      ENDIF ELSE BEGIN 

          config_file=config_file[0]
          print,'Loading config file: ',config_file
          parse_config, config_file, options, values
          !myidl_def.config_file_read = 1b

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; define system dependent variables here if relevant. Leave 
          ;; undefined if they aren't; also, don't redefine once defined
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          IF NOT !MYIDL_DEF.MYBINDIR_DEFINED THEN BEGIN 
              w=where(strlowcase(options) EQ 'mybindir',nw)
              IF nw NE 0 THEN BEGIN 
                  mybindir = values[w[0]]
                  defsysv, '!MYBINDIR', mybindir
                  !myidl_def.mybindir_defined = 1b
              ENDIF 
          ENDIF 
          
          IF NOT !MYIDL_DEF.LENSOUT_DIR_DEFINED THEN BEGIN 
              w=where(strlowcase(options) EQ 'lensout_dir',nw)
              IF nw NE 0 THEN BEGIN 
                  lensout_dir = values[w[0]]
                  defsysv, '!LENSOUT_DIR', lensout_dir
                  !myidl_def.lensout_dir_defined = 1b
              ENDIF 
          ENDIF 
          
          IF NOT !MYIDL_DEF.COMPEA4_SOFILE_DEFINED THEN BEGIN 
              w=where(strlowcase(options) EQ 'compea4_sofile',nw)
              IF nw NE 0 THEN BEGIN 
                  compea4_sofile = values[w[0]]
                  defsysv, '!COMPEA4_SOFILE', compea4_sofile
                  !myidl_def.compea4_sofile_defined = 1b
              ENDIF 
          ENDIF 
          
          IF NOT !MYIDL_DEF.COMPEA4_ENTRY_DEFINED THEN BEGIN 
              w=where(strlowcase(options) EQ 'compea4_entry',nw)
              IF nw NE 0 THEN BEGIN 
                  compea4_entry = values[w[0]]
                  defsysv, '!COMPEA4_ENTRY', compea4_entry
                  !myidl_def.compea4_entry_defined = 1b
              ENDIF 
          ENDIF 
          
          IF NOT !MYIDL_DEF.JACK_SOFILE_DEFINED THEN BEGIN 
              w=where(strlowcase(options) EQ 'jack_sofile',nw)
              IF nw NE 0 THEN BEGIN 
                  jack_sofile = values[w[0]]
                  defsysv, '!JACK_SOFILE', jack_sofile
                  !myidl_def.jack_sofile_defined = 1b
              ENDIF 
          ENDIF 
          
          IF NOT !MYIDL_DEF.JACK_ENTRY_DEFINED THEN BEGIN 
              w=where(strlowcase(options) EQ 'jack_entry',nw)
              IF nw NE 0 THEN BEGIN 
                  jack_entry = values[w[0]]
                  defsysv, '!JACK_ENTRY', jack_entry
                  !myidl_def.jack_entry_defined = 1b
              ENDIF 
          ENDIF 
          
          IF NOT !MYIDL_DEF.GAULEG_SOFILE_DEFINED THEN BEGIN 
              w=where(strlowcase(options) EQ 'gauleg_sofile',nw)
              IF nw NE 0 THEN BEGIN 
                  gauleg_sofile = values[w[0]]
                  defsysv, '!GAULEG_SOFILE', gauleg_sofile
                  !myidl_def.gauleg_sofile_defined = 1b
              ENDIF 
          ENDIF 

          IF NOT !MYIDL_DEF.GAULEG_ENTRY_DEFINED THEN BEGIN 
              w=where(strlowcase(options) EQ 'gauleg_entry',nw)
              IF nw NE 0 THEN BEGIN 
                  gauleg_entry = values[w[0]]
                  defsysv, '!GAULEG_ENTRY', gauleg_entry
                  !myidl_def.gauleg_entry_defined = 1b
              ENDIF 
          ENDIF 
          
          IF NOT !MYIDL_DEF.WWW_DIR_DEFINED THEN BEGIN 
              w=where(strlowcase(options) EQ 'www_dir',nw)
              IF nw NE 0 THEN BEGIN 
                  www_dir = values[w[0]]
                  defsysv, '!WWW_DIR', www_dir
                  !myidl_def.www_dir_defined = 1b
              ENDIF 
          ENDIF 

      ENDELSE ;; found and read config file

  ENDIF ;; config already read?

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

      ;; text plusminus
      defsysv, '!plusminus', string(177b)

      ;; typical shape noise
      defsysv, '!shapenoise', 0.32

      ;; boundary flags
;      defsysv, '!MINLAMFLAG', 2b^0
;      defsysv, '!MAXLAMFLAG', 2b^1
;      defsysv, '!MINETAFLAG', 2b^2
;      defsysv, '!MAXETAFLAG', 2b^3

      defsysv, '!FLAGS_MASKED', '1'X
      defsysv, '!FLAGS_QUAD1_MASKED', '2'X
      defsysv, '!FLAGS_QUAD2_MASKED', '4'X
      defsysv, '!FLAGS_QUAD3_MASKED', '8'X
      defsysv, '!FLAGS_QUAD4_MASKED', '10'X

      defsysv, '!FLAGS_QUAD1_MASKED_MONTE', '20'X
      defsysv, '!FLAGS_QUAD2_MASKED_MONTE', '40'X
      defsysv, '!FLAGS_QUAD3_MASKED_MONTE', '80'X
      defsysv, '!FLAGS_QUAD4_MASKED_MONTE', '100'X

      defsysv, '!FLAGS_MASKED_BASIC','200'X
      defsysv, '!FLAGS_MASKED_SIMPLE','400'X
      defsysv, '!FLAGS_MASKED_BOUND','1000'X
      defsysv, '!FLAGS_MASKED_COMBINED','2000'X

      ;; don't use this any more, use absmag cuts
      ;; from read_cuts.pro
;      maxlum = [10.0, 10.0, 15.0, 20.0, 30.0]
;      defsysv, '!maxlum', maxlum

      ;; some physical parameters
      defsysv,'!h', 0.71
      defsysv,'!herr', 0.04
      
      defsysv,'!omegam', 0.27
      defsysv,'!omegamerr',0.02

      ;; convert rho_crit to M_{\sun} Mpc^{-3}
      rhounits = double(1.e18)
      rhocrit = double(2.77545e-7) ; Msolar/pc^3
      rhocrit = rhocrit*rhounits ; Msolar/Mpc^3

      defsysv, '!rhocrit', rhocrit

      ;; delta chi-squared levels
      defsysv,'!siglevels1', [1.00, 4.00, 9.00] ;levels for single parameters
      defsysv,'!siglevels2', [2.30, 6.17, 11.8] ;levels for joint probability

      ;; hard cut on smear polarizability and galaxy probability
      defsysv, '!hardrcut', 0.8
      defsysv, '!hardprobcut', 0.8

      defsysv, '!PROB_PSF_DEFAULT', -9999.
      defsysv, '!PROBFLAG_NOITER_REQUESTED', 2b^0
      defsysv, '!PROBFLAG_NOITER', 2b^1
      defsysv, '!PROBFLAG_NOPROB', 2b^2

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

  defsymbols
  defsysv, '!COLORSP', !COLORS+!csym.prime
  kpcxtitle='Projected Radius [h!U'+!csym.minus+'1!N kpc]'
  kpcxtitle2='R [h!U'+!csym.minus+'1!N kpc]'

  mpcxtitle='Projected Radius [h!U'+!csym.minus+'1!N Mpc]'
  mpcxtitle2='R [h!U'+!csym.minus+'1!N Mpc]'
  
  shytitle = !csym.gamma+'!DT!N'
  orthoytitle = !csym.gamma+'!D'+!csym.times+'!N /10!U'+!csym.minus+'3!N'

  sigytitle = '!S'+!csym.sigma_cap+'!R!A'+!csym.minus+'!N ('+!csym.ltequal+'R) '+$
    !csym.minus+' !S'+!csym.sigma_cap+'!R!A'+!csym.minus+$
    '!N (R) [h M'+sunsymbol()+' pc!U'+!csym.minus+'2!N]'
;  deltaytitle=!csym.delta_cap+!csym.sigma_cap+'!D'+!csym.plus+$
;    '!N [h M'+sunsymbol()+' pc!U'+!csym.minus+'2!N]'
  deltaytitle=!csym.delta_cap+!csym.sigma_cap+$
    ' [h M'+sunsymbol()+' pc!U'+!csym.minus+'2!N]'

  rdeltaytitle=!csym.delta_cap+!csym.sigma_cap+'!D'+!csym.times+$
    '!N [h M'+sunsymbol()+' pc!U'+!csym.minus+'2!N]'
;  randdeltaytitle=!csym.delta_cap+!csym.sigma_cap+'!S!D'+!csym.plus+'!R!Urand'+$
;    '!N [h M'+sunsymbol()+' pc!U'+!csym.minus+'2!N]'
  randdeltaytitle=!csym.delta_cap+!csym.sigma_cap+'!Urand'+$
    '!N [h M'+sunsymbol()+' pc!U'+!csym.minus+'2!N]'

  lumytitle = 'L/area [10!U10!N L'+sunsymbol()+' Mpc!U'+!csym.minus+'2!N]'

  wgmytitle = 'w!Dgm!N  [h!U'+!csym.minus+'1!N Mpc]'

  xigmxtitle = 'r [h!U'+!csym.minus+'1!N Mpc]'
  xigmytitle = !csym.xi+'!Dgm!N(r) '+!csym.times+' '+!csym.omega_cap+'!Dm!N /0.27'

  rhoytitle = $
    !csym.rho+'(r)'+!csym.minus+'!S'+!csym.rho+'!R!A'+!csym.minus+'!N'+$
    '  [10!U12!N M'+sunsymbol()+'  Mpc!U'+!csym.minus+'3!N ]'

  rhoytitle2 = $
    !csym.rho+'(r)'+!csym.minus+'!S'+!csym.rho+'!R!A'+!csym.minus+'!N'+$
    '  [M'+sunsymbol()+' pc!U'+!csym.minus+'3!N ]'

;  rhoytitle = $
;    !csym.rho+'(r)'+!csym.minus+'!S'+!csym.rho+'!R!A'+!csym.minus+'!N'+$
;    '  [h!U2!N M'+sunsymbol()+'  Mpc!U'+!csym.minus+'3!N ]'

  mytitle = 'M(<r) [10!U12!N h!U'+!csym.minus+'1!N M'+sunsymbol()+' ]'
  mytitle2 = 'M(<r) [10!U14!N h!U'+!csym.minus+'1!N M'+sunsymbol()+' ]'
  m200ytitle = 'M!D200!N [10!U14!N h!U'+!csym.minus+'1!N M'+sunsymbol()+' ]'
  r200ytitle = 'r!D200!N [h!U'+!csym.minus+'1!N Mpc]'
  mVirYtitle = 'M!DVir!N [10!U14!N h!U'+!csym.minus+'1!N M'+sunsymbol()+' ]'
  rVirYtitle = 'r!DVir!N [h!U'+!csym.minus+'1!N Mpc]'

  defsysv, '!kpcxtitle',       kpcxtitle
  defsysv, '!kpcxtitle2',      kpcxtitle2
  defsysv, '!mpcxtitle',       mpcxtitle
  defsysv, '!mpcxtitle2',      mpcxtitle2

  defsysv, '!sigytitle',       sigytitle
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
  return

END 
