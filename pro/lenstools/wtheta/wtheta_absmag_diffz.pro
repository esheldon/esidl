PRO wtheta_absmag_diffz, lz, cindex, mag, gmr, absmag, lum_solar, kcorr=kcorr;, old=old

  IF n_params() EQ 0 THEN BEGIN 
      print,'-Syntax: wtheta_absmag_diffz, lz, cindex, mag, gmr, absmag, lum_solar, kcorr=kcorr'
      return
  ENDIF 

  ;; This requires each object have its own z

  ;IF keyword_set(old) THEN GOTO,jump

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; First check/create common block
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  COMMON wtheta_absmag_block, karray, gmr_color, zval, default
  
  nobj=n_elements(mag)
  nz = n_elements(lz)
  IF nz NE nobj THEN message,'Must have a z for each object'

  IF n_elements(karray) EQ 0 THEN BEGIN 
      
      default = -1000.

      ;;Make an array in type redshift space which describes the k corrections
      kfile = '/net/cheops2/home/esheldon/kcorr/fukugita_k_corrections_array.fit'
      IF NOT fexist(kfile) THEN BEGIN 
          kfile = '/sdss3/data5/mckay/fukugita_k_corrections_array.fit'
          IF NOT fexist(kfile) THEN BEGIN 
              message,'Cannot find k-correction file'
          ENDIF 
      ENDIF 
      karray=mrdfits(kfile,0)
      ;;This array has three indices:
      ;;	First: Color u,g,r,i,z
      ;;	Second: Redshift bin
      ;;	Third: z, then six types (E,SO,Sab,Sbc,Scd,Im)
      ;;
      ;;	For example karray(2,*,0) is the redshift values for the r data
      ;;		    karray(2,*,2) is the r' k corrections for SO galaxies

      ;;First set up the color information required to type the galaxies
      gmr_color=karray(1,*,1:6)-karray(2,*,1:6)
      zval=reform( karray(1,*,0) )
      gmr_color=reform(gmr_color)
      ;;Now correct for zero redshift color
      gmr_color(*,0)=gmr_color(*,0)+0.77
      gmr_color(*,1)=gmr_color(*,1)+0.68
      gmr_color(*,2)=gmr_color(*,2)+0.66
      gmr_color(*,3)=gmr_color(*,3)+0.52
      gmr_color(*,4)=gmr_color(*,4)+0.48
      gmr_color(*,5)=gmr_color(*,5)+0.20

  ENDIF 
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; First find k-corrections
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  ;;First use color, and redshift to pick type
  type = replicate(default, nobj)
;  absmag = fltarr(nobj)
  absmag = replicate(abs(default), nobj)

  ;; loop over objects
  FOR i=0L, nobj-1 DO BEGIN 
      ;;Find color of each type at this redshift
      
      tz=lz[i]
      ce=interpol(gmr_color(*,0),zval,tz)
      cso=interpol(gmr_color(*,1),zval,tz)
      csab=interpol(gmr_color(*,2),zval,tz)
      csbc=interpol(gmr_color(*,3),zval,tz)
      cscd=interpol(gmr_color(*,4),zval,tz)
      cim=interpol(gmr_color(*,5),zval,tz)
      carr=[ce,cso,csab,csbc,cscd,cim]
  
      IF gmr[i] ne 0.0 THEN type[i] = interpol(findgen(6), carr, gmr[i]) > 0. < 5.

  ENDFOR 

  kcorr = replicate(default, nobj)
;  kcorr = fltarr(nobj)
  tkarray = reform( karray[cindex,*,1:6] )


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; now absmag
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  w=where(type NE default,nw)
  
  IF nw NE 0 THEN BEGIN 

      kcorr[w] = -1.0*interpolate(tkarray, lz[w]/0.05, type[w])

      tmp=angdist_lambda(lz[w], dlum=dlum) ;Mpc
      dlum = dlum*1.e6          ;pc
      absmag[w] = mag[w] - 5.*alog10(dlum/10.) + kcorr[w]
      
      IF n_params() EQ 6 THEN BEGIN 
;          lum_solar = fltarr(nobj)
          lum_solar = replicate(default, nobj)

          sun=[6.38,5.06,4.64,4.53,4.52]
          lum_solar[w]=10.^((absmag[w]-sun[cindex])/(-2.5))
          lum_solar=lum_solar/1.e10
      ENDIF
  ENDIF 

END 
