PRO wtheta_absmag, lz, cindex, mag, gmr, absmag, lum_solar,kcorr=kcorr;, old=old

  IF n_params() EQ 0 THEN BEGIN 
      print,'-Syntax: wtheta_absmag, lz, cindex, mag, gmr, absmag, lum_solar, kcorr=kcorr'
      return
  ENDIF 

  ;; This only works with same lz for each input object.
  ;; Use wtheta_absmag_diffz for each having different z

  ;IF keyword_set(old) THEN GOTO,jump

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; First check/create common block
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  COMMON wtheta_absmag_block, karray, gmr_color, zval, default
  
  nobj=n_elements(mag)
  
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
  type = replicate(2.0, nobj)
  
  ;;Find color of each type at this redshift
  ce=interpol(gmr_color(*,0),zval,lz)
  cso=interpol(gmr_color(*,1),zval,lz)
  csab=interpol(gmr_color(*,2),zval,lz)
  csbc=interpol(gmr_color(*,3),zval,lz)
  cscd=interpol(gmr_color(*,4),zval,lz)
  cim=interpol(gmr_color(*,5),zval,lz)
  carr=[ce,cso,csab,csbc,cscd,cim]
  
  tkarray = reform( karray[cindex,*,1:6] )

  w=where(gmr NE 0.0, nw)
;help,carr,findgen(6)
  IF nw NE 0 THEN type[w] = interpol(findgen(6), carr, gmr[w]) > 0. < 5.

  kcorr = replicate(default, nobj)
  kcorr[w] = -1.0*interpolate(tkarray, replicate(lz,nw)/0.05, type[w])

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; now absmag
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  tmp=angdist_lambda(lz, dlum=dlum) ;Mpc
  dlum = dlum*1.e6              ;pc
  absmag = replicate(abs(default), nobj)

  absmag[w] = mag[w] - 5.*alog10(dlum/10.0) + kcorr[w]
  
  IF n_params() EQ 6 THEN BEGIN 
      sun=[6.38,5.06,4.64,4.53,4.52]
      lum_solar = replicate(default, nobj)
      lum_solar[w]=10.^((absmag[w]-sun[cindex])/(-2.5))

      w2=where(absmag EQ 0,nw2)
      IF nw2 NE 0 THEN lum_solar[w2]=0.

      lum_solar=lum_solar/1.e10
  ENDIF

END 
