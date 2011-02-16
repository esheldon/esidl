PRO run_spectra_sub, stripes, type

  randNum = 20 + lindgen(10)

  spectra_sub, 'lumthreebinnum_matchLRGNew', randNum

return

;  hirata=1
;  recorr=1

;  stripes1 = [9,10,11,12,13,14,15]
;  stripes2 = [27,28,29,30,31,32,33,34,35,36,37]

;  spectra_sub_meane,stripes1,'lumthreebin',recorr=recorr,hirata=hirata
;  spectra_sub_meane,stripes2,'lumthreebin',recorr=recorr,hirata=hirata

;  spectra_sub_meane,stripes1,'gmr',recorr=recorr,hirata=hirata
;  spectra_sub_meane,stripes2,'gmr',recorr=recorr,hirata=hirata

;  spectra_sub_meane,stripes1,'eclass',recorr=recorr,hirata=hirata,/nooverwrite
;  spectra_sub_meane,stripes2,'eclass',recorr=recorr,hirata=hirata,/nooverwrite

;  spectra_sub_meane,stripes1,'lumthreebinnum',recorr=recorr,hirata=hirata
;  spectra_sub_meane,stripes2,'lumthreebinnum',recorr=recorr,hirata=hirata

;  spectra_sub_meane,stripes1,'vdis',recorr=recorr,hirata=hirata
;  spectra_sub_meane,stripes2,'vdis',recorr=recorr,hirata=hirata

;  spectra_sub_meane,stripes1,'earlylumtwobin',recorr=recorr,hirata=hirata
;  spectra_sub_meane,stripes2,'earlylumtwobin',recorr=recorr,hirata=hirata

;  spectra_sub_meane,stripes1,'earlylumthreebin',recorr=recorr,hirata=hirata
;  spectra_sub_meane,stripes2,'earlylumthreebin',recorr=recorr,hirata=hirata

;  spectra_sub_meane,stripes1,'redlumtwobin',recorr=recorr,hirata=hirata
;  spectra_sub_meane,stripes2,'redlumtwobin',recorr=recorr,hirata=hirata

;  spectra_sub_meane,stripes1,'redlumthreebin',recorr=recorr,hirata=hirata
;  spectra_sub_meane,stripes2,'redlumthreebin',recorr=recorr,hirata=hirata

;  spectra_sub_meane,stripes1,'zsub',recorr=recorr,hirata=hirata
;  spectra_sub_meane,stripes2,'zsub',recorr=recorr,hirata=hirata

;  spectra_sub_meane,stripes1,'vlim1',recorr=recorr,hirata=hirata
;  spectra_sub_meane,stripes2,'vlim1',recorr=recorr,hirata=hirata
;  spectra_sub_meane,stripes1,'vlim2',recorr=recorr,hirata=hirata
;  spectra_sub_meane,stripes2,'vlim2',recorr=recorr,hirata=hirata
;  spectra_sub_meane,stripes1,'vlim3',recorr=recorr,hirata=hirata
;  spectra_sub_meane,stripes2,'vlim3',recorr=recorr,hirata=hirata

;  spectra_sub_meane,stripes1,'mgtr',recorr=recorr,hirata=hirata
;  spectra_sub_meane,stripes2,'mgtr',recorr=recorr,hirata=hirata
  

return

  nst = n_elements(stripes)

  FOR i=0L, nst-1 DO BEGIN 

      stripe = stripes[i]

      spectra_sub_meane, stripe, type, ext='N1.fit'
      spectra_sub_meane, stripe, type, ext='N2.fit'
      spectra_sub_meane, stripe, type, ext='N3.fit'
      spectra_sub_meane, stripe, type, ext='N4.fit'
      
      IF (stripe EQ 10) OR (stripe EQ 12) OR (stripe EQ 82) THEN BEGIN 
          spectra_sub_meane, stripe, type, ext='N5.fit'
          spectra_sub_meane, stripe, type, ext='N6.fit'
          spectra_sub_meane, stripe, type, ext='N7.fit'
          spectra_sub_meane, stripe, type, ext='N8.fit'
      ENDIF 
      
      IF stripe EQ 12  THEN BEGIN 
          spectra_sub_meane, stripe, type, ext='N50.fit'
      ENDIF 

  ENDFOR 

END 
