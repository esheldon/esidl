pro shape_smooth,struct,sms,stars=stars

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;The idea is to generate an array which is the median of the desired 
;parameter in a "per field" way....
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if n_params() eq 0 then begin
   print,'-syntax shape_smooth,struct,sms'
   return
endif

COMMON seed,seed

field0=min(struct.field)
fieldl=max(struct.field)

nfields=max(struct.field) - min(struct.field) +1

name='other'                          
smf=create_struct(name=name, $
                  'nobj',lonarr(5), $
                  'e1',dblarr(5),'e2',dblarr(5),$
                  'e1sig',dblarr(5), 'e2sig', dblarr(5), $
                  'e1err',dblarr(5),'e2err',dblarr(5),$
                  'e1n1',dblarr(5),'e2n1',dblarr(5),$
                  'e1n2',dblarr(5),'e2n2',dblarr(5),$
                  'rawe1',dblarr(5),'rawe2',dblarr(5),$
                  'rawe1sig',dblarr(5),'rawe2sig',dblarr(5),$
                  'rawe1err',dblarr(5),'rawe2err',dblarr(5),$
                  'ra',double(0.),'dec',double(0.),$
                  'rarange',dblarr(2),'decrange',dblarr(2))
sms=replicate(smf,nfields)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Now randomize the fields both togther and in an uncorrelated way
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

nobj=n_elements(struct.e1(2))

; All colors sorted same.
ssort=randomu(seed,nobj)
si=sort(ssort)
nf_all=struct(si).field

nf=indgen(5,nobj)

; Each color sorted differently
ssort=randomu(seed,nobj)
si=sort(ssort)
nf(0,*)=struct(si).field

ssort=randomu(seed,nobj)
si=sort(ssort)
nf(1,*)=struct(si).field

ssort=randomu(seed,nobj)
si=sort(ssort)
nf(2,*)=struct(si).field

ssort=randomu(seed,nobj)
si=sort(ssort)
nf(3,*)=struct(si).field

ssort=randomu(seed,nobj)
si=sort(ssort)
nf(4,*)=struct(si).field

FOR c=1,3,1 DO BEGIN
                                ; If not looking at stars, get
                                ; stuff that can be dilution corrected.
  IF keyword_set(stars) THEN BEGIN 
    wc = where(struct.e1[c] NE 1.e10 AND struct.e2[c] NE 1.e10, nwc)
  ENDIF ELSE BEGIN 
    wc = where(struct.r[c] NE -1. AND struct.r[c] NE 1.0 AND $
               struct.e1[c] NE 1.e10 AND struct.e2[c] NE 1.e10, nwc)
  ENDELSE 

  IF nwc GT 1 THEN BEGIN 

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Correct only galaxies for dilution.
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    IF keyword_set(stars) THEN R=0. ELSE R=struct[wc].r[c]
    e1corr=struct[wc].e1[c]/(1.-R)
    e2corr = struct[wc].e2[c]/(1.-R)
    errcorr = struct[wc].momerr[c]/(1.-R)
    ellip = sqrt(e1corr^2 + e2corr^2)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Raw shapes
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    rawe1=((struct[wc].ixx[c]-struct[wc].iyy[c])/ $
           (struct[wc].ixx[c]+struct[wc].iyy[c]))
    rawe2=((2.0*struct[wc].ixy[c])/ $
           (struct[wc].ixx[c]+struct[wc].iyy[c]))
    rawerr=struct[wc].momerr[c]
    rawellip=sqrt( rawe1^2 + rawe2^2 )

    ; loop over all the fields
    FOR i=field0,fieldl,1 DO BEGIN

      index=i-field0
      ; Get position information only once
      IF c EQ 2 THEN BEGIN 
        w = where(struct.field eq i)
        med=median(struct[w].ra)
        sms[index].ra=med
        med=median(struct[w].dec)
        sms[index].dec=med
        
        sms[index].rarange[0]=min(struct[w].ra)
        sms[index].rarange[1]=max(struct[w].ra)
  
        sms[index].decrange[0]=min(struct[w].dec)
        sms[index].decrange[1]=max(struct[w].dec)
      ENDIF 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; Shapes for stuff in this field with good measurements.
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      w=where(struct[wc].field EQ i AND ellip LE 1., nw)

                                ;make default value stand out.
      IF nw LE 1 THEN BEGIN
        sms[index].nobj[c]=  -10
        sms[index].e1[c]=    -10.0
        sms[index].e1err=    -10.0
        sms[index].e2[c]=    -10.0
        sms[index].e2err[c]= -10.0
        sms[index].rawe1[c]= -10.0
        sms[index].rawe2[c]= -10.0
        
      ENDIF ELSE BEGIN

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ; Number of good objects found in this field
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        sms[index].nobj[c]=nw

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ; Average the Corrected Shapes
        ; Use intrinsic spread of 0.32
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        spread = 0.32

;        stde1=sqrt( (moment(e1corr[w]))[1] )
;        stde2=sqrt( (moment(e2corr[w]))[1] )

        sig=sqrt(spread^2 + errcorr[w]^2)

        wmom, e1corr[w], sig, wmean, wsig, werr
        sms[index].e1[c] = wmean
        sms[index].e1err[c] = werr
        sms[index].e1sig[c] = wsig
        
        wmom, e2corr[w], sig, wmean, wsig, werr
        sms[index].e2[c] = wmean
        sms[index].e2err[c] = werr
        sms[index].e2sig[c] = wsig

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ; Average the Raw shapes
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;        sms[index].rawe1[c] = median(rawe1[w])
;        sms[index].rawe2[c] = median(rawe2[w])

        ; estimate intrinsic spread
        v1 = (moment(rawe1[w]))[1] 
        v2 = (moment(rawe2[w]))[1] 

        vr = rawerr[w]^2
        v1int = v1 - median(vr)
        v2int = v2 - median(vr)

        wmom, rawe1[w], sqrt( v1int + vr), wmean, wsig, werr
        sms[index].rawe1[c] = wmean
        sms[index].rawe1sig[c] = wsig
        sms[index].rawe1err[c] = werr

        wmom, rawe2[w], sqrt( v2int + vr), wmean, wsig, werr
        sms[index].rawe2[c] = wmean
        sms[index].rawe2sig[c] = wsig
        sms[index].rawe2err[c] = werr
        
      ENDELSE

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; Shapes in sorted fields. All colors sorted the same.
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      w = where(nf_all eq i and ellip LE 1.,nw)
      IF nw LE 1 THEN BEGIN
        sms[index].e1n1[c]= -10.0
        sms[index].e2n1[c]= -10.0
      ENDIF ELSE BEGIN
        sms[index].e1n1[c] = median(e1corr[w])
        sms[index].e2n1[c] = median(e2corr[w])
      ENDELSE 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; Shapes in sorted fields. Each color sorted differently
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      w = where(nf(c,*) eq i and ellip LE 1.,nw)
      IF nw LE 1 THEN BEGIN
        sms[index].e1n2[c]= -10.0
        sms[index].e2n2[c]= -10.0
      ENDIF ELSE BEGIN
        sms[index].e1n2[c] = median(e1corr[w])
        sms[index].e2n2[c]=median(e2corr[w])
      ENDELSE
    ENDFOR
  ENDIF 
ENDFOR 

return
END  





