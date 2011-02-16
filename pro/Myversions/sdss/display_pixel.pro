PRO display_pixel,pixnum,resolution=resolution,over_plot=over_plot,$
		xtitle=xtitle, ytitle=ytitle, $
		lamrange=lamrange,etarange=etarange,rarange=rarange,$
		color=color,radec=radec,_extra=_extra,decrange=decrange,$
		flip_eta=flip_eta,aitoff=aitoff,fill=fill

IF n_elements(color) EQ 0 THEN color = !p.color

IF n_elements(resolution) EQ 0 THEN BEGIN
    resolution = 256
    IF n_elements(silent) EQ 0 THEN BEGIN
        print,'Setting resolution to default 256 (64 pixels across stripe)'
    ENDIF
ENDIF

pix_bound,pixnum,lammin,lammax,etamin,etamax,resolution=resolution,$
          /survey,/silent

IF keyword_set(radec) THEN BEGIN
    csurvey2eq,lammax,etamin,ra1,dec1
    csurvey2eq,lammin,etamin,ra2,dec2
    csurvey2eq,lammin,etamax,ra3,dec3
    csurvey2eq,lammax,etamax,ra4,dec4

    lammin = 0 & lammax = 0 & etamin = 0 & etamax = 0
ENDIF 

IF keyword_set(aitoff) THEN BEGIN
    IF keyword_set(radec) THEN BEGIN
        myaitoff,ra1,dec1,x1,y1
        myaitoff,ra2,dec2,x2,y2
        myaitoff,ra3,dec3,x3,y3
        myaitoff,ra4,dec4,x4,y4
        
        min_x = min([x1,x2,x3,x4])
        max_x = max([x1,x2,x3,x4])
        min_y = min([y1,y2,y3,y4])
        max_y = max([y1,y2,y3,y4])
        
        x1 = 0 & x2 = 0 & x3 = 0 & x4 = 0
        y1 = 0 & y2 = 0 & y3 = 0 & y4 = 0
    ENDIF ELSE BEGIN
        aitoff,etamin,lammax,x_min,y_max
        aitoff,etamax,lammin,x_max,y_min
        
        min_x = min(x_min)
        max_x = max(x_max)
        min_y = min(y_min)
        max_y = max(y_max)
        
        x_min = 0 & x_max = 0 & y_min = 0 & y_max = 0
    ENDELSE
ENDIF

w = where(etamin GT etamax,n_flip)

IF n_flip GT 0 THEN BEGIN
    etamax(w) = etamax(w) + 360.0
ENDIF

IF keyword_set(flip_eta) THEN BEGIN
    w = where(etamax LT 0.0,n_flip)
    IF n_flip GT 0 THEN BEGIN
        etamax(w) = etamax(w) + 360.0
        etamin(w) = etamin(w) + 360.0
    ENDIF 
ENDIF

use_fixed_bound = 0

IF keyword_set(radec) THEN BEGIN 
    IF n_elements(decrange) EQ 0 THEN BEGIN
        IF keyword_set(aitoff) THEN BEGIN
            min_dec = min_y
            max_dec = max_y
        ENDIF ELSE BEGIN
            min_dec = min([dec1,dec2,dec3,dec4])
            max_dec = max([dec1,dec2,dec3,dec4])
        ENDELSE
    ENDIF ELSE BEGIN
        min_dec = decrange(0)
        max_dec = decrange(1)
        use_fixed_bound = 1
    ENDELSE
    
    IF n_elements(rarange) EQ 0 THEN BEGIN
        IF keyword_set(aitoff) THEN BEGIN
            min_ra = min_x
            max_ra = max_x
        ENDIF ELSE BEGIN
            min_ra = min([ra1,ra2,ra3,ra4])
            max_ra = max([ra1,ra2,ra3,ra4])
        ENDELSE
    ENDIF ELSE BEGIN
        min_ra = rarange(0)
        max_ra = rarange(1)
        use_fixed_bound = 1
    ENDELSE
ENDIF ELSE BEGIN
    IF n_elements(lamrange) EQ 0 THEN BEGIN
        IF keyword_set(aitoff) THEN BEGIN
            min_lam = min_y
            max_lam = max_y
        ENDIF ELSE BEGIN
            min_lam = min(lammin)
            max_lam = max(lammax)
        ENDELSE
    ENDIF ELSE BEGIN
        min_lam = lamrange(0)
        max_lam = lamrange(1)
        use_fixed_bound = 1
    ENDELSE
    
    IF n_elements(etarange) EQ 0 THEN BEGIN
        IF keyword_set(aitoff) THEN BEGIN
            min_eta = min_x
            max_eta = max_x
        ENDIF ELSE BEGIN
            min_eta = min(etamin)
            max_eta = max(etamax)
        ENDELSE
    ENDIF ELSE BEGIN
        min_eta = etarange(0)
        max_eta = etarange(1)
        use_fixed_bound = 1
    ENDELSE
ENDELSE    

IF keyword_set(aitoff) EQ 0 OR use_fixed_bound EQ 1 THEN BEGIN
    IF keyword_set(radec) THEN BEGIN 
        pix2ang,pixnum,ra,dec,resolution=resolution,/radec,/silent
        
        w = where(dec LE max_dec AND dec GE min_dec AND $
                  ra LE max_ra AND ra GE min_ra,n_pass)
        
        IF n_pass EQ 0 THEN BEGIN
            print,'No pixels in specified range.  Bailing...'
            RETURN
        ENDIF
        
        ra = 0 & dec = 0
        
        ra1 = ra1(w)
        ra2 = ra2(w)
        ra3 = ra3(w)
        ra4 = ra4(w)

        dec1 = dec1(w)
        dec2 = dec2(w)
        dec3 = dec3(w)
        dec4 = dec4(w)
    ENDIF ELSE BEGIN
        pix2ang,pixnum,lam,eta,resolution=resolution,/survey,/silent
        
        w = where(lam LE max_lam AND lam GE min_lam AND $
                  eta LE max_eta AND eta GE min_eta,n_pass)
        
        IF n_pass EQ 0 THEN BEGIN
            print,'No pixels in specified range.  Bailing...'
            RETURN
        ENDIF
        
        lam = 0 & eta = 0
        
        lammin = lammin(w)
        lammax = lammax(w)
        etamin = etamin(w)
        etamax = etamax(w)
    ENDELSE
ENDIF

if n_elements(xtitle) eq 0 and n_elements(ytitle) eq 0 then begin
	IF keyword_set(radec) THEN BEGIN
		xtitle = textoidl('\alpha')
		ytitle = textoidl('\delta')
	ENDIF ELSE BEGIN
		IF keyword_set(aitoff) THEN BEGIN
			ytitle = textoidl('\lambda')
			xtitle = textoidl('\eta')
		ENDIF ELSE BEGIN
			ytitle = textoidl('\lambda')
			xtitle = textoidl('\eta')
		ENDELSE
	ENDELSE
endif

IF keyword_set(radec) THEN BEGIN
    IF keyword_set(aitoff) THEN BEGIN
        IF n_elements(over_plot) EQ 0 THEN BEGIN
            pplot,[min_ra,max_ra],[min_dec,max_dec],xtitle=xtitle,$
                 ytitle=ytitle,/nodata,/ynozero,_extra=_extra
        ENDIF
        
        FOR i=0L,n_elements(ra1)-1 DO BEGIN
            myaitoff,ra1(i),dec1(i),x1,y1
            myaitoff,ra2(i),dec2(i),x2,y2
            myaitoff,ra3(i),dec3(i),x3,y3
            myaitoff,ra4(i),dec4(i),x4,y4
            tmp_y = [y1,y2,y3,y4,y1]
            tmp_x = [x1,x2,x3,x4,x1]
            IF abs(x1 - x2) LT 10.0 AND abs(x3 - x4) LT 10.0 THEN BEGIN
                IF keyword_set(fill) THEN BEGIN
                    polyfill,tmp_x,tmp_y,color=color
                ENDIF ELSE BEGIN
                    oplot,tmp_x,tmp_y,color=color
                ENDELSE
            ENDIF
        ENDFOR
    ENDIF ELSE BEGIN
        IF n_elements(over_plot) EQ 0 THEN BEGIN
            pplot,[min_ra,max_ra],[min_dec,max_dec],xtitle=xtitle,$
                 ytitle=ytitle,/nodata,/ynozero,_extra=_extra
        ENDIF
        
        FOR i=0L,n_elements(ra1)-1 DO BEGIN
            IF keyword_set(fill) THEN BEGIN
                IF ra1(i) LT min_ra THEN ra1(i) = min_ra
                IF ra2(i) LT min_ra THEN ra2(i) = min_ra
                IF ra3(i) LT min_ra THEN ra3(i) = min_ra
                IF ra4(i) LT min_ra THEN ra4(i) = min_ra
                IF dec1(i) LT min_dec THEN dec1(i) = min_dec
                IF dec2(i) LT min_dec THEN dec2(i) = min_dec
                IF dec3(i) LT min_dec THEN dec3(i) = min_dec
                IF dec4(i) LT min_dec THEN dec4(i) = min_dec

                IF ra1(i) GT max_ra THEN ra1(i) = max_ra
                IF ra2(i) GT max_ra THEN ra2(i) = max_ra
                IF ra3(i) GT max_ra THEN ra3(i) = max_ra
                IF ra4(i) GT max_ra THEN ra4(i) = max_ra
                IF dec1(i) GT max_dec THEN dec1(i) = max_dec
                IF dec2(i) GT max_dec THEN dec2(i) = max_dec
                IF dec3(i) GT max_dec THEN dec3(i) = max_dec
                IF dec4(i) GT max_dec THEN dec4(i) = max_dec

                tmp_ra = [ra1(i),ra2(i),ra3(i),ra4(i),ra1(i)]
                tmp_dec = [dec1(i),dec2(i),dec3(i),dec4(i),dec1(i)]
                
                polyfill,tmp_ra,tmp_dec,color=color
            ENDIF ELSE BEGIN
                tmp_ra = [ra1(i),ra2(i),ra3(i),ra4(i),ra1(i)]
                tmp_dec = [dec1(i),dec2(i),dec3(i),dec4(i),dec1(i)]
                
                oplot,tmp_ra,tmp_dec,color=color
            ENDELSE
        ENDFOR
    ENDELSE
ENDIF ELSE BEGIN
    IF keyword_set(aitoff) THEN BEGIN
        IF n_elements(over_plot) EQ 0 THEN BEGIN
            pplot,[min_eta,max_eta],[min_lam,max_lam],xtitle=xtitle,$
                 ytitle=ytitle,/nodata,/ynozero,_extra=_extra
        ENDIF
        
        FOR i=0L,n_elements(lammin)-1 DO BEGIN
            aitoff,etamin(i),lammin(i),x1,y1
            aitoff,etamax(i),lammin(i),x2,y2
            aitoff,etamax(i),lammax(i),x3,y3
            aitoff,etamin(i),lammax(i),x4,y4
            tmp_y = [y1,y2,y3,y4,y1]
            tmp_x = [x1,x2,x3,x4,x1]
            IF abs(x1 - x2) LT 10.0 AND abs(x3 - x4) LT 10.0 THEN BEGIN
                IF keyword_set(fill) THEN BEGIN
                    polyfill,tmp_x,tmp_y,color=color
                ENDIF ELSE BEGIN
                    oplot,tmp_x,tmp_y,color=color
                ENDELSE
            ENDIF
        ENDFOR
        
        IF keyword_set(bound) THEN BEGIN
            tmp_lam = [min_lam,max_lam,max_lam,min_lam,min_lam]
            tmp_eta = [max_eta,max_eta,min_eta,min_eta,max_eta]
            oplot,tmp_eta,tmp_lam,color=color
        ENDIF
    ENDIF ELSE BEGIN
        IF n_elements(over_plot) EQ 0 THEN BEGIN
            pplot,[min_lam,max_lam],[min_eta,max_eta],xtitle=xtitle,$
                 ytitle=ytitle,/nodata,/ynozero,_extra=_extra
        ENDIF
        
        FOR i=0L,n_elements(lammin)-1 DO BEGIN
            IF keyword_set(fill) THEN BEGIN
                IF lammin(i) LT min_lam THEN lammin(i) = min_lam
                IF lammax(i) GT max_lam THEN lammax(i) = max_lam
                IF etamin(i) LT min_eta THEN etamin(i) = min_eta
                IF etamax(i) GT max_eta THEN etamax(i) = max_eta
                
                tmp_lam = [lammin(i),lammax(i),lammax(i),lammin(i),lammin(i)]
                tmp_eta = [etamax(i),etamax(i),etamin(i),etamin(i),etamax(i)]
                
                polyfill,tmp_lam,tmp_eta,color=color
            ENDIF ELSE BEGIN
                tmp_lam = [lammin(i),lammax(i),lammax(i),lammin(i),lammin(i)]
                tmp_eta = [etamax(i),etamax(i),etamin(i),etamin(i),etamax(i)]
                
                oplot,tmp_lam,tmp_eta,color=color
            ENDELSE
        ENDFOR
    ENDELSE
ENDELSE

RETURN

END
