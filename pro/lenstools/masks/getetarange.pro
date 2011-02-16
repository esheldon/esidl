PRO getetarange, cat, binstruct, oplotstruct=oplotstruct

  IF n_params() LT 1 THEN BEGIN
      print,'-Syntax: getetarange, cat, binstruct'
      return
  ENDIF 

  d2r = !dpi/180.0d0            ;Change degrees to radians.
  r2d = 180./!dpi               ;radians to degrees

  binst=create_struct('clambda', 0d, $
                      'maxeta', 0d, $
                      'mineta', 0d)

  IF NOT tag_exist(cat,'clambda') THEN BEGIN 
      eq2csurvey, cat.ra, cat.dec, lambda, eta
  ENDIF ELSE BEGIN 
      lambda = cat.clambda
      eta = cat.ceta
  ENDELSE 

  binsize=0.3                   ;degrees

  minlam = min(lambda, max=maxlam)

  hist = histogram(lambda, binsize=binsize, min=minlam, max=maxlam,$
                   rever=rev_ind)
  whist=where(hist NE 0, nhist)

  IF nhist NE 0 THEN BEGIN 
      binstruct = replicate(binst, nhist)

      FOR i=0L, nhist-1 DO BEGIN 
          
          binnum=whist[i]
          w=rev_ind( rev_ind[binnum]:rev_ind[binnum+1]-1 )

          binstruct[i].clambda = mean_check(lambda[w])
          mineta = min(eta[w], max=maxeta)
          binstruct[i].maxeta=maxeta
          binstruct[i].mineta=mineta
      ENDFOR 
  ENDIF 

  simpctable
  plot,lambda,eta,psym=3,/ynozero
  oplot,[binstruct.clambda,binstruct.clambda],[binstruct.maxeta,binstruct.mineta],$
    psym=4,color=!red

  IF n_elements(oplotstruct) NE 0 THEN BEGIN 
      myusersym,'fill_circle'
      oplot, oplotstruct.clambda, oplotstruct.ceta, psym=8, symsize=0.7, color=!orange
  ENDIF 

  ;; now let the user adjust the points

  reply=' '
  read,'Do you want to adjust the points (y/n)?',reply
  IF strlowcase(reply) EQ 'n' THEN return
delvarx,w1,w2,w
  continue1 = 1
  WHILE continue1 DO BEGIN 

      ;; get a point to change
      print,'Click near a point'
      cursor, lamlam, etaeta, /data, /down
      etaeta=double(etaeta)
      lamlam=double(lamlam)

      ld = abs(binstruct.clambda - lamlam)
      w=where(ld EQ min(ld))
      ed1 = abs(binstruct[w].maxeta - etaeta)
      ed2 = abs(binstruct[w].mineta - etaeta)
      IF ed1 LT ed2 THEN chmax=1 ELSE chmax=0

      gcirc, 0, etaeta*d2r, lamlam*d2r, $
                binstruct.maxeta*d2r, binstruct.clambda*d2r, $
                dist1
      gcirc, 0, etaeta*d2r, lamlam*d2r, $
                binstruct.mineta*d2r, binstruct.clambda*d2r, $
                dist2

      min1=min(dist1) & w1=where(dist1 EQ min1)
      min2=min(dist2) & w2=where(dist2 EQ min2)
      minn = min([min1, min2])
      IF minn EQ min1 THEN doone=1 ELSE doone=0

      IF chmax THEN BEGIN
          w=w1
          print,w
          oplot, [binstruct[w].clambda], [binstruct[w].maxeta], $
            psym=4, color=!cyan
      ENDIF ELSE  BEGIN
          w=w2
          print,w
          oplot, [binstruct[w].clambda], [binstruct[w].mineta], $
            psym=4, color=!cyan
      ENDELSE 

      ;; now change the point
      continue2 = 1
      WHILE continue2 DO BEGIN

          print,'Click on a new eta position'
          cursor, newlam, neweta, /data, /down
          newlam=binstruct[w].clambda
          neweta=double(neweta)

          oplot,[newlam],[neweta],psym=4,color=!yellow
          
          reply=' '
          read,'keep this point(y/n)',reply
          IF strlowcase(reply) EQ 'y' THEN continue2=0 $
          ELSE oplot,[newlam],[neweta],psym=4,color=!black

      ENDWHILE 

      IF chmax THEN binstruct[w].maxeta = neweta $
      ELSE binstruct[w].mineta = neweta

      reply=' '
      read,'Do another(y/n)?',reply
      IF strlowcase(reply) EQ 'n' THEN continue1=0

  ENDWHILE 

  plot,lambda,eta,psym=3,/ynozero
  oplot,[binstruct.clambda,binstruct.clambda],[binstruct.maxeta,binstruct.mineta],$
    psym=4,color=!red
      
END 
