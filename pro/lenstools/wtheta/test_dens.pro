PRO test_dens, meanr, ndens, ndenserr, areafrac, npair, darea, hit, tarea, tdensity, convert=convert, rec=rec, bx=bx, by=by, fx=fx, fy=fy, dfrom_etarange=dfrom_etarange, lameta=lameta,circle=circle


  ptime, ttt, /savetime
  d2r=!dpi/180d

  rmax = 600./3600.
  rmax=2.
  rmin = 10./3600.
  binsize = 200./3600.
  nbin = long( (rmax-rmin)/binsize ) + 1

  darea = fltarr(nbin)
  R1 = darea & R2 = R1
  FOR i=0L, nbin-1 DO BEGIN 
      R1[i] = rmin + i*binsize
      R2[i] = rmin + (i+1)*binsize
      darea[i] = !pi*(R2[i]^2 - R1[i]^2) ;Mpc^2
  ENDFOR 

  npair = fltarr(nbin)
  ndsum = fltarr(nbin)
  asum = fltarr(nbin)
  rsum = fltarr(nbin)
  hit=0.

  ;; generate bx,by, fx, fy
  IF n_elements(bx) EQ 0 THEN BEGIN 

      maxx = 249.95702
      maxy = 1.2583971

      minx = 130.26052
      miny = -1.2586306
      
      nb = 100000L
      nf = 1000L
      
      tarea = (maxx-minx)*(maxy-miny)
      tdensity = nb/tarea

      bx = arrscl( double(randomu(seed, nb)), minx, maxx, arrmin=0., arrmax=1. )
      by = arrscl( double(randomu(seed, nb)), miny, maxy, arrmin=0., arrmax=1. )
      s=sort(bx)
      bx = bx[s]
      
      fx = arrscl( double(randomu(seed, nf)), minx+rmax, maxx-rmax, arrmin=0., arrmax=1. )
      fy = arrscl( double(randomu(seed, nf)), miny, maxy, arrmin=0., arrmax=1. )
      s=sort(fx)
      fx = fx[s]
      
      dfrommax = maxy - fy
      dfrommin = fy - miny
      ;; convert to lam,eta
      IF keyword_set(lameta) OR keyword_set(convert) THEN BEGIN 
          eq2survey, bx, by, tbx, tby
          s=sort(tbx)
          bx = tbx[s]
          by = tby[s]
          
          eq2survey, fx, fy, tfx, tfy
          s=sort(tfx)
          fx = tfx[s]
          fy = tfy[s]

          IF keyword_set(dfrom_etarange) THEN BEGIN 
              read_etarange, 10, etarange
              wtheta_dfromedge, fx, fy, 1.0, etarange, dfrommax, dfrommin
              dfrommax=dfrommax/d2r & dfrommin=dfrommin/d2r
          ENDIF 
      ENDIF 
      wf=lindgen(nf)

  ENDIF ELSE BEGIN 

      ;; input bx,by
      nb = n_elements(bx)
      nf = n_elements(fx)
      ;; case where input is already lambda/eta
      ;; should be sorted
      IF keyword_set(lameta) THEN BEGIN 

          ;; convert to ra/dec to get dfrom*
          survey2eq, bx, by, tbra, tbdec
          tmaxy = max(tbdec, min=tminy)
          tmaxx = max(tbra, min=tminx)

          tarea = (tmaxx-tminx)*(tmaxy-tminy)
          tdensity = nb/tarea

          IF keyword_set(dfrom_etarange) THEN BEGIN 
              read_etarange, 10, etarange
              wtheta_dfromedge, fx, fy, 1.0, etarange, dfrommax, dfrommin
              dfrommax=dfrommax/d2r & dfrommin=dfrommin/d2r
          ENDIF ELSE BEGIN 
              survey2eq, fx, fy, tfra, tfdec
              dfrommax = tmaxy - tfdec
              dfrommin = tfdec - tminy
          ENDELSE 
      ENDIF ELSE BEGIN 
          ;; input is ra/dec
          maxy = max(by, min=miny)
          maxx = max(bx, min=minx)
          dfrommax = maxy - fy
          dfrommin = fy - miny

          tarea = (maxx-minx)*(maxy-miny)
          tdensity = nb/tarea

          ;; input is ra/dec, but convert to lam,eta
          IF keyword_set(convert) THEN BEGIN 
              
              eq2survey, bx, by, tbx, tby
              s=sort(tbx)
              bx = tbx[s]
              by = tby[s]
              
              eq2survey, fx, fy, tfx, tfy
              s=sort(tfx)
              fx = tfx[s]
              fy = tfy[s]
              
              IF keyword_set(dfrom_etarange) THEN BEGIN 
                  read_etarange, 10, etarange
                  wtheta_dfromedge, fx, fy, 1.0, etarange, dfrommax, dfrommin
                  dfrommax=dfrommax/d2r & dfrommin=dfrommin/d2r
              ENDIF 
          ENDIF 
          
      ENDELSE 

      maxy = max(by, min=miny)
      maxx = max(bx, min=minx)

      wf=where( (fx GT minx+rmax) AND (fx LT maxx-rmax) AND $
                (dfrommin GT 0.) AND (dfrommax GT 0.), nf )


  ENDELSE 

  plot,bx,by,psym=3
  oplot,fx[wf],fy[wf],psym=7,color=!yellow
  ;IF keyword_set(circle) THEN tvcircle, rmax, fx, fy, /data, color=!yellow

  ww=lindgen(nb)
  FOR ii=0L, nf-1 DO BEGIN 
      index=wf[ii]

      IF ((ii MOD 1000) EQ 0) AND (ii NE 0) THEN print,'.',format='(a,$)'

      mxx = fx[index]+rmax
      mnx = fx[index]-rmax

      binary_search, bx, mnx, i1, /round
      binary_search, bx, mxx, i2, /round
      CASE 1 OF
          (i1 NE -1) AND (i2 NE -1): BEGIN
              w=ww[i1:i2]
              num=n_elements(w)
          END 
          (i1 EQ -1) AND (i2 NE -1): BEGIN
              w=ww[0:i2]
              num=n_elements(w)
          END 
          (i2 EQ -1) AND (i1 NE -1): BEGIN
              w=ww[i1:nb-1]
              num=n_elements(w)
          END 
          ELSE: BEGIN
              w=-1
              num=0
          END 
      ENDCASE 

      IF w[0] NE -1 THEN BEGIN 

          IF keyword_set(lameta) OR keyword_set(convert) THEN BEGIN 
              gcirc, 0, fy[index]*d2r, fx[index]*d2r, $
                by[w]*d2r, bx[w]*d2r, $
                dis
              R=dis/d2r
          ENDIF ELSE BEGIN 
              IF keyword_set(rec) THEN BEGIN 
                  dx = fx[index] - bx[w]
                  dy = fy[index] - by[w]
                  R = sqrt(dx^2 + dy^2)
              ENDIF ELSE BEGIN 
                  gcirc, 0, fx[index]*d2r, fy[index]*d2r, $
                    bx[w]*d2r, by[w]*d2r, $
                    dis
                  R=dis/d2r
              ENDELSE 
          ENDELSE 
          
          hist=histogram(R, binsize=binsize, min=rmin, $
                         max=rmax,rever=rev_ind)

;          ntt=long(total(hist))
;          w=where( (R GE rmin) AND (R LE rmax), ntmp)
;          IF ntt NE ntmp THEN print,ntt,ntmp

          numbin=n_elements(hist)
          IF numbin NE nbin THEN BEGIN
              print,'What!!!'
              print,numbin,nbin
              message,' '
          ENDIF 

          whist = where(hist NE 0, nhist)
          ng=0
          ;; Check if there are any in this annulus rmin-rmax
          IF nhist NE 0 THEN BEGIN 

              FOR j=0L, nhist-1 DO BEGIN 

                  binnum = whist[j]

                  barea = binarea(dfrommax[index], $
                                  R1[binnum], $
                                  R2[binnum],$
                                  d_edge2=dfrommin[index],hit=thit)

;                  IF barea LT darea[binnum] THEN print,'!'

                  wb=rev_ind( rev_ind(binnum):rev_ind(binnum+1)-1 )

                  ;; no weights for now
                  np = n_elements(wb)
                  npair[binnum] = npair[binnum] + np
                  ndsum[binnum] = ndsum[binnum] + np/barea
                  asum[binnum] = asum[binnum] + barea/darea[binnum]
                  rsum[binnum] = rsum[binnum] + total(R[wb])
                  ;IF thit GT 1 THEN print,'!'
                  hit=hit+(thit GT 0)

              ENDFOR 
              setzero, R, dx, dy
          ENDIF 

      ENDIF 

  ENDFOR 
  
  ndens = ndsum/nf
  areafrac = asum/nf
  meanr = rsum/npair
  ndenserr = sqrt(npair)/(areafrac*darea)/nf

  print
  ploterror,meanr,ndens,ndenserr,psym=1
  oplot,[0,100],[tdensity,tdensity],color=!red

  ptime, systime(1)-ttt

return
END 
