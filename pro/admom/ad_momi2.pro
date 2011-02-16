PRO setbad, ixx, iyy, ixy, rho4, uncert

  ixx=0. & iyy=0. & ixy=0. & rho4=0. & uncert=9.999

  return
END 


PRO ad_momi2, im, incenx, inceny, shiftmax, $
             sky, skysig, $
             ixx, iyy, ixy, rho4, uncert, $
             wcenx, wceny, numiter, whyflag, isum, interp

  IF n_params() LT 5 THEN BEGIN
      print,'-Syntax: ad_momi, im, incenx, inceny, shiftmax, sky, skysig, ixx, iyy, ixy [, rho4, uncert, wcenx, wceny, numiter, whyflag]'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; why flags
  ;;
  ;; 1: sum=0 in centroid calculation
  ;; 2: centroid moved by more than shiftmax
  ;; 3: sum=0 in moments calculation
  ;; 4: ixx <= 0 and iyy <= 0
  ;; 5: detm <= detmtol
  ;; 6: detn <= 0
  ;; 7: second detw check (ixx, iyy) has either one <=0
  ;; 8: maxit reached
  ;; 9: detw=0
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; THIS VERSION DIFFERS FROM AD_MOMI.PRO IN THAT IT DOES _NOT_
  ;; USE THE REBIN FUNCTION TO INTERPOLATE. ALSO DOES LOOPS
  ;; INSTEAD OF USING THE TOTAL FUNCTION WHEN INTERPOLATING
  ;; This interpolates better than ad_momi.pro, but much,much slower
  ;; due to not using rebin

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; notes:  ixx,iyy,ixy must already have default values
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Constant parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  tol = 0.01
  detmtol = 1.e-7
  maxit = 30
  maxexp = 9.                   ;maximum exponent

  interpfac = 4.0

  xinterp = 3.0
  xinterp2 = xinterp^2
  nobj = n_elements(incenx)

  s=size(im)
  nx = s[1]
  ny = s[2]

  index = lindgen(nx*ny)
  xindex = index MOD nx
  yindex = index/nx


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; arrays and such
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  rho4 = fltarr(nobj)
  uncert = fltarr(nobj)
  numiter = lonarr(nobj)
  whyflag = lonarr(nobj)

  wcenx = fltarr(nobj)
  wceny = fltarr(nobj)

  isum = fltarr(nobj)

  interp = bytarr(nobj)

  w=fltarr(2,2)
  m=w
  n=w
                                ;loop over centers
  FOR kk=0, nobj-1 DO BEGIN 

      w[0,0] = max( [1.0, ixx[kk]] ) ;starting moments for weight
      w[1,1] = max( [1.0, iyy[kk]] )
      w[0,1] = ixy[kk]
                                ;;;;;;;;;;;;;;;;;;;;
      e1old = 10.               ;initialize 
      e2old = 10.               ;;;;;;;;;;;;;;;;;;;;;
      xcen = incenx[kk]         
      ycen = inceny[kk]
      xcenorig = xcen           ;starting centroids
      ycenorig = ycen

      tixx = ixx[kk]
      tiyy = iyy[kk]
      tixy = ixy[kk]
      trho4 = 0.
      tuncert = 9.999

      imom = 0                  ;number of interations

      WHILE (imom LT maxit) DO BEGIN ; iterate

          interpflag = 0
          imom = imom+1
                                ;initialize
          sum   = 0.
          sums4 = 0.

          detw = w[0,0]*w[1,1] - w[0,1]^2
          IF detw LE 0. THEN BEGIN
              setbad, tixx, tiyy, tixy, trho4, tuncert
              whyflag[kk] = 9
              GOTO, nextobj     ;saves imbedded if's
          ENDIF 

          w1 = w[0,0]/detw
          w2 = w[1,1]/detw      ;now axis ratios
          w12 = w[0,1]/detw

          IF (w[0,0] LT xinterp OR w[1,1] LT xinterp OR $
              detw LT xinterp2) THEN BEGIN 

              ;; keep track of whether there was *ever*
              ;; interpolation

              interp[kk] = 1
              interpflag = 1
          ENDIF 

          grad = 4.*sqrt( max( [w[0,0], w[1,1]] ) ) ;use 4 sigma neighborhood

          ix1 = long( max( [xcen-grad-0.5, 0.0] ) )
          iy1 = long( max( [ycen-grad-0.5, 0.0] ) )
          ix2 = long( min( [xcen+grad+0.5, float(nx-1)] ) )
          iy2 = long( min( [ycen+grad+0.5, float(ny-1)] ) )

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; First find weighted centroid
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          IF interpflag THEN BEGIN ;interpolate undersampled data
              ;print,'hello1'
              sumx=0.
              sumy=0.
              FOR i=ix1, ix2 DO BEGIN
                  xx = i-xcen
                  xx2 = xx*xx
                  xl = xx-0.375
                  xh = xx+0.375
                  FOR j=iy1, iy2 DO BEGIN 
                      yy=j-ycen
                      yy2=yy*yy
                      yl=yy-0.375
                      yh=yy+0.375
                      
                      expon=xl*xl*w2+yl*yl*w1-2.*xl*yl*w12
                      expon=max( [expon, xh*xh*w2+yh*yh*w1-2.*xh*yh*w12] )
                      expon=max( [expon, xl*xl*w2+yh*yh*w1-2.*xl*yh*w12] )
                      expon=max( [expon, xh*xh*w2+yl*yl*w1-2.*xh*yl*w12] )
                      IF expon LE 9. THEN BEGIN 
                          FOR xxx=xl,xh,0.2499 DO BEGIN
                              xx2=xxx*xxx
                              FOR yyy=yl,yh,0.2499 DO BEGIN
                                  yy2=yyy*yyy
                                  expon=(xx2*w2+yy2*w1-2.*xxx*yyy*w12)
                                  weight=exp(-0.5*expon)
                                  ymod=(im[i,j]-sky[kk])*weight/interpfac
                                  sumx=sumx+ymod*float(i)
                                  sumy=sumy+ymod*float(j)
                                  sum = sum+ymod
                              ENDFOR 
                          ENDFOR
                      ENDIF 
                  ENDFOR 
              ENDFOR 
              IF sum LE 0. THEN BEGIN
                  setbad, tixx, tiyy, tixy, trho4, tuncert
                  whyflag[kk] = 1
                  GOTO, nextobj ;saves imbedded if's
              ENDIF 
              xcen = sumx/sum
              ycen = sumy/sum
          ENDIF ELSE BEGIN 

              wobj = where(xindex GE ix1 AND xindex LE ix2 AND  $
                           yindex GE iy1 AND yindex LE iy2, nwobj)
              IF nwobj EQ 0 THEN BEGIN
                  print,'What 1!'
                  return
              ENDIF 
              xx = xindex[wobj]-xcen
              xx2 = xx^2
              yy = yindex[wobj]-ycen
              yy2 = yy^2

              expon = xx2*w2 + yy2*w1 - 2.*xx*yy*w12
              wexp = where(expon LE 9. AND expon GT 0., nwexp)
              IF nwexp NE 0 THEN BEGIN
                  useind=wobj[wexp]
                  weight = exp(-0.5*expon[wexp])
                  ymod = (im[useind] - sky[kk])*weight

                  sum  = total( ymod )
                  IF (sum LE 0.) THEN BEGIN 
                      setbad, tixx, tiyy, tixy, trho4, tuncert
                      whyflag[kk] = 1
                      GOTO, nextobj ;saves imbedded if's
                  ENDIF 
                                ;;;;;;; New Centers! ;;;;;;;
                  xcen = total( ymod*float(xindex[useind]) )/sum
                  ycen = total( ymod*float(yindex[useind]) )/sum
              ENDIF ELSE BEGIN
                  print,'What 3!'
                  return
              ENDELSE 
          ENDELSE 
          diffx = abs(xcen - xcenorig)
          diffy = abs(ycen-ycenorig)
          IF (diffx GT shiftmax) OR $
                         (diffy GT shiftmax) THEN BEGIN
              setbad, tixx, tiyy, tixy, trho4, tuncert
              print,ntostr([xcen,xcenorig,ycen,ycenorig])
              whyflag[kk] = 2
              GOTO, nextobj ;saves imbedded if's
          ENDIF 

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Now find quadrupole moments around weighted centers
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          
          sum=0.
          ;;;; recalculate the neighborhood (not done by phil)
;          ix1 = long( max( [xcen-grad-0.5, 0.0] ) )
;          iy1 = long( max( [ycen-grad-0.5, 0.0] ) )
;          ix2 = long( min( [xcen+grad+0.5, float(nx-1)] ) )
;          iy2 = long( min( [ycen+grad+0.5, float(ny-1)] ) )

;          wobj = where(xindex GE ix1 AND xindex LE ix2 AND  $
;                       yindex GE iy1 AND yindex LE iy2, nwobj)
;          IF nwobj EQ 0 THEN BEGIN
;              print,'What 4!'
;              return
;          ENDIF 

          IF interpflag THEN BEGIN ;interpolate

              sumxx=0.
              sumyy=0.
              sumxy=0.
              FOR i=ix1, ix2 DO BEGIN
                  xx = i-xcen
                  xx2 = xx*xx
                  xl = xx-0.375
                  xh = xx+0.375
                  FOR j=iy1, iy2 DO BEGIN 
                      yy=j-ycen
                      yy2=yy*yy
                      yl=yy-0.375
                      yh=yy+0.375
                      
                      expon=xl*xl*w2+yl*yl*w1-2.*xl*yl*w12
                      expon=max( [expon, xh*xh*w2+yh*yh*w1-2.*xh*yh*w12] )
                      expon=max( [expon, xl*xl*w2+yh*yh*w1-2.*xl*yh*w12] )
                      expon=max( [expon, xh*xh*w2+yl*yl*w1-2.*xh*yl*w12] )
                      IF expon LE 9. THEN BEGIN 
                          FOR xxx=xl,xh,0.2499 DO BEGIN
                              xx2=xxx*xxx
                              FOR yyy=yl,yh,0.2499 DO BEGIN
                                  yy2=yyy*yyy
                                  expon=(xx2*w2+yy2*w1-2.*xxx*yyy*w12)
                                  weight=exp(-0.5*expon)
                                  ymod=(im[i,j]-sky[kk])*weight/interpfac^2
                                  sumxx=sumxx+xx2*ymod
                                  sumyy=sumyy+yy2*ymod
                                  sumxy=sumxy+xxx*yyy*ymod
                                  sums4=sums4+expon^2*ymod
                                  sum = sum+ymod
                              ENDFOR 
                          ENDFOR
                      ENDIF 
                  ENDFOR 
              ENDFOR 
              IF sum LE 0. THEN BEGIN
                  setbad, tixx, tiyy, tixy, trho4, tuncert
                  whyflag[kk] = 3
                  GOTO, nextobj ;saves imbedded if's
              ENDIF 
              m[0,0]=sumxx/sum
              m[1,1]=sumyy/sum
              m[0,1]=sumxy/sum
          ENDIF ELSE BEGIN 
              xx = xindex[wobj]-xcen
              xx2 = xx^2
              yy = yindex[wobj]-ycen
              yy2 = yy^2
              xy = xx*yy

              expon = xx2*w2 + yy2*w1 - 2.*xy*w12
              wexp = where(expon LE 9. AND expon GT 0., nwexp)
              IF nwexp NE 0 THEN BEGIN
                  useind=wobj[wexp]
                  weight = exp(-0.5*expon[wexp])
                  ymod = (im[useind] - sky[kk])*weight

                  sum  = total( ymod )
                  IF (sum LE 0.) THEN BEGIN 
                      setbad, tixx, tiyy, tixy, trho4, tuncert
                      whyflag[kk] = 3
                      GOTO, nextobj ;saves imbedded if's
                  ENDIF 
                                ;;;;;;; Moments ;;;;;;;;
                  m[0,0] = total( ymod*xx2[wexp] )/sum
                  m[1,1] = total( ymod*yy2[wexp] )/sum
                  m[0,1] = total( ymod*xy[wexp] )/sum
                  sums4  = total( ymod*expon[wexp]^2 )

              ENDIF  ELSE BEGIN
                  print,'What 6!'
                  return
              ENDELSE 
          ENDELSE 

          IF (m[0,0] LE 0.) AND (m[1,1] LE 0.) THEN BEGIN
              setbad, tixx, tiyy, tixy, trho4, tuncert
              whyflag[kk] = 4
              GOTO, nextobj ;saves imbedded if's
          ENDIF 

          td = w[0,0] + w[1,1]  ;phil uses w here to cause 2 iterations
          e1 = (w[0,0]-w[1,1])/td  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          e2 = 2.*w[0,1]/td        ;; if converged, then do this stuff
          IF (abs(e1-e1old) LT tol) AND (abs(e2-e2old) LT tol) THEN BEGIN
              whyflag[kk] = 0

              ;; use temporary because will have to deal with errors below
              ;; Note under the convergence, the w's, moments of the weight,
              ;; will be twice the moments of the object.
              tixx = w[0,0]     
              tiyy = w[1,1]
              tixy = w[0,1]
              ;;tixx = m[0,0]     
              ;;tiyy = m[1,1]
              ;;tixy = m[0,1]
              trho4 = sums4/sum

              incenx[kk] = xcen ;might want to remove this
              inceny[kk] = ycen
              wcenx[kk] = xcen
              wceny[kk] = ycen
              numiter[kk] = imom

              ;; Calculate uncertainty
              detw = ( w[0,0]*w[1,1] - w[0,1]^2 )^.25
          
              sdiff = 4.*sum-sums4
              IF sdiff GT 0. THEN BEGIN 
                  tuncert = 4.*sqrt(!pi)*skysig[kk]*detw/sdiff
              ENDIF ELSE BEGIN
                  tuncert = 9.999
              ENDELSE 
              isum[kk] = sum    ;can recover sum(wi*I*xx) by multiplying
              GOTO, nextobj     ;ixx by isum

          ENDIF ELSE BEGIN      ;not yet converged
              detm = m[0,0]*m[1,1] - m[0,1]^2
              IF detm LE detmtol THEN BEGIN
                  setbad, tixx, tiyy, tixy, trho4, tuncert
                  whyflag[kk] = 5
                  GOTO, nextobj ;saves imbedded if's
              ENDIF 
              
              detm = 1./detm
              detw = 1./detw

              n[0,0] =  m[1,1]*detm - w[1,1]*detw
              n[1,1] =  m[0,0]*detm - w[0,0]*detw
              n[0,1] = -m[0,1]*detm + w[0,1]*detw
              detn   =  n[0,0]*n[1,1] - n[0,1]^2

              IF detn LE 0. THEN BEGIN 
                  setbad, tixx, tiyy, tixy, trho4, tuncert
                  whyflag[kk] = 6
                  GOTO, nextobj ;saves imbedded if's
              ENDIF 

              detn = 1./detn
              w[0,0] =  n[1,1]*detn
              w[1,1] =  n[0,0]*detn
              w[0,1] = -n[0,1]*detn
              e1old  = e1
              e2old  = e2
          ENDELSE 
          IF (w[0,0] LT 0.) OR (w[1,1] LT 0) THEN BEGIN 
              setbad, tixx, tiyy, tixy, trho4, tuncert
              whyflag[kk] = 7
              GOTO, nextobj     ;saves imbedded if's
          ENDIF 
      ENDWHILE                  ;While not maxit

      IF imom EQ maxit THEN BEGIN 
          setbad, tixx, tiyy, tixy, trho4, tuncert
          numiter[kk] = maxit
          whyflag[kk] = 8
          GOTO, nextobj         ;saves imbedded if's
      ENDIF 

      nextobj:                  ;used GOTO to save imbedded if's

      ixx[kk] = tixx
      iyy[kk] = tiyy
      ixy[kk] = tixy
      rho4[kk] = trho4
      uncert[kk] = tuncert

  ENDFOR                        ;loop over centers

  return
END 
