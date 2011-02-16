PRO setbad, ixx, iyy, ixy, rho4, uncert

  ixx=0. & iyy=0. & ixy=0. & rho4=0. & uncert=9.999

  return
END 


PRO ad_momi, im, incenx, inceny, shiftmax, $
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
  ;; 7: either ixx, iyy <=0 (second detw check)
  ;; 8: maxit reached
  ;; 9: detw=0
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; notes:  ixx,iyy,ixy must already have default values
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; THIS VERSION USES THE REBIN FUNCTION TO DO THE 
  ;; INTERPOLATION

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Constant parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  tol = 0.01
  detmtol = 1.e-7
;  maxit = 100
  maxit=30
  maxexp = 9.                   ;maximum exponent
  minexp = 0.

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

          wobj = where(xindex GE ix1 AND xindex LE ix2 AND  $
                       yindex GE iy1 AND yindex LE iy2, nwobj)
          IF nwobj EQ 0 THEN BEGIN
              print,'What 1!'
              return
          ENDIF 

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; First find weighted centroid
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          IF interpflag THEN BEGIN ;interpolate undersampled data

              object = im[ix1:ix2, iy1:iy2]
              objsz = size(object)
              objsx = objsz[1]
              objsy = objsz[2]
              newsx = interpfac*objsx
              newsy = interpfac*objsy

                                ;using /sample uses neighbor sampling
              object = rebin(temporary(object), newsx, newsy, /sample)

              newindex = lindgen(newsx*newsy)
              newxindex = newindex MOD newsx
              newyindex = newindex/newsx

              tcx = xcen - ix1
              tcy = ycen - iy1
              newcenx = (tcx + .5)*interpfac - .5
              newceny = (tcy + .5)*interpfac - .5

              xx = newxindex - newcenx
              xx2 = xx^2
              yy = newyindex - newceny
              yy2 = yy^2
                                ;Just interpolating the weight, not the 
                                ;image!
              expon = xx2*w2 + yy2*w1 - 2.*xx*yy*w12
              wexp = where(expon LE maxexp AND expon GT minexp, nwexp)
              IF nwexp NE 0 THEN BEGIN
                  weight = exp(-0.5*expon[wexp])
                  ymod = (object[wexp] - sky[kk])*weight

                  sum  = total( ymod )
                  IF (sum LE 0.) THEN BEGIN 
                      setbad, tixx, tiyy, tixy, trho4, tuncert
                      whyflag[kk] = 1
                      GOTO, nextobj ;saves imbedded if's
                  ENDIF 
                  newcenx = total( ymod*float(newxindex[wexp]) )/sum
                  newceny = total( ymod*float(newyindex[wexp]) )/sum

                  xcen = (newcenx + .5)/interpfac - .5 + ix1
                  ycen = (newceny + .5)/interpfac - .5 + iy1

              ENDIF ELSE BEGIN
                  print,'What 2!'
                  return
              ENDELSE 
          ENDIF ELSE BEGIN 
              xx = xindex[wobj]-xcen
              xx2 = xx^2
              yy = yindex[wobj]-ycen
              yy2 = yy^2

              expon = xx2*w2 + yy2*w1 - 2.*xx*yy*w12
              wexp = where(expon LE maxexp AND expon GT minexp, nwexp)
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

          IF nwobj EQ 0 THEN BEGIN
              print,'What 4!'
              return
          ENDIF 

          IF interpflag THEN BEGIN ;interpolate

              object= im[ix1:ix2, iy1:iy2]
              objsz = size(object)
              objsx = objsz[1]
              objsy = objsz[2]
              newsx = interpfac*objsx
              newsy = interpfac*objsy
                                ;using /sample sets uses neighbor sampling
              object = rebin(temporary(object), newsx, newsy, /sample)

              newindex = lindgen(newsx*newsy)
              newxindex = newindex MOD newsx
              newyindex = newindex/newsx

              tcx = xcen - ix1
              tcy = ycen - iy1
              newcenx = (tcx + .5)*interpfac - .5
              newceny = (tcy + .5)*interpfac - .5

              xx = newxindex - newcenx
              xx2 = xx^2
              yy = newyindex - newceny
              yy2 = yy^2
              xy = xx*yy

              expon = xx2*w2 + yy2*w1 - 2.*xy*w12
              wexp = where(expon LE maxexp AND expon GT minexp, nwexp)
              IF nwexp NE 0 THEN BEGIN
                  weight = exp(-0.5*expon[wexp])
                  ymod = (object[wexp] - sky[kk])*weight/interpfac^2

                  sum  = total( ymod )
                  IF (sum LE 0.) THEN BEGIN 
                      setbad, tixx, tiyy, tixy, trho4, tuncert
                      whyflag[kk] = 1
                      GOTO, nextobj ;saves imbedded if's
                  ENDIF 
                  
                  m[0,0] = total( ymod*xx2[wexp] )/sum
                  m[1,1] = total( ymod*yy2[wexp] )/sum
                  m[0,1] = total( ymod*xy[wexp] )/sum
                  sums4  = total( ymod*expon[wexp]^2 )

              ENDIF ELSE BEGIN
                  print,'What 5!'
                  return
              ENDELSE 
          ENDIF ELSE BEGIN 
              xx = xindex[wobj]-xcen
              xx2 = xx^2
              yy = yindex[wobj]-ycen
              yy2 = yy^2
              xy = xx*yy

              expon = xx2*w2 + yy2*w1 - 2.*xy*w12
              wexp = where(expon LE maxexp AND expon GT minexp, nwexp)
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
