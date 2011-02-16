pro test_error_formulae_rdiff, ntrial, $
                               nobjx, nobjy, xsize, ysize, $
                               sigma, aratio, theta, counts, $
                               e1, e2, e1err, e2err, e1err2, e2err2, $
                               momerr, inpute1, inpute2, $
                               s2n,maxerr=maxerr,etot=etot
  
  IF n_params() LT 9 THEN BEGIN ;Help message
      print,'-Syntax: test_error_formulae, ntrial, counts, sigma, nobjx, nobjy, xsize, ysize, aratio, theta, e1, e2, e1err, e2err, momerr, inpute1, inpute2, s2n'
      return 
  ENDIF 
  
  ;; e1,e2 input
  setup_mystuff
  bindir=!mybindir

  AMOMENT_FAINT =  '200000'X    ; too faint for adaptive moments 
  AMOMENT_SHIFT =  '400000'X    ; centre moved too far while
                                ;determining adaptive moments 
  AMOMENT_MAXITER = '800000'X   ; Too many iterations while
                                ;determining adaptive moments 
  num = 1. - aratio^2
  den = 1. + aratio^2

  inpute1 = (num/den)*cos(2.*theta*!pi/180.)
  inpute2 = (num/den)*sin(2.*theta*!pi/180.)

  sky = 1000.
  skysig = sqrt(sky)

  image = replicate(sky, xsize, ysize)

  ;; calculate x-positions
  
  ;; cneters
  col = arrscl( findgen(nobjx+2), 0., xsize-1 )
  row = arrscl( findgen(nobjy+2), 0., ysize-1 )

  ;; all would be col[0:nobjx+1]
  col = col[1:nobjx]
  row = row[1:nobjy]
  
  ;; add objects 
  delvarx,sendrow, sendcol
  ;;counts = 10000.
  objsize = [51,51]
 
  makegauss, object, objsize, sigma, $
    aratio = aratio, counts=counts, $
    theta=theta,cen=cen

  s2n = counts/sqrt(counts + 16.*!pi*sigma^2*sky)

  print
  print,'Using objects with e1: ',inpute1,' and e2: ',inpute2
  print,'S/N',s2n
  print,'Sky: ',sky,' Skysig: ',skysig
  IF keyword_set(maxerr) THEN print,'Maxerr'
  IF keyword_set(maxerr) AND keyword_set(etot) THEN print,'Using etot'
  print

  FOR ix=0L, nobjx-1 DO BEGIN 
      xmin = (col[ix] - (objsize[0]-1)/2.0)
      xmax = (col[ix] + (objsize[0]-1)/2.0)
      FOR iy=0L, nobjy-1 DO BEGIN 
          ymin = (row[iy] - (objsize[1]-1)/2.0)
          ymax = (row[iy] + (objsize[1]-1)/2.0)

          ;;cen=[(objsize[0]-1.)/2.,(objsize[1]-1.)/2.]
          ;;cen[0] = cen[0]+randomn(seed)/2.
          ;;cen[1] = cen[1]+randomn(seed)/2.
          ;;makegauss, object, objsize, 4., $
          ;;  aratio = aratio, counts=counts, $
          ;;  theta=theta,cen=cen

          FOR iy2=0L, objsize[1]-1 DO BEGIN 
              yy=ymin+iy2
              IF yy LE ysize-1 THEN BEGIN 
                  FOR ix2=0L, objsize[0]-1 DO BEGIN 
                      xx = xmin+ix2
                      IF xx LE xsize-1 THEN BEGIN 
                          image[xx,yy] = image[xx,yy] + object[ix2,iy2]
                      ENDIF 

                  ENDFOR 
              ENDIF 
          ENDFOR 
          add_arrval, col[ix], sendcol
          add_arrval, row[iy], sendrow
          print,'.',format='(a,$)'
      ENDFOR 
  ENDFOR 
  print
  print

  ;; send object to admom

  ;;rdis,image,/silent
  ;;sao,image
  ;;delvarx,image
  ;;return
  
  delvarx, e1, e2, e1err, e2err, e1err2, e2err2,momerr
  FOR i=0L, ntrial-1 DO BEGIN 

      IF ((i+1) MOD 20) EQ 0 THEN print,'Trial: '+ntostr(i+1)+'/'+ntostr(ntrial)

      sendimage=image

      add_noise,sendimage
      
      nrow = long(ysize)
      ncol = long(xsize)
      
      nobj = long(nobjx)*long(nobjy)
      
      ixxpyy = fltarr(nobj)
      ixxmyy = ixxpyy
      i2xy = ixxmyy
      ixxpyyErr = ixxmyy
      ixxmyyErr = ixxmyy
      i2xyErr = ixxmyy

      tmomerr = ixxmyy
      rho4 = ixxmyy
      whyflag = lonarr(nobj)
      
      ;;colprint,sendcol,sendrow

      tmp = call_external(value = [0B,0B,0B,0B,0B,0B,0B,0B,0B,$
                                   0B,0B,0B,0B,0B,0B,0B,0B],$
                          bindir+'admom_rdiff.so','admom_rdiff', $
                          sendimage,$ 
                          nrow, ncol,$
                          sendrow, sendcol, $
                          nobj, sky, skysig,$
                          ixxpyy, ixxmyy, i2xy, $
                          ixxpyyErr, ixxmyyErr, i2xyErr, $
                          tmomerr, rho4, whyflag)
            
      w=where( ((whyflag AND AMOMENT_FAINT) EQ 0) AND $
               ((whyflag AND AMOMENT_SHIFT) EQ 0) AND $
               ((whyflag AND AMOMENT_MAXITER) EQ 0), nw)
      

      te1=fltarr(nobj)
      te2=fltarr(nobj)
      te1err=te1
      te2err=te1
      te1err2=te2
      te2err2=te2

      te1[w] = ixxmyy[w]/ixxpyy[w]
      te2[w] = i2xy[w]/ixxpyy[w]
      
      IF NOT keyword_set(maxerr) THEN BEGIN 
          te1err[w] = ixxmyyErr[w]/ixxpyy[w]*sqrt(1. + 2.*te1[w]^2 )
          te2err[w] = i2xyErr[w]/ixxpyy[w]*sqrt(1. + 2.*te2[w]^2 )
      ENDIF ELSE BEGIN
          IF NOT keyword_set(etot) THEN etot2=1.0 $
          ELSE etot2=te1[w]^2 + te2[w]^2
          te1err[w] = ixxmyyErr[w]/ixxpyy[w]*sqrt(1. + 2.*etot2)
          te2err[w] = i2xyErr[w]/ixxpyy[w]*sqrt(1. + 2.*etot2)
      ENDELSE 
      ;; these should be the same unless e1/e2 are very small
      te1err2[w] = abs(te1[w])*sqrt( (ixxmyyErr[w]/ixxmyy[w])^2 + $
                                     (ixxpyyErr[w]/ixxpyy[w])^2 )
      te2err2[w] = abs(te2[w])*sqrt( (i2xyErr[w]/i2xy[w])^2 + $
                                     (ixxpyyErr[w]/ixxpyy[w])^2 )

      add_arrval, te1, e1
      add_arrval, te2, e2
      add_arrval, te1err, e1err
      add_arrval, te2err, e2err
      add_arrval, te1err2, e1err2
      add_arrval, te2err2, e2err2
      add_arrval, tmomerr, momerr

      setzero, sendimage

  ENDFOR 



END 
