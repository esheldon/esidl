pro test_error_formulae, ntrial, $
                         nobjx, nobjy, xsize, ysize, $
                         sigma, aratio, theta, counts, $
                         e1, e2, momerr, inpute1, inpute2, $
                         s2n
  
  IF n_params() LT 9 THEN BEGIN ;Help message
      print,'-Syntax: test_error_formulae, ntrial, nobjx, nobjy, xsize, ysize, sigma, aratio, theta, counts, e1, e2, momerr, inpute1, inpute2, s2n'
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
  
  delvarx, e1, e2, momerr
  FOR i=0L, ntrial-1 DO BEGIN 

      IF ((i+1) MOD 20) EQ 0 THEN print,'Trial: '+ntostr(i+1)+'/'+ntostr(ntrial)

      sendimage=image

      add_noise,sendimage
      
      nrow = long(ysize)
      ncol = long(xsize)
      
      nobj = long(nobjx)*long(nobjy)
      
      ixx = fltarr(nobj)
      iyy = ixx
      ixy = ixx

      tmomerr = ixx
      rho4 = ixx
      whyflag = lonarr(nobj)
      
      ;;colprint,sendcol,sendrow

      tmp = call_external(value = [0B,0B,0B,0B,0B,0B,0B,0B,0B,$
                                   0B,0B,0B,0B,0B,0B,0B,0B,0B],$
                          bindir+'admom.so','admom', $
                          sendimage,$ 
                          ncol, nrow,$
                          sendcol, sendrow, $
                          nobj, sky, skysig,$
                          ixx, iyy, ixy, $
                          tmomerr, rho4, whyflag)
            
      w=where( ((whyflag AND AMOMENT_FAINT) EQ 0) AND $
               ((whyflag AND AMOMENT_SHIFT) EQ 0) AND $
               ((whyflag AND AMOMENT_MAXITER) EQ 0), nw)
      
      size = ixx[w] + iyy[w]
      top = (ixx[w] - iyy[w])
      te1=fltarr(nobj)
      te2=fltarr(nobj)


      te1[w] = top/size
      te2[w] = 2.*ixy[w]/size
      
      add_arrval, te1, e1
      add_arrval, te2, e2
      add_arrval, tmomerr, momerr

      setzero, sendimage

  ENDFOR 



END 
