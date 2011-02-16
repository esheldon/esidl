FUNCTION shapelet::init
  return,1
END 


;; Utility functions

FUNCTION shapelet::nmax2ncoeff, nmax
  ncoeff = long((nmax + 1) * (nmax + 2) / 2)
  return,ncoeff
END 

FUNCTION shapelet::decomp_struct, nmax, dtype, nmax_psf=nmax_psf

  ncoeff = self->nmax2ncoeff(nmax)
  decomp = { dtype:dtype, $
             nmax:nmax,                                $
             scale:-1d,                                $
             coeff:dblarr(ncoeff),                      $
             n1:intarr(ncoeff),                        $
             n2:intarr(ncoeff),                        $
             flux:0d,                                  $
             $
             center: [0d, 0d],                         $ 
             size: [0,0],                              $
             $
             C_p:9999d                                 $
           }


  IF n_elements(nmax_psf) NE 0 THEN BEGIN 
      psf_ncoeff = self->nmax2ncoeff(nmax_psf)
      decomp = $
        create_struct(decomp, $
                      'psf_nmax',  nmax_psf,           $
                      'psf_scale', -1d,                $
                      'psf_coeff',  dblarr(psf_ncoeff), $
                      'psf_n1',    intarr(psf_ncoeff), $
                      'psf_n2',    intarr(psf_ncoeff) )
  ENDIF 

  IF dtype EQ 'approx' THEN BEGIN 
      decomp = create_struct(decomp,'ortho',-1d)
  ENDIF 
  IF dtype EQ 'pcr' THEN BEGIN 
      decomp = create_struct(decomp,'npc',0)
  ENDIF

  return,decomp

END 


function shapelet::positive, x
  return, x > 0
end


FUNCTION shapelet::trim_image, image, sky, skysig, nskysig, $
                 smoothing_window=smoothing_window, $
                 slack=slack, $
                 minx=minx, maxx=maxx, $
                 miny=miny, maxy=maxy, $
                 smoothed_image=smoothed_image, $
                 status=status

  on_error, 2
  status = 1

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: '
      print,'  new = shapelet->trim_image(image, sky, skysig, nskysig, '
      print,'                             smoothing_window=, '
      print,'                             slack=, '
      print,'                             minx=, maxx=, miny=, maxy=, '
      print,'                             smoothed_image=,'
      print,'                             status=)'
      print
      message,'Halting'
  ENDIF 

  ;; The smoothing window
  IF n_elements(window) EQ 0 THEN window=4

  ;; Get image dimensions
  sz=size(image)

  nx = sz[1]
  ny = sz[2]
  ntot = sz[4]

  ;; A 2-d index into the image
  index = lindgen(ntot)
  x = index MOD nx
  y = index/nx

  ;; Smooth the image. /edge_truncate is crucial
  timage = smooth(image, window, /edge_truncate)

  w=where(timage GE sky+nskysig*skysig, nw)

  IF nw EQ 0 THEN return, -1


  IF arg_present(smoothed_image) THEN BEGIN 
      smoothed_image=temporary(timage)
  ENDIF ELSE BEGIN 
      timage = 0
  ENDELSE 

  ;; Get the boundary
  minx = min(x[w], max=maxx)
  miny = min(y[w], max=maxy)

  ;; Add slack pixels
  IF n_elements(slack) NE 0 THEN BEGIN 
      minx = minx - slack[0] > 0
      maxx = maxx + slack[0] < nx-1

      miny = miny - slack[0] > 0
      maxy = maxy + slack[0] < ny-1
  ENDIF 

  new_image=image[minx:maxx, miny:maxy]

  status = 0
  return, new_image

END 



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Function to compute the shapelet coefficients if some new scale were
; to be used, from the current coefficients and shapelet scale. See
; the appendix of Refregier (MNRAS, 2003) for more detail.
;
; Author : Brandon C. Kelly, Steward Obs., Feb 2005
;
; Inputs :
;   COEFS - The shapelet coefficients for shapelet of scale SCALE1.
;   SCALE1 - The scale corresponding to COEFS
;   SCALE2 - The new scale that one wants the shapelet coefficients
;            for.
;   N10, N20 - The N1 and N2 arrays that contain the order of the
;              shapelet coefficients. More specifically,
;              (N10[j],N20[j]) is the order of the jth shapelet. These
;              vectors should be the same size as COEFS.
; Output :
;   The shapelet coefficients corresponding to shapelets of scale
;   SCALE2.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function shapelet::scale_transform, coefs, scale1, scale2, n10, n20

  if n_params() lt 5 then begin
      print, 'Syntax- new_coefs = shapelet->scale_transform( coefs, scale1, '
      print, '                                               scale2, n10, n20 )'
      return, 0
  endif

  nmax = max(n10 + n20)         ;get max shapelet order

  b1 = (scale1^2 - scale2^2) / (scale1^2 + scale2^2)
  b2 = 2.0 * scale1 * scale2 / (scale1^2 + scale2^2)

  T = dblarr(nmax+1, nmax+1)    ;the 1-d transformation matrix

;create the transformation matrices
  for n1 = 0, nmax do begin
      
      for n2 = 0, nmax do begin
                                ;n1 and n2 must be both odd or even,
          IF ((n1 MOD 2) EQ (n2 MOD 2)) THEN BEGIN
              
              case 1 of         ;set up l vector
                  min([n1, n2]) eq 0 : l = 0
                  min([n1, n2]) eq 1 : l = 1
                  (n1 mod 2) eq 0 : l = 2 * indgen( min([n1, n2]) / 2 + 1 ) ;l must be even
                  (n1 mod 2) eq 1 : l = 2 * indgen( ceil( min([n1, n2]) / 2.0 ) ) + 1 ;l must be odd
              endcase
              
              entry = (-1.0)^( (n2 - l) / 2.0 ) * sqrt( factorial(n1) * factorial(n2) ) / $
                ( factorial( (n1 - l) / 2.0 ) * factorial( (n2 - l) / 2.0 ) * factorial(l) ) * $
                (b1 / 2.0)^((n1 + n2) / 2.0 - l) * b2^(l + 0.5)
              
              T[n2,n1] = total(entry, /nan) ;n2 indexes the column, n1 the row

          ENDIF
          
      endfor
      
  endfor

;get transformation matrix for 2-d from the 1-d matrices
  T1 = T[*,n10]
  T1 = T1[n10,*]
  T2 = T[*,n20]
  T2 = T2[n20,*]

  T = T1 * T2

;now transform the coefficients to new scale
  new_coefs = T ## transpose(coefs)

  return, reform(new_coefs)
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Function to compute the flux of a galaxy from its shapelet
; coefficients.
;
; Author : B.C. Kelly, Steward Obs., Feb 2005
;
; Inputs :
;   COEFS - The shapelet coefficients.
;   SCALE - The scale of the shapelet functions.
;   N1, N2 - The N1 and N2 arrays that contain the order of the
;            shapelet coefficients. More specifically,
;            (N1[j],N2[j]) is the order of the jth shapelet. These
;            vectors should be the same size as COEFS.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function shapelet::flux, coefs, scale, n1, n2

  if n_params() lt 4 then begin
      print, 'Syntax- flux = shapelet->flux( coefs, scale, n1, n2 )'
      return, 0
  endif

  m = n_elements(coefs)

  sum = 0d
  for i = 0, m - 1 do begin
      
      if (n1[i] mod 2 eq 0) and (n2[i] mod 2 eq 0) then sum = sum + $
        2.0^((2.0 - n1[i] - n2[i]) / 2.0) * sqrt( factorial(n1[i]) * factorial(n2[i]) ) / $
        ( factorial(n1[i] / 2.0) * factorial(n2[i] / 2.0) ) * coefs[i]
      
  endfor

  flux = sqrt(!pi) * scale * sum

  return, flux
end



pro shapelet::split_array, ind, ind1, ind2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Return the IND2 indices of an array when the IND1 indices are
; removed from the original indices IND.  Note that
; IND = sort([IND1,IND2])
;
; Author : B.C.Kelly, Steward Observatory, April 2004
;
; Inputs :
;   IND - The original indices.
;   IND1 - The indices to remove from IND
;
; Outputs : 
;   IND2 - The remaining indices, after IND1 has been removed from
;          IND.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if n_params() lt 3 then begin
    print, 'Syntax- shapelet->split_array, ind, ind1, ind2'
    return
endif

ind2 = [-1]
for i = 0, n_elements(ind) - 1 do begin
    dum = where(ind1 eq ind[i], nmatch)
    if nmatch eq 0 then ind2 = [ind2, i]
endfor

ind2 = n_elements(ind2) gt 1 ? ind2[1:*] : ind2

return
end




PRO shapelet::resize_image, image, x0, sz

; Author: Brandon Kelly June 2002
; Routine to resize and image, centered around a point x0, with a given
; size (sz)

if n_params() lt 3 then begin
    print, 'Syntax- shapelet->resize_image, image, x0, sz'
    return
endif

sz1=sz[0]
sz2=sz[1]
nx=size(image)

IF x0[0]-sz1/2 LE 0 THEN BEGIN
    newx=fix(sz1/2-x0[0])
    image1=fltarr(nx[1]+newx, nx[2])
    if n_elements(image[*,0]) ne nx[1] or $
      n_elements(image[0,*]) ne nx[2] then begin
        print, 'Size of Image does not match induces for Image 1, returning...'
        return
    endif else image1[newx : nx[1]+newx-1, 0 : nx[2]-1]=image
ENDIF
IF x0[0]-sz1/2 GT 0 THEN BEGIN
    newx=fix(x0[0]-sz1/2)
    if n_elements(image[*,0]) lt nx[1] or $
      n_elements(image[0,*]) lt nx[2] then begin
        print, 'Size of Image is too small compared to Image 1, returning...'
        return
    endif else image1=image[newx : nx[1]-1, 0:nx[2]-1]
ENDIF
IF nx[1]-1-x0[0]-sz1/2 LE 0 THEN BEGIN
    nx1=size(image1)
    newx=fix(sz1/2-(nx[1]-1-x0[0]))
    image2=fltarr(nx1[1]+newx,nx[2])
    if n_elements(image1[*,0]) ne nx1[1] or $
      n_elements(image1[0,*]) ne nx1[2] then begin
        print, 'Size of Image 1 does not match induces for Image 2, returning...'
        return
    endif else image2[0:nx1[1]-1, 0:nx1[2]-1]=image1
ENDIF
IF nx[1]-1-x0[0]-sz1/2 GT 0 THEN BEGIN
    nx1=size(image1)
    newx=fix(nx[1]-1-x0[0]-sz1/2)
    if n_elements(image1[*,0]) lt nx1[1]-newx or $
      n_elements(image1[0,*]) lt nx1[2] then begin
        print, 'Size of Image 1 is too small compared to Image 2, returning...'
        return
    endif else image2=image1[0:nx1[1]-newx-1,0:nx1[2]-1]
ENDIF
IF x0[1]-sz2/2 LE 0 THEN BEGIN
    nx2=size(image2)
    newy=fix(sz2/2-x0[1])
    image3=fltarr(nx2[1], nx2[2]+newy)
    if n_elements(image2[*,0]) ne nx2[1] or $
      n_elements(image2[0,*]) ne nx2[2] then begin
        print, 'Size of Image 2 does not match induces for Image 3, returning...'
        return
    endif else image3[0 : nx2[1]-1, newy : newy+nx2[2]-1]=image2
ENDIF
IF x0[1]-sz2/2 GT 0 THEN BEGIN
    nx2=size(image2)
    newy=fix(x0[1]-sz2/2)
    if n_elements(image2[*,0]) lt nx2[1] or $
      n_elements(image2[0,*]) lt nx2[2] then begin
        print, 'Size of Image 2 is too small compared to Image 3, returning...'
        return
    endif else image3=image2[0 : nx2[1]-1, newy : nx2[2]-1]
ENDIF
IF nx[2]-1-x0[1]-sz2/2 LE 0 THEN BEGIN
    nx3=size(image3)
    newy=fix(sz2/2-(nx[2]-1-x0[1]))
    image4=fltarr(nx3[1],nx3[2]+newy)
    if n_elements(image3[*,0]) ne nx3[1] or $
      n_elements(image3[0,*]) ne nx3[2] then begin
        print, 'Size of Image 3 does not match induces for Image 4, returning...'
        return
    endif else image4[0:nx3[1]-1, 0:nx3[2]-1]=image3
ENDIF
IF nx[2]-1-x0[1]-sz2/2 GT 0 THEN BEGIN
    nx3=size(image3)
    newy=fix(nx[2]-1-x0[1]-sz2/2)
    if n_elements(image3[*,0]) lt nx3[1]-1 or $
      n_elements(image3[0,*]) lt nx3[2]-newy-1 then begin
        print, 'Size of Image 3 is too small compared to Image 4, returning...'
        return
    endif else image4=image3[0:nx3[1]-1, 0:nx3[2]-newy-1]
ENDIF
image=image4

return
END




;+
; NAME:
;   CMREPLICATE
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;
; PURPOSE:
;   Replicates an array or scalar into a larger array, as REPLICATE does.
;
; CALLING SEQUENCE:
;   ARRAY = CMREPLICATE(VALUE, DIMS)
;
; DESCRIPTION: 
;
;   The CMREPLICATE function constructs an array, which is filled with
;   the specified VALUE template.  CMREPLICATE is very similar to the
;   built-in IDL function REPLICATE.  However there are two
;   differences:
;
;      * the VALUE can be either scalar or an ARRAY.
;
;      * the dimensions are specified as a single vector rather than
;        individual function arguments.
;
;   For example, if VALUE is a 2x2 array, and DIMS is [3,4], then the
;   resulting array will be 2x2x3x4.
;
; INPUTS:
;
;   VALUE - a scalar or array template of any type, to be replicated.
;           NOTE: These two calls do not produce the same result:
;                  ARRAY = CMREPLICATE( 1,  DIMS)
;                  ARRAY = CMREPLICATE([1], DIMS)
;           In the first case the output dimensions will be DIMS and
;           in the second case the output dimensions will be 1xDIMS
;           (except for structures).  That is, a vector of length 1 is
;           considered to be different from a scalar.
;
;   DIMS - Dimensions of output array (which are combined with the
;          dimensions of the input VALUE template).  If DIMS is not
;          specified then VALUE is returned unchanged.
;
; RETURNS:
;   The resulting replicated array.  
;
; EXAMPLE:
;   x = [0,1,2]
;   help, cmreplicate(x, [2,2])
;     <Expression>    INT       = Array[3, 2, 2]
;   Explanation: The 3-vector x is replicated 2x2 times.
;
;   x = 5L
;   help, cmreplicate(x, [2,2])
;     <Expression>    LONG      = Array[2, 2]
;   Explanation: The scalar x is replicated 2x2 times.
;
; SEE ALSO:
;
;   REPLICATE
;
; MODIFICATION HISTORY:
;   Written, CM, 11 Feb 2000
;   Fixed case when ARRAY is array of structs, CM, 23 Feb 2000 
;   Apparently IDL 5.3 can't return from execute().  Fixed, CM, 24 Feb
;     2000
;   Corrected small typos in documentation, CM, 22 Jun 2000
;
;-
; Copyright (C) 2000, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
;-
function shapelet::cmreplicate, array, dims

  on_error, 2
  if n_params() EQ 0 then begin
      message, 'RARRAY = shapelet->CMREPLICATE(ARRAY, DIMS)', /info
      print
      message,'Halting'
  endif

      
  if n_elements(dims) EQ 0 then return, array
  if n_elements(array) EQ 0 then $
    message, 'ERROR: ARRAY must have at least one element'

  ;; Construct new dimensions, being careful about scalars
  sz = size(array)
  type = sz(sz(0)+1)
  if sz(0) EQ 0 then return, make_array(value=array, dimension=dims)
  onedims = [sz(1:sz(0)), dims*0+1] ;; For REFORM, to extend # of dims.
  newdims = [sz(1:sz(0)), dims    ] ;; For REBIN, to enlarge # of dims.
  nnewdims = n_elements(newdims)
  
  if nnewdims GT 8 then $
    message, 'ERROR: resulting array would have too many dimensions.'

  if type NE 7 AND type NE 8 AND type NE 10 AND type NE 11 then begin
      ;; Handle numeric types

      ;; Argghh!  Why doesn't REBIN take an *array* of dimensions!
      ;; *Instead* we need to run EXECUTE(), with a string like this:
      ;; rebin(array1, newdims(0), newdims(1), ...)
      ;; That's what the following format string does.
      fmt = '('+strtrim(nnewdims,2)+'("newdims(",I0,")",:,","))'
      arglist = string(lindgen(nnewdims), format=fmt)
      cmd = 'array1 = rebin(reform([array], onedims),'+arglist+')'
      array1 = 0
      retval = execute(cmd)
      if retval EQ 1 then return, array1
      
      ;; If execution reaches here then an error occurred.
      message, 'ERROR: array could not be resized'
      return, 0L

  endif else begin
      ;; Handle strings, structures, pointers, and objects separately
      
      ;; Handle structures, which are never scalars
      if type EQ 8 AND sz(0) EQ 1 AND n_elements(array) EQ 1 then $
        return, reform(make_array(value=array, dimension=dims), dims, /over)

      nold = n_elements(array)
      nadd = 1L
      for i = 0L, n_elements(dims)-1 do nadd = nadd * round(dims(i))
      if type EQ 8 then $
        array1 = make_array(value=array(0), dimension=[nold,nadd]) $
      else $
        array1 = make_array(type=type, dimension=[nold,nadd], /nozero)
      array1 = reform(array1, [nold,nadd], /overwrite)
      array2 = reform([array], n_elements(array))

      ;; Efficient copying, done by powers of two
      array1(0,0) = temporary(array2)
      stride = 1L   ;; stride increase by a factor of two in each iteration
      i = 1L & nleft = nadd - 1
      while nleft GT stride do begin
          array1(0,i) = array1(*,0:stride-1)  ;; Note sneaky IDL optimization
          i = i + stride & nleft = nleft - stride
          stride = stride * 2
      endwhile
      if nleft GT 0 then array1(0,i) = array1(*,0:nleft-1)

      return, reform(array1, newdims, /overwrite)
  endelse

end




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Function to construct the Gauss-Hermite (shapelet) basis. The rows of
; this matrix contain the basis functions, e.g., GH[*,j] contains the
; jth shapelet mode sampled at the NX data points. Note that NX and X0
; are 2-element vectors.
;
; Author : Brandon C. Kelly, Steward Obs., Feb 2005
;
; Inputs :
;    NX - The dimensions of the shapelet basis array. The shapelet
;         basis images will be NX[0] x NX[1] arrays.
;    X0 - The centroid of the shapelets, a 2-element vector.
;    SCALE - The scale of the shapelets.
;    NMAX - The maximum order of the shapelets. All shapelet states
;           for which n1 + n2 =< NMAX will be generated.
;
; Output :
;    The shapelet basis, an [NX[0] * NX[1], M] array, where M is the
;    number of shapelet states. 2-d arrays of the shapelet functions
;    can be generated from the rows of GH as reform(GH[*,j], NX[0],
;    NX[1]).
;
; Optional Outputs :
;   N1, N2 - M-element array containing the order of the
;            shapelets. GH[*,j] contains the shapelet of order
;            (N1[j],N2[j]).
;
; Functions Called :
;   PHIN
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function shapelet::construct_gh_matrix, nx, x0, scale, nmax, n1=n1, n2=n2, $
                 status=status

  status = 1
  if n_params() lt 4 then begin
      print, 'Syntax- GH = shapelet->construct_gh_matrix( nx, x0, scale, nmax, n1=, n2=, status=)'
      return, 0
  endif

  m = long((nmax + 1) * (nmax + 2) / 2) ;number of shapelet states

  n1 = intarr(m)
  n2 = n1

  i = 0L
  for n = 0, nmax do begin      ;creat the index vectors, n1 and n2
      for n1i = 0, n do begin
          n1[i] = n1i
          n2[i] = n - n1i
          i = i + 1
      endfor
  endfor

  x1 = findgen(nx[0]) - x0[0]   ;the absissca in the horizontal direction
  x2 = findgen(nx[1]) - x0[1]   ;the absissca in the vertical direction

  ngh = long(nx[0])*nx[1]
  GH = dblarr(ngh, m)   ;the shapelet basis
;get the basis vectors and store them in GH
  FOR j = 0, m - 1 DO BEGIN 
      phi_n = self->phin( n1[j], n2[j], x1, x2, [scale,scale] )
      IF n_elements(phi_n) NE ngh THEN BEGIN 
          message,'phi_n is wrong size',/inf
          return, -1
      ENDIF 
      GH[*,j] = temporary(phi_n)
  ENDFOR 
  status = 0
  return, GH
end



function shapelet::hermite, n, x
  
; July 99 - Written by A. Refregier
;
; PURPOSE: compute the hermite polynomial Hn(x) of order n.
; INPUT: n,x: output Hn(x) where x can be a scalar, a vector or an array
; OUTPUT: Hn(x)
  
  nx = n_elements(x)
  
  case n of
     0: begin
          sx = size(x)
          case sx(0) of
             0: h = 1.
             1: h = replicate(1., sx(1))
             2: h = replicate(1., sx(1), sx(2))
          else: begin
                  print, 'hermite: x must be a scalar, a vector or an array'
                end
          endcase
        end
     1: h = 2.*x
     2: h = 4.*x^2-2.
     3: h = 8*x^3-12.*x
     4: h = 16*x^4-48.*x^2+12.
     5: h = 32*x^5-160.*x^3+120.*x
     6: h = 64.*x^6-480.*x^4+720.*x^2-120.
  else: begin
          c = self->hermitecof(n)     
          sx = size(x)
          if sx(0) eq 0 then h = total(c*x^findgen(n+1)) else begin
             h = c(0)+c(1)*x
             for i=2, n do h = h+c(i)*x^i
          endelse
        end
  endcase
     
  
  return, h
  end  


function shapelet::hermitecof, n
  
; July 99 - Written by A. Refregier
;
; PURPOSE: compute the polynomial coefficients for Hn(x), the Hermite
; polynomial of order n
; INPUT: n: order of the Hermite polynomial
; OUTPUT: hermitecof: coefficients ci for Hn(x)=Sum_i ci*x^i, i=0,..,n
  
c = dblarr(n+1)
for s=0, n/2 do $ 
  c(n-2*s) = (-1)^s*2.^(n-2*s)*factorial(n)/(factorial(s)*factorial(n-2*s))

return, c
end
  

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Function to compute the inter-quartile range of a data set. This can
; give a more robust estimation of the standard deviation of a data
; set that is assumed to be approximately normal, as
;              
;                    IQR / 1.34 = Sigma
;                 
; Author : Brandon C. Kelly, Steward Obs., Sept. 2004
;
; INPUTS :
;    X - The data.
;
; OUTPUT :
;    The inter-quartile range of the data, IQR.
;
; OPTIONAL OUTPUT :
;    SIGMA - This will be a robust estimate of the standard deviation
;            assuming that the data set is drawn from a normal
;            distribution. SIGMA will be the smaller of the sample
;            standard deviation and IQR / 1.34.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function shapelet::iqr, x, sigma

  if n_params() eq 0 then begin
      print, 'Syntax- Result = shapelet->iqr( x, [sigma])'
      return, 0
  endif
  
  nx = n_elements(x)
  
  sorted = sort(x)
  
  iqr = x[sorted[3 * nx / 4]] - x[sorted[nx / 4]]
  
  sigma = min( [stddev(x, /nan), iqr / 1.34] )
  
  return, iqr
end

FUNCTION shapelet::phin, n1, n2, x1, x2, scale

; This function returns the basis functions phi_n(x) for the hermite
; polynomials.  It is linked to code written by A. Refregier.
; 
; Author:  Brandon Kelly May 2002
;
; Input:
; n1 = first order indicator
; n2 = second order indicator
; x1 = array containing the values for x_1
; x2 = array containing the values for x_2
; scale= a two-element vector containing the x and y scales
;
; Output: the basis functions phi_n(x)

  s1 = size(x1)
  s2 = size(x2)
  phi_n1 = dblarr(s1[1])
  h_n1 = self->hermite(n1, x1 / scale[0])
  phi_n1 = $
    double( (2.0^n1 * sqrt(!pi) * factorial(n1) * scale[0])^(-0.5) * h_n1 * $
            exp( -(x1 / scale[0])^2 / 2.0) )
  
  phi_n2 = dblarr(s2[1])
  h_n2 = self->hermite(n2, x2 / scale[1])
  phi_n2 = $
    double((2.0^n2 * sqrt(!pi) * factorial(n2) * scale[1])^(-0.5) * h_n2 * $
           exp( -(x2 / scale[1])^2 / 2.0) )
  
  phi = transpose(phi_n2) ## phi_n1
  
  return, phi
end
















;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Function to compute Mallow's C_p, an estimate of the chi-square
; between the true function and the fit. This function assumes that
; both y and yfit have been normalized by the covariance matrix of y.
;
; Author : Brandon C. Kelly, Steward Obs., Feb 2005
;
; Inputs :
;    Y - The data, normalized by the covariance matrix so that it is
;        normally distributed with covariance equal to the identity
;        matrix.
;    YFIT - The fit to Y, normalized by the covariance matrix of Y.
;    M - The number of parameters used in the linear model :
;                  YFIT = S ## Y ,
;        For S some matrix that converts Y to YFIT.
;
; Output :
;    The mean-squared error of YFIT with the true function that
;    generated, i.e., the mean of the distribution that generated Y.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function shapelet::cp_mse, y, yfit, m

  if n_params() ne 3 then begin
      print, 'Syntax- C_p = shapelet->cp_mse( y, yfit, m )'
      return, -1
  endif

  n = n_elements(y)
  cp = mean( (y - yfit)^2, /nan ) - 1 + 2d * m / n
  
  return, cp
end


















;; Decomposition function


FUNCTION shapelet::approx_decomp, $
                 image_struct, image_var, nmax, $
                 ghfit=ghfit, ghmatrix=gh, $
                 max_size=max_size, status=status

  on_error, 2
  status=1

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: '
      print,'  decomp = shapelet->approx_decomp(image_struct, image_var, nmax,'
      print,'                                   max_size=, '
      print,'                                   ghfit=, ghmatrix=, status='
      print
      message,'Halting'
  ENDIF 


  x0 = image_struct.center
  scale = image_struct.scale
  skysig = image_struct.skysig
  image = image_struct.image

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; We will limit the maximum image size due to memory considerations
  ;; E.S.S.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(max_size) EQ 0 THEN max_size = 500

  max_size = max_size[0]
  msstr=ntostr(max_size)
  msstr='['+msstr+', '+msstr+']'

  nx0 = size(image, /dim)
  IF nx0[0] GT max_size OR nx0[1] GT max_size THEN BEGIN 
      print
      print,'Image Size is larger than allowed size: '+msstr
      print,'Exiting approx_decomp'
      print
      message,'Halting'
  ENDIF 

  szfact = 20.0             ;number of scales that the image will be resized to
;  szfact=5
  osize = long(szfact * scale) ; resize the image to this size, add blank sky

  osize = [osize, osize]

  ;; choose the larger of osize and the orignal size
  IF osize[0] LT nx0[0] THEN osize[0] = nx0[0] 
  IF osize[1] LT nx0[1] THEN osize[1] = nx0[1]

  IF osize[0] GT max_size OR osize[1] GT max_size THEN BEGIN 
      print
      msstr=ntostr(max_size[0])
      msstr='['+msstr+', '+msstr+']'
      print,'Expanded Size is larger than allowed size: '+msstr
      print,'Exiting approx_decomp'
      print
      message,'Halting'
  ENDIF 

  ;; making a copy which we will modify
  image0 = image
  image_var0 = image_var
  x00 = x0

  ;; Expand the size if needed
  IF osize[0] NE nx0[0] OR osize[1] NE nx0[1] THEN BEGIN
      self->resize_image, image, x0, osize ;resize the image
      self->resize_image, image_var, x0, osize
      
      nx = size(image, /dim)
      x0 = nx / 2.0             ;new centroid

      w = where(image eq 0, nw) ;add sky noise to blank sky that was added
      if nw ne 0 then begin
          image[w] = image[w] + skysig * randomn(seed,nw)
          image_var[w] = skysig^2
      ENDIF
  ENDIF ELSE nx = nx0
                                ;construct the shapelet
                                ;(Gauss-Hermite) basis
  GH = self->construct_gh_matrix(nx, x0, scale, nmax, n1=n1, n2=n2)

;GH matrix is approximately orthogonal, use this to estimate the
;shapelet coefficients
  coefs = GH ## transpose(image[*])
  coefs = transpose(coefs)

  m = n_elements(coefs)         ;number of parameters

;now find best scale on a grid of scale values from 0.5*scale to 1.5*scale
  nscales = 10
  gammas = scale / 2.0 + dindgen(nscales) / (nscales - 1.0) * scale

  error = fltarr(nscales)

  for k = 0, nscales - 1 do BEGIN
                                ;transform the shapelet coefficients
                                ;to have scale gammas[k]
      coefs_k = self->scale_transform( coefs, scale, gammas[k], n1, n2 )
                                ;construct new basis with this scale
      GH_k = self->construct_gh_matrix(nx0, x00, gammas[k], nmax)
      ghfit = coefs_k ## GH_k
      ghfit = reform(ghfit, nx0[0], nx0[1])

      C_p = self->cp_mse( image0 / sqrt(image_var0), ghfit / sqrt(image_var0), m ) ;estimate of the MSE
      error[k] = C_p

  endfor

;now use a spline to interpolate MSE to get best scale
  nscales2 = 200
  gammas2 = scale / 2.0 + dindgen(nscales2) / (nscales2 - 1.0) * scale
  error2 = interpol(error, gammas, gammas2, /spline)
  minerr = min(error2, minind)  ;find minimum error
;get shapelet coefficients for this scale
  coefs_k = self->scale_transform( coefs, scale, gammas2[minind], n1, n2 )
  GH_k = self->construct_gh_matrix(nx0, x00, gammas2[minind], nmax)

;find MSE in model with best scale
  ghfit = coefs_k ## GH_k
  ghfit = reform(ghfit, nx0[0], nx0[1])

;use Mallow's C_p statistic to estimate the mean squared error (MSE)
  C_p = self->cp_mse( image0 / sqrt(image_var0), ghfit / sqrt(image_var0), m )

  image = image0                ;resore inputs
  x0 = x00
  decomp = {coeff:coefs_k, scale:gammas2[minind], n1:n1, n2:n2, C_p:C_p}

  status=0
  return, decomp
end





;Function to estimate the shapelet coefficients via principal
;component regression.
FUNCTION shapelet::pcr_decomp, image_struct, nmax, $
                 image_var=image_var, $
                 psf_struct=psf_struct, $
                 ghfit=ghfit, ghmatrix=gh, svthresh=svthresh

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: '
      print,'  shapelet->pcr_decomp(image_struct, nmax, '
      print,'                       image_var=, '
      print,'                       psf_struct=,'
      print,'                       ghfit=, ghmatrix=, svthresh=)'
      print
      message,'Halting'
  ENDIF 

  x0    = image_struct.center
  scale = image_struct.scale
  image = image_struct.image

  nx = size(image, /dim)
  if n_elements(image_var) eq 0 then image_var = replicate(1.0, nx[0], nx[1])
  if n_elements(svthresh) eq 0 then begin
      svthresh = 0.0 
      auto = 1
  endif else auto = 0

  ;; construct the shapelet (Gauss-Hermite) basis
  IF n_elements(psf_struct) NE 0 THEN BEGIN 
      scale = sqrt(self->positive(scale^2 - psf_struct.scale^2)) > 1.0
  ENDIF 
  GH = self->construct_gh_matrix(nx, x0, scale, nmax, n1=n1, n2=n2)

  m = n_elements(n1)            ;number of coefficients to estimate

  image = (image / sqrt( image_var ))[*] ;normalize by the noise and make into a vector

  if n_elements(psf_struct) ne 0 then begin
                                ;convolve shapelet basis with the PSF
                                ;and normalize by the noise. this
                                ;creates a new basis, E, which we
                                ;decompose the project the image onto


      psf = psf_struct.image
      psf = psf/total(psf)

      E = dblarr(nx[0] * nx[1], m)
      phi = reform(GH[*,0], nx[0], nx[1]) ;get 2-d representation of first shapelet
      E[0,0] = (convolve(phi, psf, ft_psf = psf_fft))[*] ;do convolutions with FFT, faster
      for j = 1, m - 1 do begin
          phi = reform(GH[*,j], nx[0], nx[1])
          E[0,j] = (convolve(phi, psf, ft_psf = psf_fft))[*] ;zero-index here for speed,
                                ;this is really a vector
      endfor
      
  endif else E = GH             ;no input PSF, just use regular shapelet basis

  Estar = E / self->cmreplicate( (sqrt( image_var ))[*], m ) 
  Etrans = transpose(Estar)

;get principal components from eigendecomposition of Gram matrix
  Gram = Estar ## Etrans        ;get Gram matrix for the shapelet basis

  eval = eigenql(Gram, eigenvect=V, /double) ;get eigendecomposition of Gram matrix
  eval = eval > 0               ;make sure eigenvalues are all positive
  d = sqrt(eval)           ;the singular values are the sqrt of the eigenvalues
                                ;find where singular values are
                                ;greater than the threshold
  keep = where( d / max(d) ge svthresh, nkeep )

  V = (transpose(V))[keep,*]    ;IDL has that funny column-major format, 
                                ;so make the PC-basis in normal matrix
                                ;format
;construct the principal component basis, in normal matrix format
;the PCs are the columns of Z, ie., PC j is z[j,*]
  Z = Etrans ## (V / self->cmreplicate(d, m)) ;note that Z is same as U in SVD
  theta = (image ## Z)          ;get the coefficients in the PC basis

;ad hoc method of optimizing the MSE between the true thetas and the
;estimated ones. first sort those thetas with d/max(d) > 0.1 in
;descending order and add them in one-at-a-time, then add in the
;remaining thetas ordered by their singular values. this will ensure
;that thetas with small sing. values are only added in at the end, as
;these are the thetas that contribute to the high variability of the
;betas.
  sortind = where(d / max(d) ge 0.1, complement=singind)
  sorted = reverse(sort(abs(theta[sortind]))) ;sort the thetas[d/max(d) > 0.1] in descending order
  sorted = [sorted, singind]    ;add in induces for the thetas[d/max(d) < 0.1]

  if auto then begin        ;choose the number of PCs to use by minimizing SURE

      risk = dblarr(m)
      risk[0] = (1.0 + total( (theta[sorted[1:*]]^2 - 1.0), /nan )) / m 

      for j = 1, m - 2 do $
        risk[j] = risk[j-1] + (1.0 - (theta[sorted[j]]^2 - 1.0)) / m

      risk[m-1] = 1.0

      minerr = min(risk, p)  ;find minimum error, p+1 is the number of PCs used
      p = p + 1 ;now p is the number of PCs used, ie, the number of parameters estimated

  endif else p = n_elements(keep)
                                ;transform back to get the shapelet coefs
  coefs = V[sorted[0:p-1],*] ## transpose(theta[sorted[0:p-1]] / d[sorted[0:p-1]])
  coefs = transpose(coefs)

  if n_elements(psf) eq 0 then begin
                                ;now find best scale on a grid of
                                ;scale values from 0.5*scale to
                                ;1.5*scale
      nscales = 10
      gammas = scale / 2.0 + dindgen(nscales) / (nscales - 1.0) * scale
      
      error = fltarr(nscales)
      
      for k = 0, nscales - 1 do begin
                                ;transform the shapelet coefficients
                                ;to have scale gammas[k]
          coefs_k = self->scale_transform( coefs, scale, gammas[k], n1, n2 )
                                ;construct new basis with this scale
          GH_k = self->construct_gh_matrix(nx, x0, gammas[k], nmax)
          GH_k = GH_k / self->cmreplicate( (sqrt( image_var ))[*], m )   
          
          ghfit = coefs_k ## GH_k
          
          C_p = self->cp_mse( image, ghfit, p )
          error[k] = C_p
          
      endfor
                                ;now use a spline to interpolate MSE
                                ;to get best scale
      nscales2 = 200
      gammas2 = scale / 2.0 + dindgen(nscales2) / (nscales2 - 1.0) * scale
      error2 = interpol(error, gammas, gammas2, /spline)
      minerr = min(error2, minind) ;find minimum error
                                ;get shapelet coefficients for this scale
      coefs_k = self->scale_transform( coefs, scale, gammas2[minind], n1, n2 )
      GH_k = self->construct_gh_matrix(nx, x0, gammas2[minind], nmax)
      
                                ;find MSE in model with best scale
      ghfit = coefs_k ## (GH_k / self->cmreplicate( (sqrt( image_var ))[*], m ))
      coefs = coefs_k

      scale = gammas2[minind]
      E = GH_k
      
  endif else ghfit = coefs ## Estar

;use Mallow's C_p statistic to estimate the mean squared error (MSE)
  C_p = self->cp_mse( image, ghfit, p )

  ghfit = coefs ## E ;get reconstructed image, not normalized to noise or corrected for PSF

  ghfit = reform(ghfit, nx[0], nx[1])

  image = reform(image, nx[0], nx[1]) * sqrt(image_var) ;restore input

  decomp = {coeff:coefs, scale:scale, n1:n1, n2:n2, C_p:C_p, npc:p}

  return, decomp
end



;Function to compute the shapelet coefficients by minimizing
;chi-squared of the true image and the reconstructed one. This is done
;by finding the p < m nested subset of the shapelet basis
;that minimize Mallow's C_p. This is the most accurate decomposition
;technique but is also the slowest.
function shapelet::nss_decomp, image_struct, nmax, $
                 image_var=image_var, $
                 ghfit=ghfit, ghmatrix=gh, $
                 status=status

  status=1
  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: '
      print,'  shapelet->nss_decomp(image_struct, nmax, '
      print,'                       image_var=, '
      print,'                       ghfit=, ghmatrix=, status=)'
      print
      message,'Halting'
  ENDIF 

  x0    = image_struct.center
  scale = image_struct.scale
  image = image_struct.image

  nimvar = n_elements(image_var)
  IF nimvar NE 0 THEN BEGIN
      IF nimvar NE 1 AND nimvar NE n_elements(image) THEN BEGIN 
          message,'image_var must be a scalar or same size as image'
      ENDIF 
  ENDIF ELSE BEGIN 
      image_var = image_struct.skysig^2
  ENDELSE 

  nx = size(image, /dim)
                                ;construct the shapelet
                                ;(Gauss-Hermite) basis
  GH = self->construct_gh_matrix(nx, x0, scale, nmax, n1=n1, n2=n2)

  m = n_elements(n1)            ;number of coefficients to estimate

  image = (image / sqrt( image_var ))[*] ;normalize by the noise and make into a vector
  GH = GH / self->cmreplicate( (sqrt( image_var ))[*], m )

;get best nested subset via forward selection
  p = 0L
  ghfit = 0d
  set_IA = lindgen(m)           ;the inactive set
  minerr = 1d30
  coefs = dblarr(m)
  error = dblarr(m)

  message,'Finding best nested subset',/inf
  repeat begin
;      print, float(p) / m
      if p gt 0 then self->split_array, lindgen(m), set_A, set_IA ;get new inactive set
      if p lt m - 1 then begin
          corr = GH[*,set_IA] ## transpose(image - ghfit) ;correlation of the shapelet states
                                ;with the image and the residual of
                                ;the image and the current best fit
          cmax = max(abs(corr), jmax) ;find maximum correlation among the inactive set
          jmax = set_IA[jmax]   ;convert jmax to index in the full set
      endif else jmax = set_IA[0]
      if p eq 0 then begin
          set_A = [jmax]        ;initiate the active set
          Gram_inv = 1.0 / ( GH[*,set_A] ## transpose(GH[*,set_A]) )
      endif else begin
          set_A = [set_A, jmax] ;add state with max corr to active set
          Gram_inv = invert( GH[*,set_A] ## transpose(GH[*,set_A]), /double )
      endelse
      p = p + 1
                                ;estimate the coefficients for this subset
      coefs_k = transpose(Gram_inv ## GH[*,set_A] ## transpose(image))
      
      ghfit = coefs_k ## GH[*,set_A] ;image reconstructed from shapelet coefficients

      error[p-1] = self->cp_mse( image, ghfit, p ) ;estimate of the MSE

      if error[p-1] lt minerr then begin
          coefs[set_A] = coefs_k
          minerr = error[p-1]
          nss = set_A           ;the best nested subset
      endif

  endrep until p eq m         ;stop when full least squares solution is reached

  if n_elements(nss) ne m then $
    self->split_array, lindgen(m), nss, set_IA ;get the inactive set for the best NSS

;now find best scale on a grid of scale values from 0.5*scale to
;1.5*scale
  nscales = 20
  gammas = scale / 2.0 + dindgen(nscales) / (nscales - 1.0) * scale

  m = n_elements(nss)  ;the number of shapelet states in the best nested subset
  error = fltarr(nscales)

  message,'Finding best scale',/inf
  for k = 0, nscales - 1 do begin
                                ;transform the shapelet coefficients
                                ;to have scale gammas[k]
      coefs_k = self->scale_transform( coefs, scale, gammas[k], n1, n2 )
                                ;construct new basis with this scale
      GH_k = self->construct_gh_matrix(nx, x0, gammas[k], nmax)
;    GH_k = GH_k / self->cmreplicate( (sqrt( image_var ))[*], m )

      ghfit = coefs_k[nss] ## GH_k[*,nss] ;only use those coefs in the best NSS
      ghfit = ghfit / (sqrt( image_var ))[*]

      C_p = self->cp_mse( image, ghfit, m )
      error[k] = C_p

  endfor

;now use a spline to interpolate MSE to get best scale
  nscales2 = 200
  gammas2 = scale / 2.0 + dindgen(nscales2) / (nscales2 - 1.0) * scale
  error2 = interpol(error, gammas, gammas2, /spline)
  minerr = min(error2, minind)  ;find minimum error
;get shapelet coefficients for this scale
  coefs_k = self->scale_transform( coefs, scale, gammas2[minind], n1, n2 )
  coefs_k[set_IA] = 0d          ;set to zero those coefs not in the active set
  GH_k = self->construct_gh_matrix(nx, x0, gammas2[minind], nmax)

;find MSE in model with best scale
  ghfit = coefs_k[nss] ## GH_k[*,nss]
  ghfit = ghfit / (sqrt( image_var ))[*]

;use Mallow's C_p statistic to estimate the mean squared error (MSE)
  C_p = self->cp_mse( image, ghfit, m )

  ghfit = ghfit * (sqrt( image_var ))[*]
  ghfit = reform(ghfit, nx[0], nx[1])

; note needed since image is a copy anyway
;  image = reform(image, nx[0], nx[1]) * sqrt(image_var) ;restore input
  decomp = {coeff:coefs_k, scale:gammas2[minind], n1:n1, n2:n2, C_p:minerr}

  status=0
  return, decomp
end







;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; ROUTINE TO RECONSTRUCT AN IMAGE, GIVEN ITS SHAPELET COEFFICIENTS AND
; SCALE. 
;
; AUTHOR : BRANDON C. KELLY, STEWARD OBS., AUG 2005
;
; INPUTS :
;
;   COEF - THE SHAPELET COEFFICIENTS.
;   SCALE - THE SHAPELET SCALE.
;   NX - THE IMAGE SIZE, A 2-ELEMENT VECTOR. THE RESULTING IMAGE WILL
;        BE AN [NX[0], NX[1]] ARRAY.
;   NMAX - THE MAXIMUM SHAPELET ORDER, N1 + N2 <= NMAX
; 
; OPTIONAL INPUTS :
;
;   X0 - THE SHAPELET CENTER. THE DEFAULT IS THE CENTER OF THE IMAGE.
;   N1 - THE HORIZONTAL ORDER OF THE SHAPELETS. IF NOT INPUT, THEN THE
;        ORDERING PERFORMED AUTOMATICALLY BY CONSTRUCT_GH_MATRIX.PRO
;        WILL BE ASSUMED. IF THE SHAPELET INFORMATION IS THAT OUTPUT
;        FROM BRANDON KELLY'S SHAPELET_DECOMP.PRO, THEN DO NOT INPUT
;        N1 OR N2. NOTE THAT COEF[J] WILL CORRESPOND TO THE SHAPELET
;        OF ORDER |N1[J], N2[J]>.
;   N2 - SAME AS N1, BUT FOR THE VERTICAL ORDER.
;
; OUTPUT :
;
;   THE IMAGE RECONSTRUCTED FROM ITS SHAPELET COEFFICIENTS.
;
; ROUTINES CALLED :
;   CONSTRUCT_GH_MATRIX
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION shapelet::recon, coef, scale, imsize, nmax, center=center, $
                 n1=n1, n2=n2, status=status

  status=1
  IF n_params() LT 4 THEN BEGIN
      print, 'Syntax- result = shapelet_recon( coef, scale, imsize, nmax, center=center, n1=n1, n2=n2, status=)'
      return, 0
  ENDIF
  
  IF n_elements(center) EQ 0 THEN center = imsize / 2
  
  ;; construct the shapelet (Gauss-Hermite) basis
  GH = self->construct_gh_matrix(imsize, center, scale, nmax, n1=n1, n2=n2, $
                                 status=ghstatus)
  
  IF ghstatus NE 0 THEN BEGIN 
      return, -1
  ENDIF 
  ;; reconstruct the image
  image = coef ## GH            
  ;; the image is a vector, make into an array
  image = reform(image, imsize[0], imsize[1]) 
  
  status = 0
  return, image
END














FUNCTION shapelet::poisson_var, image, skysig
  muhat = $
    0.5 * sqrt( 4.0 * self->positive( image^2 - skysig^2 ) + 1.0 ) - 0.5
  return,muhat
END 

FUNCTION shapelet::decomp_struct, nmax, dtype, nmax_psf=nmax_psf

  ncoeff = self->nmax2ncoeff(nmax)
  decomp = { dtype:dtype, $
             scale:-1d,                                $
             nmax:nmax,                                $
             ncoeff: ncoeff,                           $
             coeff:dblarr(ncoeff),                      $
             n1:intarr(ncoeff),                        $
             n2:intarr(ncoeff),                        $
             flux:0d,                                  $
             $
             center: [0d, 0d],                         $ 
             size: [0,0],                              $
             $
             C_p:9999d                                 $
           }


  IF n_elements(nmax_psf) NE 0 THEN BEGIN 
      psf_ncoeff = self->nmax2ncoeff(nmax_psf)
      decomp = $
        create_struct(decomp, $
                      'psf_nmax',  nmax_psf,           $
                      'psf_ncoeff', psf_ncoeff,        $
                      'psf_scale', -1d,                $
                      'psf_coeff',  dblarr(psf_ncoeff), $
                      'psf_n1',    intarr(psf_ncoeff), $
                      'psf_n2',    intarr(psf_ncoeff) )
  ENDIF 

  IF dtype EQ 'approx' THEN BEGIN 
      decomp = create_struct(decomp,'ortho',-1d)
  ENDIF 
  IF dtype EQ 'pcr' THEN BEGIN 
      decomp = create_struct(decomp,'npc',0)
  ENDIF

  return,decomp

END 


FUNCTION shapelet::make_gauss, sigma

  im_size = 2*4*[sigma, sigma]
  makegauss, gauss, im_size, sigma
  return, gauss

END 
FUNCTION shapelet::make_exp, rzero

  im_size = 2*4*[rzero, rzero]
  make_exp, exp, im_size, rzero
  return, exp

END 

FUNCTION shapelet::make_image, type, scale, counts

  CASE type OF
      'gauss': image=self->make_gauss(scale)
      'exp': image = self->make_exp(scale)
  ENDCASE 

  image = image/total(image)*counts
  return, image

END 

FUNCTION shapelet::addnoise, input_image, s2n, noise=noise

  image = input_image
  ;; s2n = sqrt( sum of s2n^2 over pixels )
  ;; assume sky noise dominated, so it's just 
  ;; a fixed number.  
  ;; s2n^2 = sum( (im_i/noise_i)^2 )
  ;;     = 1/noise^2*sum( im_i^2 )
  ;; -> noise = sqrt( total(im^2) )/s2n

  timage2 = total(image^2)

  noise = sqrt(timage2)/s2n

  npix = n_elements(image)
  image[*] = image[*] + noise*randomu(seed, npix, /normal)

;  print,'s2n = ',s2n
;  print,'total(image) = ',total(image)
;  print,'noise = ',noise
;  print,'calculated s2n = ',sqrt( total( image^2/noise^2 ) )
;  print,'From original = ',sqrt( timage2/noise^2 )

  return, image

END 

PRO shapelet::plot_crossection, image, _extra=_extra
  im_size = size(image, /dim)
  cen = (im_size-1.0)/2.0

  image_cross = image[cen[0], *]
  image_cross = image_cross*im_size[1]/max(image_cross)*0.75
  pplot, lindgen(im_size[0]), image_cross, /overplot, _extra=_extra
END 
PRO shapelet::plot_image, image, title=title
  tvasinh, image, title=title
  self->plot_crossection, image
END 
PRO shapelet::test_decomp, type, galtype, correct_psf=correct_psf

  nmax = 10
  nmax_psf=10

  psf_sigma = 2.0
  psf = self->make_gauss(psf_sigma)
  psf_size = size(psf, /dim)
  psf_cen = (psf_size-1.0)/2.0


  counts = 100.0
  CASE galtype OF
      'gauss': BEGIN 
          imscale = 10.0
          o_image = self->make_image('gauss', imscale, counts)
      END 
      'exp': BEGIN 
          imscale = 7.0
          o_image = self->make_image('exp', imscale, counts)
      END 
  ENDCASE 
  im_size = size(o_image, /dim)
  cen = (im_size - 1.0)/2.0

  ;; PSF convolved galaxy
;  c_image = convol(o_image, psf, edge_truncate=1)
  oc_image = convolve(o_image, psf)

  !p.multi = [0,2,3]
  !p.charsize = 2


  ;; Plot the original galaxy
  self->plot_image, o_image,  title='Original '+galtype+' galaxy'
  self->plot_image, oc_image, title='Convolved image'

  ;; add noise
  s2n = 100.0
  image = self->addnoise(o_image, s2n, noise=old_noise)
  c_image = self->addnoise(oc_image, s2n, noise=noise)

  self->plot_image, image,   title='Galaxy+noise'
  self->plot_image, c_image, title='Convolved+noise'
  self->plot_image, psf, title='psf'


  ;; Now the decomposition

  isize = sqrt( imscale^2 + psf_sigma^2 )
  imstruct = self->image_struct(c_image, cen, noise, isize)
  psfstruct= self->image_struct(psf, psf_cen, 1.0, psf_sigma)

  tm = systime(1)
  decomp = self->decomp(imstruct, type, nmax, $
                        correct_psf=correct_psf, $
                        psf_struct=psfstruct, nmax_psf=nmax_psf, $
                        recon=recon)
  ptime, systime(1)-tm


  title = 'Reconstruction'
  IF keyword_set(correct_psf) AND type EQ 'pcr' THEN BEGIN 
      title=title + ' (deconvolved)'
      leg = 'Original'
  ENDIF ELSE BEGIN 
      leg = 'Convolved'
  ENDELSE 
  self->plot_image, recon, title=title

  IF keyword_set(correct_psf) AND type EQ 'pcr' THEN BEGIN 
      self->plot_crossection, o_image,color=!red
      legend,$
        ['Recon','Original'], $
        lines=[0,0],colors=[!p.color, !red], $
        /right, box=0, charsize=1

  ENDIF ELSE BEGIN 
      self->plot_crossection, o_image, color=!darkGreen,line=2
      self->plot_crossection, oc_image,color=!red
      legend,$
        ['Recon','Convolved','Original'], $
        lines=[0,0,2],colors=[!p.color, !red,!darkGreen], $
        /right, box=0, charsize=1

  ENDELSE 


  !p.multi=0

  

END 




; Makes a copy of the image
FUNCTION shapelet::image_struct, image, center, skysig, scale
  
  sz = size(image, /dim)

  return,$
    {size:   sz,$
     center: center, $
     skysig: skysig, $
     scale:  scale,  $
     image:  image}

  
END 


FUNCTION shapelet::decomp, image_struct, dtype_in, nmax, $
                 psf_struct=psf_struct, $
                 nmax_psf=nmax_psf, $
                 correct_psf=correct_psf, $
                 $
                 max_size=max_size, $ ; max size for approx method
                 $
                 svthresh=svthresh, $
                 $              ;Optional outputs
                 recon=recon, $
                 ghmatrix=ghmatrix, $
                 status=status

  on_error, 2
  status = 1
  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: '
      print,'  shapelet->decomp(image_struct, dtype, nmax, '
      print,'                   psf_struct=, nmax_psf=,'
      print,'                   /correct_psf, '
      print,'                   max_size=, ; max size for approx method'
      print,'                   svthresh=, ; svd threshold for pcr'
      print,'                   recon=,    ; reconstructed image'
      print,'                   ghmatrix=,'
      print,'                   status='
      print
      message,'Halting'
  ENDIF 

  dtype = strlowcase(strtrim(dtype_in))
  if dtype NE 'nss' AND dtype ne 'ols' and dtype ne 'approx' and dtype ne 'pcr' then begin
      print, 'If input, DTYPE must be one of the following : '
      print, '   APPROX, PCR'
      message,'Halting'
  ENDIF
  IF n_elements(svthresh) GT 0 THEN BEGIN 
      IF (svthresh GT 1 OR svthresh LT 0) THEN BEGIN 
          print, 'SVTHRESH must be less than or equal to 1 ' + $
            'and greater or equal to 0.'
          message,'Halting'
      ENDIF 
  ENDIF


  npsf = n_elements(psf_struct)
  IF npsf NE 0 THEN BEGIN 
      wpsf = where(psf_struct.image EQ 0, nwpsf)

      IF nwpsf NE 0 THEN BEGIN 
          psf_struct.image[wpsf] = $
            psf_struct.image[wpsf] + $
            image_struct.skysig * randomn(seed, nwpsf)
      ENDIF 

  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; create the output struct
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  decomp = self->decomp_struct(nmax, dtype, nmax_psf=nmax_psf)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; estimate of the poisson variance
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; Same size as image
  muhat = self->poisson_var(image_struct.image, image_struct.skysig)  
  ;;estimate of the variance in the image
  image_var = muhat + image_struct.skysig^2  

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; decompose image based on method. Only support approx and pcr for now
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  case dtype of
      'approx' : BEGIN 
          decomp0=$
            self->approx_decomp(image_struct, image_var, nmax, $
                                max_size=max_size, $
                                ghfit=recon, $
                                ghmatrix=ghmatrix, $
                                status=dstatus)
          IF dstatus NE 0 THEN BEGIN 
              delvarx, recon, ghmatrix
              return, -1
          ENDIF 
          IF npsf NE 0 THEN BEGIN 
              dpsf = self->pcr_decomp(psf_struct, nmax_psf)
          ENDIF 
      END 
      'nss': BEGIN 
          decomp0=$
            self->nss_decomp(image_struct, nmax, $
                             image_var=image_var, $
                             ghfit=recon, $
                             ghmatrix=ghmatrix, $
                             status=dstatus)
          IF dstatus NE 0 THEN BEGIN 
              delvarx, recon, ghmatrix
              return, -1
          ENDIF 
          IF npsf NE 0 THEN BEGIN 
              dpsf = self->pcr_decomp(psf_struct, nmax_psf)
          ENDIF 
      END 
      'pcr' : begin
          IF npsf GT 0 THEN BEGIN               
              dpsf=self->pcr_decomp(psf_struct, nmax_psf)
          ENDIF 
          
          ;; An intrinsic idl proc stupidly crashes instead of
          ;; exiting with an error status.  We can catch the error
          ;; by using execute
          
          command = $
            'decomp0 = '+$
            '  self->pcr_decomp(image_struct, nmax, '+$
            '   image_var=image_var, ghfit=recon, ghmatrix=ghmatrix,'+$
            '   svthresh=svthresh'
          IF keyword_set(correct_psf) THEN $
            command = command+', psf_struct=psf_struct'
          
          command = command + ')'
          IF NOT execute(command) THEN BEGIN 
              delvarx, recon, ghmatrix
              return, -1
          ENDIF 
          
      END 
      ELSE: message,'Unknown dtype: '+string(dtype)
  ENDCASE 


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; store shapelet decomposition data into output structure
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  decomp.scale = decomp0.scale
  decomp.coeff = decomp0.coeff
  decomp.n1 = decomp0.n1
  decomp.n2 = decomp0.n2
  decomp.C_p = decomp0.C_p
  if npsf NE 0 then begin
      decomp.psf_scale = dpsf.scale
      decomp.psf_coeff = dpsf.coeff
      decomp.psf_n1 = dpsf.n1
      decomp.psf_n2 = dpsf.n2
      
      psf_flux = self->flux( dpsf.coeff, dpsf.scale, $
                             dpsf.n1, dpsf.n2 )
      
      ;;normalize PSF to have flux of unity
      decomp.psf_coeff = decomp.psf_coeff / psf_flux 
  endif
  
  ;; normalize coefficients to have flux
  ;; of unity (in counts)
  flux = self->flux( decomp0.coeff, decomp0.scale, $
                     decomp0.n1, decomp0.n2 )
  decomp.coeff = decomp.coeff / flux
  decomp.flux = flux
  
  decomp.center = image_struct.center
  decomp.size = image_struct.size
  
  IF dtype EQ 'approx' THEN BEGIN 
      ortho = 1.0d - total( total(GHMatrix^2, 1) ) / decomp.ncoeff
      decomp.ortho = ortho 
  ENDIF 

  IF dtype EQ 'pcr' THEN BEGIN 
      decomp.npc = decomp0.npc
  ENDIF 



  status = 0
  return,decomp

END 
























FUNCTION shapelet::atlas_decomp, gal, dtype, nmax, clr,$
                 nmax_psf=nmax_psf, $
                 correct_psf=correct_psf, $
                 $
                 max_size=max_size, $
                 $
                 svthresh=svthresh, $
                 $
                 trim_atlas=trim_atlas, $
                 trim_nsigma=trim_nsigma, $
                 trimmed_image=trimmed_image, $
                 slack=slack, $
                 $
                 $              ; optional outputs
                 psFieldInfo=psFieldInfo, $
                 recon=recon, $
                 image_struct=image_struct, $
                 ghmatrix=ghmatrix, $
                 status=status

;  on_error,2
  status = 1
  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: '
      print,'  shapelet->atlas_decomp(galstruct, dtype, nmax, clr, '
      print,'                         nmax_psf=,'
      print,'                         /correct_psf, '
      print,'                         max_size=, ; max size for approx method'
      print,'                         svthresh=, ; svd threshold for pcr'
      print,'                         /trim_atlas, trim_nsigma=, '
      print,'                         trimmed_image=, '
      print,'                         image_struct=, '
      print,'                         recon=,    ; reconstructed image'
      print,'                         ghmatrix=,'
      print,'                         status='
      print
      message,'Halting'
  ENDIF 

  delvarx, recon, ghmatrix, psFieldInfo

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Check the basic arguments
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(gal) GT 1 THEN message,'One object at a time'
  nclr = n_elements(clr)
  IF n_elements(clr) GT 1 THEN message,'One clr at a time'
  IF clr LT 0 OR clr GT 4 THEN message,'clr must be win [0,4]'
  band = !colors[clr]



  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Check optional args
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  ;; Should we trim the atlas images?
  IF n_elements(trim_nsigma) NE 0 THEN BEGIN 
      IF NOT keyword_set(trim_atlas) THEN trim_atlas=1
  ENDIF 

  IF n_elements(slack) EQ 0 THEN slack = 2

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Get the atlas image
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  CASE clr OF
      0: read_atlas, gal, $
        row0=row0, col0=col0, imu=image, status=astatus, /silent
      1: read_atlas, gal, $
        row0=row0, col0=col0, img=image, status=astatus, /silent
      2: read_atlas, gal, $
        row0=row0, col0=col0, imr=image, status=astatus, /silent
      3: read_atlas, gal, $
        row0=row0, col0=col0, imi=image, status=astatus, /silent
      4: read_atlas, gal, $
        row0=row0, col0=col0, imz=image, status=astatus, /silent
  ENDCASE 




  ;; Was the atlas image found?
  IF astatus EQ 0 THEN BEGIN 
              
;      atlas_dir, gal.run, gal.rerun, gal.camcol, atldir
              
      ;; get the PSF data
              
;      psfield_name, gal.run, gal.camcol, gal.field, name
;      psfield = atldir + name
;      psfield = psfield[0]
;      read_psf, atldir + name, psp ;get psf data

      psp = !sdss->psfield_read(gal.run, gal.camcol, gal.field)
      psf = !sdss->psfrec(*psp[clr], gal.rowc[clr], gal.colc[clr])

;      psf_rec, psp, gal.rowc, gal.colc, psfp, clr ;get psf image
;      psf = *psfp[0]

      skysig = (*psp[5]).skysig[clr] ;get the sky variance
      
                                ;get the galaxy center
      rowctr = gal.rowc[clr] - row0[clr]
      colctr = gal.colc[clr] - col0[clr]
      
      ;; Trim image of zeros if necessary.  refind the center
      IF keyword_set(trim_atlas) THEN BEGIN 
          
          nskysig = 2
          trimmed_image = $
            self->trim_image(image, 0.0, skysig, nskysig, $
                             smoothing_window=4, slack=slack, $
                             minx=minx, miny=miny, status=trstatus)
          IF trstatus NE 0 THEN BEGIN 
              delvarx, recon, ghmatrix
              return, -1
          ENDIF 
          
          image=0
          IF arg_present(trimmed_image) THEN BEGIN 
              image = trimmed_image
          ENDIF ELSE BEGIN 
              image = temporary(trimmed_image)
          ENDELSE 
          colctr = colctr - minx
          rowctr = rowctr - miny
          
      ENDIF 
      
      ;; Center in atlas image
      x0 = [colctr, rowctr]
      
      ;; assume PSF center is center of PSF image
      psfx0 = ( size(psf, /dim) -1.0 ) / 2.0
      
      ;; initial guess for the galaxy scale
      ;; based adaptively-weighted moments
      scale = sqrt( gal.m_rr_cc[clr] / 2.0 )
      psf_scale = sqrt( gal.m_rr_cc_psf[clr] / 2.0 )

      ;; User can return some psf info. E.S.S. 
      ;; Whoops, only returned for last object
      IF arg_present(psFieldInfo) THEN BEGIN 
          psFieldInfo = {                  $
                          psp1:*psp[0],    $
                          psp2:*psp[1],    $
                          psp3:*psp[2],    $
                          psp4:*psp[3],    $
                          psp5:*psp[4],    $
                          psp6:*psp[5],    $
                          psf:psf,         $
                          skysig:skysig    $
                        }
      ENDIF 
      
      ptr_free, psp

      image_struct = self->image_struct(image, x0, skysig, scale)
      psf_struct = self->image_struct(psf, psfx0, 1.0, psf_scale)


  ENDIF ELSE BEGIN 
      message,'Atlas image not found',/inf
      return,-1
  END 

  IF n_elements(nmax_psf) EQ 0 THEN nmax_psf = 10
  decomp = self->decomp(image_struct, dtype, nmax, $
                        psf_struct=psf_struct, nmax_psf=nmax_psf, $
                        correct_psf=correct_psf, $
                        max_size=max_size, svthresh=svthresh, $
                        recon=recon, ghmatrix=ghmatrix, status=dstatus)


  IF dstatus NE 0 THEN return,-1




  ;; Add in extra tags for atlas images. This is info that was
  ;; already required to run this program, so we can just send 
  ;; the output struct back into this program for further tests

  addstruct = { $
                run:         gal.run, $
                rerun:       gal.rerun, $
                camcol:      gal.camcol, $
                field:       gal.field, $
                id:          gal.id, $
                rowc:        gal.rowc, $
                colc:        gal.colc, $
                m_rr_cc:     gal.m_rr_cc, $
                m_rr_cc_psf: gal.m_rr_cc_psf, $
                clr:         clr $
              }
  
  decomp = create_struct(addstruct,decomp)
  status = 0
  return,decomp


END 










FUNCTION shapelet::cleanup
  return,1
END 

PRO shapelet__define

  struct = {$
             shapelet, $
             dummy: 0 $
           }

END 








