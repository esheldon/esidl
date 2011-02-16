;+
; NAME:
;       QGAUSS2D
;
; PURPOSE:
;       This function computes the double integral of a bivariate
;       function F(x,y) using either a dy-dx or dx-dy order of integration.
;       Modified from INT_2d to do gauss-legendre integration and to
;       also take array input rather than just functions.  
;
; CATEGORY:
;       Numerical Analysis.
;
; CALLING SEQUENCE:
;    Method 1:
;       Result = QGAUSS2D(z, x, y, npts, /silent)
;    Method 2:
;       Result = QGAUSS2D(Fxy, AB_Limits, PQ_Limits, Pts, order=, /silent)
;
; INPUTS:
;   Method 1:
;       z, x, y: the surface z, and the abscissa x,y over which to
;           integrate.
;
;   Method 2:
;       Fxy:  A scalar string specifying the name of a user-supplied
;             IDL function that defines the bivariate function to be
;             integrated. The function must accept x & y and return
;             a scalar result.
;       AB_Limits:  A two-element vector containing the lower, A, and
;                   upper, B, limits of integration. These limits correspond
;                   to x when using a dy-dx order of integration or they
;                   correspond to y when using a dx-dy order of integration. 
;       PQ_Limits:  A scalar string specifying the name of a user-supplied 
;                   IDL function that defines the lower, P, and upper, Q, 
;                   limits of integration. These limits correspond to y when
;                   using a dy-dx order of integration or they correspond to
;                   x when using a dx-dy order of integration. They must 
;                   accept x (for a dy-dx order of integration) or y (for a
;                   dx-dy order of integration) and return a two-element
;                   vector result.
;       Pts:  The number of transformation points used in the
;             computation. 
;
; KEYWORD PARAMETERS:
;       ORDER:  A scalar value of either 0 or 1. If set to 0, the integral 
;               is computed using a dy-dx order of integration. If set to 1,
;               the integral is computed using a dx-dy order of integration.
;               This is not relavent for Method 1.
;
; EXAMPLE:
;   Method1: 
;       Compute integral of a surface on an x-y grid.
;       res = qgauss2d( z, x, y, 100 )
;
;   Method2:
;
;       Compute the double integral of the bivariate function
;       F(x, y) = y * cos(x^5) over the region: [0.0, x^2] for y and 
;                                               [0.0, 2.0] for x using 
;                                               a dy-dx order of integration.
;
;       ;Define the limits of integration for y as a function of x.
;         function PQ_Limits, x
;           return, [0.0, x^2]
;         end
;
;       ;Define the limits of integration for x.
;         AB_Limits = [0.0, 2.0]
;
;       ;Integrate with 48 and 96 point formulas using a dy-dx order of 
;       ;integration
;         print, QGAUSS2D('Fxy', AB_Limits, 'PQ_Limits', 48) 
;         print, QGAUSS2D('Fxy', AB_Limits, 'PQ_Limits', 96)
;
;       QGAUSS2D with 48 transformation points yields:    0.055142668
;       QGAUSS2D with 96 transformation points yields:    0.055142668
;
;       Compute the double integral of the bivariate function
;       F(x, y) = y * cos(x^5) over the equivalent region using a dx-dy order
;       of integration.
;
;       ;Define the limits of integration for x as a function of y.
;         function PQ_Limits, y
;           return, [sqrt(y), 2.0]
;         end
;
;       ;Define the limits of integration for y.
;         AB_Limits = [0.0, 4.0]
;
;       ;Integrate with 48 and 96 point formulas using a dx-dy order of
;       ;integration 
;         print, QGAUSS2D('Fxy', AB_Limits, 'PQ_Limits', 48, /order)
;         print, QGAUSS2D('Fxy', AB_Limits, 'PQ_Limits', 96, /order)
;
;       QGAUSS2D with 48 transformation points yields:    0.055142678
;       QGAUSS2D with 96 transformation points yields:    0.055142668
;
;       The exact solution (9 decimal accuracy) yields: 0.055142668
;
; MODIFICATION HISTORY:
;       Written by:  GGS, RSI, January 1994
;       Modified:    GGS, RSI, September 1994
;                    Added 96 point transformation data.               
;                    Added DOUBLE keyword. Replaced nested FOR loop with
;                    vector operations resulting in faster execution.
;       Modified:    GGS, RSI, June 1995
;                    Added ORDER keyword allowing a dx-dy order of integration.
;       Modified:    GGS, RSI, April 1996
;                    Modified keyword checking and use of double precision.
;       Modified:    Erin Scott Sheldon UofMichigan  Now calculates the
;                    weights and abscissa with gauleg  
;                    This allows _any_ integer PTS
;       Modified:    Now may integrate x,y,z surfaces given as arrays in
;                    addition to functions.
;                    All calculations in double precision
;                    Erin Sheldon NYU, July 2006
;-          



; some test functions
FUNCTION _gauss2dpqlim, x
  return, [-5.0, 5.0]
END 
FUNCTION _gauss2dfunc, x, y

  axis_ratio = 0.8
  angle = !pi/5.0
  c = cos(angle)
  s = sin(angle)
  xp =  x*c + y*s
  yp = (-x*s + y*c)/axis_ratio

  return, exp( -0.5*( xp^2 + yp^2) )

END 
FUNCTION _gauss2d_data
  n=100
  x = arrscl(findgen(n), -5.0, 5.0 )
  y = arrscl(findgen(n), -5.0, 5.0 )

  z = fltarr(n,n)
  FOR i=0L,n-1 DO BEGIN 
      FOR j=0L,n-1 DO BEGIN 
          z[i,j] = _gauss2dfunc(x[i],y[j])
      ENDFOR 
  ENDFOR 

  struct = {x:x, y:y, z:z}
  return, struct
END 









PRO _qgauss2d_gauleg, npts, silent=silent
  COMMON gauss2d_block, XXi, WWi

  ;; only need to do this if npts changes
  IF npts NE n_elements(WWi) THEN BEGIN 
      IF NOT keyword_set(silent) THEN BEGIN 
          print,'QGAUSS2D: Calculating weights and abscissa for npoints =  ',npts
      ENDIF 
      gauleg, -1., 1., npts, XXi, WWi
  ENDIF 

END 


FUNCTION _qgauss2d_interp_data, z, x, y
  COMMON gauss2d_block, XXi, WWi


  npts = n_elements(XXi)

  st = $
    { $
      zinterp: dblarr(npts,npts), $
      xf1: x[0], $
      xf2: x[0], $
      yf1: x[0], $
      yf2: x[0] $
    }

  ;; x and y must be ordered
  x1 = x[0]
  x2 = x[n_elements(x)-1]
  y1 = y[0]
  y2 = y[n_elements(y)-1]
  
  st.xf1 = (x2-x1)/2.0
  st.xf2 = (x2+x1)/2.0

  st.yf1 = (y2-y1)/2.0
  st.yf2 = (y2+y1)/2.0

  xvals = XXi*st.xf1 + st.xf2
  yvals = XXi*st.yf1 + st.yf2
  
  nx = n_elements(x)
  ny = n_elements(y)

  xi = lindgen(nx)
  yi = lindgen(ny)

  xi_interp = interpol( xi, x, xvals )
  yi_interp = interpol( yi, y, yvals )

  FOR i=0L,npts-1 DO BEGIN 
      FOR j=0L, npts-1 DO BEGIN 
          st.zinterp[i,j] = interpolate(z, xi_interp[i], yi_interp[j])
      ENDFOR 
  ENDFOR 

  return, st
END 


FUNCTION _qgauss2d_data, z, x, y, npts, silent=silent
  
  COMMON gauss2d_block, XXi, WWi

  _qgauss2d_gauleg, npts, silent=silent

  st = _qgauss2d_interp_data(z, x, y)

  Aj = 0.0

  FOR i=0L, npts-1 DO BEGIN 

      zidata = reform( st.zinterp[i, *] )
      Aj = Aj + WWi[i]*total( zidata*WWi, /double)*st.yf1
      
  ENDFOR 

  return, Aj*st.xf1
  

END 



FUNCTION _qgauss2d_func, Fxy, AB_Limits, PQ_Limits, Npts, $
                         order=order, silent=silent
  
  ON_ERROR, 2
  
  if N_ELEMENTS(AB_Limits) ne 2 then $
    MESSAGE, "AB_Limits parameter must be a two-element vector."
  
  
  ;; Rewrote: calculates weights and abscissa, but does repeat
  ;; unless called again with different value of Npts. E.S.S
  
  COMMON gauss2d_block, XXi, WWi

  _qgauss2d_gauleg, npts, silent=silent

  
  nLimit = [AB_Limits[1] - AB_Limits[0], AB_Limits[1] + AB_Limits[0]] / 2.0
  Aj = 0d
  
  if KEYWORD_SET(Order) eq 0 then begin
      ;; Compute using a dy-dx order of integration.
      for i = 0, Npts-1 do begin
          X  = nLimit[0] * XXi[i] + nLimit[1]
          fLimit = CALL_FUNCTION(PQ_Limits, X)
          kF = (fLimit[1] - fLimit[0]) / 2.0
          kS = (fLimit[1] + fLimit[0]) / 2.0
          jX = TOTAL(WWi * CALL_FUNCTION(Fxy, X, kF*XXi+kS), /Double)
          Aj = Aj + (WWi[i] * kF * jX)
      endfor
  endif else if Order eq 1 then begin 
      ;; Compute using a dx-dy order of integration.
      for i = 0, Npts-1 do begin
          Y  = nLimit[0] * XXi[i] + nLimit[1]
          fLimit = CALL_FUNCTION(PQ_Limits, Y)
          kF = (fLimit[1] - fLimit[0]) / 2.0
          kS = (fLimit[1] + fLimit[0]) / 2.0
          jY = TOTAL(WWi * CALL_FUNCTION(Fxy, kF*XXi+kS, Y), /Double)
          Aj = Aj + (WWi[i] * kF * jY)
      endfor
  endif else MESSAGE, "Order keyword must be 0 or 1."
  
  RETURN, (Aj * nLimit[0])
  
end

FUNCTION qgauss2d, par1, par2, par3, par4, order=order, silent=silent, $
                 double=double  ; ignored

  on_error, 2
  np = n_params()
  IF np NE 3 AND np NE 4 THEN BEGIN
      print,'-Syntax: '
      print,'         result=qgauss2d(z, x, y, npts, order=, /silent)'
      print,' OR'
      print,'         result=qgauss2d(funcname, AB_Limits, PQ_Limits, npts, order=, /silent)'
      print,'            Where funcname is a string name of a function '
      print,'            taking only x,y'
      print
      message,'Halting'
  ENDIF 
  
  IF size(par1, /type) EQ 7 THEN BEGIN 
      return, _qgauss2d_func(par1, par2, par3, par4, order=order, silent=silent)
  ENDIF ELSE BEGIN 
      return, _qgauss2d_data(par1, par2, par3, par4, silent=silent)
  ENDELSE 

END 


