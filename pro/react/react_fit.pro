PRO REACT_FIT,k,P,err=err, yfit, alpha=alpha, tau=tau, radius=radius, thetas=thetas, plot=plot

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;       REACT_FIT 
;
; PURPOSE:
;       Perform a non-parametric fit to a one-dimensional function.
;       Assumes the form y = f(x) + e, where y are the inputed
;       response variables, x is the inputed covariate, and e 
;       us some unknown error (which has a Gaussian dispersion).
;       The y measurements  may (or may not) have specified errors
;       which can also be inputed. These are carried through the
;       REACT technique as weights. The method is described in
;       Miller, Nichol, Genovese, and Wasserman, (2001).
;
; CALLING SEQUENCE:
;       react_fit, x, y, [err=err, yfit, alpha=alpha, tau=tau, radius=radius, thetas=thetas, /plot]
;
; INPUTS:
;       k,P: The X and Y data values.
;
; OPTIONAL INPUTS:
;       err: The error on the Y measurements
;     alpha: The confidence level (0-1). Defaults to 0.95.
;
; KEYWORD PARAMETERS:
;       /plot: Make a plot to file react_out.ps
;
; OUTPUTS:
;       yfit: The best non parametric fit. **Note** that the original
;             X and Y values have been sorted according to increasing
;             X values. These new, non-parametrically fitted values
;             of y (yfit) are functions of the SORTED X values. 
;
; OPTIONAL OUTPUTS:
;       tau:  see the paper (used in confidence radius) 
;       radius: see the paper 
;       thetas: see paper  (useful for comparision to other functions)
;
; CALLED ROUTINES:
;
;    DIAG (included in tar file)
;    GET_DIAG (included in tar file)
;    DIFF (included in tar file)
;    SHRINK (included in tar file)
;    LINKEDLIST (got from www.dfanning.com)
;
; PROCEDURE:
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;some checks;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

IF (N_params() eq 0) THEN BEGIN
  PRINT,'-Syntax: react_fit, x, y, [err=err, yfit, alpha=alpha, tau=tau, radius=radius,thetas=thetas,/plot]'
  PRINT,'   Returns yfit, a non-parametric fit to the data'
  RETURN
ENDIF 

IF (n_elements(k) ne n_elements(P)) THEN BEGIN
  PRINT, 'x vector length not equal to y vector length'
  RETURN
ENDIF
IF (n_elements(alpha) EQ 0) THEN alpha = 0.95
IF (n_elements(err) ne 0) and (n_elements(err) ne n_elements(P)) THEN BEGIN
  PRINT, 'err vector length not equal to x vector length'
  RETURN
ENDIF
IF (n_elements(err) EQ 0) THEN BEGIN
   err = 1.
ENDIF 

;;;;;;;; some parameters  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


IF KEYWORD_SET(plot) THEN BEGIN 
 red=[0,1,1,0,0,1]
 green=[0,1,0,1,0,1]
 blue=[0,1,0,0,1,0]
 tvlct, 255*red, 255*green,255*blue
 set_plot,'ps'
 device,/inches,xsize=3.5,scale_factor=1.0
 device,/inches,ysize=3.5,scale_factor=1.0
 device,color=1
 device, file="react_out.eps",/ENCAPSULATED
 !P.CHARSIZE=1.2
 !P.CHARTHICK = 2.0
 !P.THICK = 4.0
ENDIF ELSE BEGIN
 set_plot, 'x'
ENDELSE

;First, we must sort the y values
size_array = n_elements(k)
x = double(k[sort(k)])
y = double(P[sort(k)])
IF (n_elements(err) eq n_elements(k)) THEN BEGIN
 noise_cov = err[sort(k)]^2.0
ENDIF ELSE BEGIN
 noise_cov = dblarr(n_elements(k))
 noise_cov = noise_cov + 1.0
 skip_comp = round(n_elements(k)/2)
ENDELSE

; Create the cosine basis
n = n_elements(y)
a = 1/sqrt(double(n))
U = dblarr(n,n)
U[*,0] = a
a = sqrt(2.d0)*a
b = double(!PI/double(2*n))
xx = indgen(n)
FOR I = 1, n -1 DO BEGIN
  U[*,I] = a*cos((2.0*double(xx) +1.d0)*double(I)*b)
ENDFOR
; Now, U is the basis (an nxn array)


; Estimate z from the data using our new basis
z = (U)##(double(y))

;If no error given, estimate sigmasq using the *** method
IF (n_elements(err) ne n_elements(x)) THEN BEGIN
 sigmasq = total(z[skip_comp+1:n-1]*z[skip_comp+1:n-1])/double(n-skip_comp-1)
ENDIF ELSE BEGIN
 ;Or Use the real errors (as given by the user)
 diagB = (U^2) ## noise_cov
 sigmasq = 1.0
 sigmasq = sigmasq*diagB
ENDELSE

;Now, Calculate fstar
m = n_elements(z)
g = (z^2 - sigmasq)/(z^2)
w = z^2/double(m)
;print, n_elements(w)
; don't destroy g (yet)
f = g

; Our end-goal will be to make a shrunken array of
; f's.  There will be a corresponding array (same size)
; that contains the indices to which the f's will be
; assigned.  The final fstar will be its original length.
; Therefore, each position in the index array will need
; to contain multiple indexes.  The way to do this is via
; a LinkedList (see DFANNING's LinkedList routines.)
; 
; Below I create the initial LinkedList (called index)
; in which there are m indices, each containing one
; index (itself).
index = Obj_New("LINKEDLIST", 0)
FOR I = 1, m - 1 DO BEGIN  
 index -> Add, I
ENDFOR

WHILE (n_elements(f) gt 1) DO BEGIN
 d = diff(f)
 violators = where(d gt 0)
 if (n_elements(violators) le 0) then goto, leave_loop
 i = violators[n_elements(violators)-1] +1
 if (i le 0) then goto, leave_loop
 b = w[i] + w[i-1]
 a = (f[i]*w[i] + f[i-1]*w[i-1])/b
; Now, I'm shrinking the index.
 result = index->Get_Item(i,/DEREFERENCE)
 result2 = index->Get_Item(i-1,/DEREFERENCE)
 add_arr = append_arr(result, result2)
 index->Replace_Item,i-1,add_arr
 index->Delete,i
; Now, I'm shrinking the function
 f[i-1] = a
 new_f = shrink(f,i)
 f = new_f
;Now, I'm shrinking w
 w[i-1] = b
 v = shrink(w,i)
 w = v
ENDWHILE

leave_loop:
fstar = dblarr(m)
FOR I = 0, n_elements(f) -1 DO BEGIN
  vector = index->Get_Item(i,/DEREFERENCE)
  fstar[vector] = f[i]
  fstar[where(fstar le 0.0)] = 0.0
ENDFOR
riskstar = mean (sigmasq * fstar^2 + (z^2- sigmasq)*(1-fstar)^2)
xihat = fstar*z
fhat = transpose(U) ## xihat
;print, xihat
;Change name for output
thetas = xihat

IF KEYWORD_SET(plot) THEN BEGIN
 ;Now, plot out the results (nicely and in color please)
 plotsym,0,1,/FILL
 plot,x,fhat,psym=8, color = 0, /nodata
 oplot, x, y, psym=8, color = 2
 oplot, x, fhat,psym=8,color=4
 device,/close
ENDIF
yfit = transpose(fhat)

Obj_destroy, index
riskhat = mean( sigmasq * fstar^2  + (z^2 - sigmasq)*(1 - fstar)^2 )

;;The above is equivalent to the below
;;risk_hat = dblarr(41)
;;FOR I = 0, 40 DO BEGIN
;;    risk_hat[I] = risk_hat[I] + sigmasq[I] * fstar[I]^2  + (z[I]^2 - sigmasq[I])*(1 - fstar[I])^2
;;ENDFOR
;;risk_hat_tot = mean(risk_hat)
;;print, riskhat, risk_hat_tot
; Now we're after tauF which will be used
; to determine the confidence ball radius

noise_arr = diag(noise_cov)
B1 = U ## noise_arr ## transpose(U)
B2 = get_diag(B1)
B2_diag = diag(B2)
tr_row = dblarr(m,m)
FOR KK = 0, m-1 DO BEGIN  
 tr_row[*,kk] = ( ((2.d0 * fstar[kk] - 1.d0) * B1[kk,*])^2.d0 )
ENDFOR
tr = total(tr_row)
uu = B2 * (z^2 - B2) * (1.0 - fstar)^2
vv = (z * (1.0 - fstar))

;Below is is the radius^2 of the ball
tau = 2 * sqrt( total(uu)  +  transpose(vv) ## (B1 - B2_diag) ## vv  + 0.5 * tr )/sqrt(n)

;Now, calculate alpha% Confidence radii
zalpha = abs(gauss_cvf(double(alpha)))
radius = sqrt(double(m)*riskhat + sqrt(double(m))*tau*zalpha)
RETURN

END
