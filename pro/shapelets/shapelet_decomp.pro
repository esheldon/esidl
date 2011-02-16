;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; DECOMPOSE A GALAXY IMAGE INTO GAUSS-HERMITE FUNCTIONS
; (SHAPELETS). ONE CAN EITHER INPUT AN SDSS STRUCTURE CONTAINING THE
; NECESSARY INFORMATION (RUN, RERUN, CAMCOL, ETC.), OR A SINGLE
; IMAGE. A STRUCTURE CONTAINING THE SHAPELET PARAMETERS IS
; RETURNED. IF A PSF IMAGE EXISTS, THEN THE SHAPELET COEFFICIENTS ARE
; CORRECTED FOR THE PSF BY PROJECTING ONTO A NEW BASIS SPANNED BY THE
; SHAPELETS CONVOLVED WITH THE PSF.
;
; AUTHOR : BRANDON C. KELLY, STEWARD OBSERVATORY, JAN. 2005
;
; INPUTS :
;   NMAX - The maximum shapelet order (n1 + n2) to decompose the
;          image up to.
;
; OPTIONAL INPUTS :
;   GAL - The SDSS structure containing information for finding the
;         galaxy atlas image and PSF. This is passed on to GET_ATLAS,
;         so any tags that are required by GET_ATLAS must be included
;         in GAL.
;   CLR - The index of the SDSS band, u, g, r, i, or z (0, 1, 2, 3, or
;         4).
;   DTYPE - A string variable specifying the decomposition method. Can
;           be one of the following :
;
;       APPROX - Estimate the shapelet coefficients by approximating
;                the orthogonality properties of the shapelet
;                modes. This is done by expanding the galaxy image out
;                to a size of 20 shapelet scales. The image is not
;                weighted by it's standard deviation before
;                decomposing. In general this is a good approximation
;                so long as NMAX is not too large (e.g., > 25). This
;                can give more stable estimates of the shapelet
;                coefficients when NMAX is large, as compared to the
;                OLS estimate. In general, this is also the fastest
;                method, but that depends on the scale of the galaxy
;                to be decomposed. Also, the APPROX method does not
;                correct for the PSF, as convolving the shapelet
;                states with the PSF would mean that orthogonality
;                would not longer be a good approximation.
;       OLS - The shapelet coefficients are estimated using the
;             ordinary least squares method (OLS). This will give more
;             accurate estimates of the coefficients when NMAX is not
;             too large, but can give wildly variable estimates of the
;             coefficients when NMAX becomes ~ 15-20. The image is
;             weighted by its standard deviation, so the uncertainties
;             are taken into account during the decomposition. One
;             should check the shapelet coefficients on output if NMAX
;             is large to ensure the stability of the estimate; just
;             because the result may give an accurate estimate of the
;             true galaxy morphology, does not mean that the OLS
;             estimate of the shapelet coefficients gives an accurate
;             estimate of the galaxy's true shapelet
;             coefficients. Also, the shapelet states are convolved
;             with the SDSS PSF, giving an estimate of the deconvolved
;             shapelet coefficients. If the PSF is large then the OLS
;             estimate will be unstable at smaller NMAX.
;       PCR - Perform principal component regression for shapelet
;             decomposition. This is done by computing the principal
;             components of the shapelet basis matrix, and then
;             decomposing in this new basis. The number of principal
;             components used can be set manually or chosen
;             automatically. If chosen manually, then the user must
;             input a value for SVTHRESH, otherwise the number of
;             principal components is that which minimizes Stein's
;             Unbiased Risk Estimation (SURE). Once the principal
;             component expansion is estimated, the coefficients are
;             then transformed into the original shapelet basis. This
;             is the most stable and accurate method of computing the
;             shapelet coefficients, but in general the slowest. The
;             shapelet states are convolved with the SDSS PSF, giving
;             an estimate of the deconvolved shapelet coefficients.
;
;           The default method is OLS.
;
;    SVTHRESH - The theshold on the singular values to use if one
;               desires PCR. Denote the singular values by d_j. All
;               singular values for which d_j / max(d_j) >= SVTHRESH
;               is true are kept, otherwise their corresponding
;               coefficients are set to zero. Common choices are ~
;               1e-3. Note that the singular values are the
;               square-root of the eigenvalues associated with the
;               principal components, and the variance in each
;               principal component is d_j^2 / n.
;    IMAGE - A galaxy image. If input, then GAL is ignored and IMAGE
;            is decomposed instead.
;    NO_PSF_CORRECT - If set, then will not correct for effects of the
;                     PSF, i.e., if a PSF image exists and
;                     NO_PSF_CORRECT = 1 then the shapelet basis
;                     functions will not be convolved with the PSF in
;                     estimating the shapelet coefficients. The
;                     default is to convolve the shapelets with the
;                     PSF if a PSF image exists. If the shapelet
;                     states are not convolved with the PSF then the
;                     routine finds the optimal shapelet scale for
;                     decomposition.
;    PSF_IN - The PSF. This is only used if IMAGE is also input. If
;             GAL is input, then the PSF information is taken from
;             GAL.
;    NMAX_PSF - Like, NMAX, but for the PSF. The default is 10.
;
; OUTPUTS :
;    DECOMP - The decomposition structure. DECOMP contains the
;             following tags :
;  
;       RA, DEC - The RA and DEC of the galaxy(s). This is copied from
;                 the SDSS GAL structure, so these will be left blank
;                 if one is just inputing an image.
;       RUN, RERUN, CAMCOL, FIELD, ID - The SDSS parameters specifying
;                                       the galaxy, taken from the
;                                       input structure GAL.
;       BAND - The SDSS photometric band that was used, u, g, r, i, or
;              z.
;       NMAX - The maximum order used in the decomposition, max(n1 +
;              n2).
;       SCALE - The shapelet scale used in the decomposition, in units
;               of pixels. This is initially estimated from the SDSS
;               adaptively-weighted moments if GAL is supplied, or
;               from the RMS radius of the input image. After the
;               shapelet coefficients have been calculated, the scale
;               which is found to minimize the mean-squared error
;               between the original image and the reconstructed
;               image on a grid of scales between a factor of 0.5 and
;               2.0 of the original scale is chosen.
;       COEF - The shapelet coefficients.
;       N1, N2 - The shapelet order in the horizontal (x) and vertical
;                (y) directions. COEF[j] is the shapelet coefficient
;                for the |n1[j],n2[j]> state.
;       C_p - Mallow's C_p statistic. This is an unbiased estimate of
;             the loss between the true galaxy image (i.e., the image
;             before it was contaminated with noise) and the one
;             reconstructed from the shapelet coefficients. To be more
;             precise, consider the usual Gaussian model y_i = mu_i +
;             e_i, where mu_i is the true number of counts in the i^th
;             pixel, y_i is the observed number, and e_i is a random
;             deviate drawn from a normal density of mean zero and
;             standard deviation sigma_i. Sigma_i depends on the sky
;             noise and the number of counts in the pixel. Consider
;             the estimate, muhat_i, which is the shapelet fit to the
;             i^th pixel. The expected loss under this model is :
;
;             E((mu - muhat)^2 / (n * sigma^2)) = E((y - muhat)^2 /
;                  (n * sigma^2)) - 1 + 2p / n
;
;             Here, E(.) denotes expectation value, n is the number of
;             pixels, and p is the number of shapelet states
;             used. Above, I have left out the subscript i, the
;             quantities in the expecations are actually summed over
;             i. Note that the first quantity on the right is the
;             usual chi-square statistic. Mallow's C_p statistic is
;             defined as :
;
;             C_p = chi-square - 1 + 2p / n,
;
;             which is an unbiased estimate of the true loss.
;       PSF_NMAX, PSF_COEF, PSF_SCALE - The same as above, but for the
;                                       PSF. This is only calculated
;                                       if GAL or PSF is input. The
;                                       PSF coefficients are
;                                       calculated using the same
;                                       DYTYPE as the galaxy coefs,
;                                       but are ignored if APPROX is
;                                       used.
;       DTYPE - The decomposition method used, will be 'APPROX',
;               'OLS', or 'PCR'.
;       NPC - If DTYPE = 'PCR' then NPC is the number of principal
;             components used in the decomposition. This is the number
;             of parameters estimated in the fit.
;       
;
; CALLED ROUTINES :
;    GET_ATLAS, POSITIVE, CMP_DMOM, RESIZE_IMAGE, PHIN, CMREPLICATE,
;    FETCH_DIR, PSFIELD_NAME, READ_PSF, PSF_REC, IQR, SPLIT_ARRAY, CONVOLVE
;
; MODIFICATION HISTORY : 
;    Brandon Kelly (BK) - Convolve Shapelets with PSF and decompose in
;                         this basis for OLS, NSS and PCR, added PCR
;                         option (Feb. 2005)
;    BK - Changed C_p to give an estimate of the mean-squared error
;         instead of the total quadratic loss, easier to compare with
;         other decomposition with different number of pixels (Feb
;         2005)
;    BK - Changed Approx Method to return the reconstructed image
;         centered on the original input image with the original
;         dimensions. (Feb 2005)
;    BK - Changed ESS's ortho statistic so we don't have to compute
;         the entire Gram matrix. Now just calculate the diagonal of
;         the gram matrix, much faster. (Feb 2005)
;    BK - Changed OLS so the OLS solution is now done via Cholesky
;         decomposition, instead of the LU decomposition as
;         before. Cholesky decomposition is faster. (Feb 2005)
;    BK - Removed NSS option. (Feb 2005)
;    BK - Made CP_MSE, CONSTRUCT_GH_MATRIX, and
;         SHAPELET_SCALE_TRANSFORM stand-alone functions. (Feb 2005)
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;                           DECOMPOSITION METHODS
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro _shapelets_split_array, ind, ind1, ind2

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
    print, 'Syntax- split_array, ind, ind1, ind2'
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


;Function to compute the shapelet coefficients by exploiting the
;(approximate) orthogonality of the basis functions

pro approx_decomp, image, nmax, x0, scale, decomp, skysig=skysig, $
                   image_var=image_var, ghfit=ghfit, ghmatrix=gh, $
                   max_size=max_size, status=status

  status=1

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
      return
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
      return
  ENDIF 

  ;; making a copy which we will modify
  image0 = image
  image_var0 = image_var
  x00 = x0

  ;; Expand the size if needed
  IF osize[0] NE nx0[0] OR osize[1] NE nx0[1] THEN BEGIN
      resize_image, image, x0, osize ;resize the image
      resize_image, image_var, x0, osize
      
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
  GH = construct_gh_matrix(nx, x0, scale, nmax, n1=n1, n2=n2)

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
      coefs_k = shapelet_scale_transform( coefs, scale, gammas[k], n1, n2 )
                                ;construct new basis with this scale
      GH_k = construct_gh_matrix(nx0, x00, gammas[k], nmax)
      ghfit = coefs_k ## GH_k
      ghfit = reform(ghfit, nx0[0], nx0[1])

      C_p = cp_mse( image0 / sqrt(image_var0), ghfit / sqrt(image_var0), m ) ;estimate of the MSE
      error[k] = C_p

  endfor

;now use a spline to interpolate MSE to get best scale
  nscales2 = 200
  gammas2 = scale / 2.0 + dindgen(nscales2) / (nscales2 - 1.0) * scale
  error2 = interpol(error, gammas, gammas2, /spline)
  minerr = min(error2, minind)  ;find minimum error
;get shapelet coefficients for this scale
  coefs_k = shapelet_scale_transform( coefs, scale, gammas2[minind], n1, n2 )
  GH_k = construct_gh_matrix(nx0, x00, gammas2[minind], nmax)

;find MSE in model with best scale
  ghfit = coefs_k ## GH_k
  ghfit = reform(ghfit, nx0[0], nx0[1])

;use Mallow's C_p statistic to estimate the mean squared error (MSE)
  C_p = cp_mse( image0 / sqrt(image_var0), ghfit / sqrt(image_var0), m )

  image = image0                ;resore inputs
  x0 = x00
  decomp = {coef:coefs_k, scale:gammas2[minind], n1:n1, n2:n2, C_p:C_p}

  status=0
  return
end


pro approx_decomp_old, image, nmax, x0, scale, decomp, skysig=skysig, $
                   image_var=image_var, ghfit=ghfit, ghmatrix=gh, $
                   max_size=max_size, status=status

  status=1

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
      return
  ENDIF 

  szfact = 20.0             ;number of scales that the image will be resized to
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
      return
  ENDIF 

  ;; making a copy which we will modify
  image0 = image
  x00 = x0

  ;; Expand the size if needed
  IF osize[0] NE nx0[0] OR osize[1] NE nx0[1] THEN BEGIN
      resize_image, image, x0, osize ;resize the image
      resize_image, image_var, x0, osize
      
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
  GH = construct_gh_matrix(nx, x0, scale, nmax, n1=n1, n2=n2)

;GH matrix is approximately orthogonal, use this to estimate the
;shapelet coefficients
  coefs = (image[*]) ## transpose(GH)

;get reconstructed image
  ghfit = coefs ## GH
  ghfit = reform(ghfit, nx[0], nx[1]) ;make into an array

  m = n_elements(coefs)         ;number of parameters

  IF osize[0] NE nx0[0] OR osize[1] NE nx0[1] THEN BEGIN
;resize the images back to the original size
      resize_image, ghfit, x0, nx0
      resize_image, image, x0, nx0
      resize_image, image_var, x0, nx0
  ENDIF

;now find best scale on a grid of scale values from 0.5*scale to 1.5*scale
  nscales = 10
  gammas = scale / 2.0 + dindgen(nscales) / (nscales - 1.0) * scale

  error = fltarr(nscales)

  for k = 0, nscales - 1 do BEGIN
                                ;transform the shapelet coefficients
                                ;to have scale gammas[k]
      coefs_k = shapelet_scale_transform( coefs, scale, gammas[k], n1, n2 )
                                ;construct new basis with this scale
      GH_k = construct_gh_matrix(nx, x0, gammas[k], nmax)
      ghfit = coefs_k ## GH_k
      ghfit = reform(ghfit, nx[0], nx[1])
      IF osize[0] NE nx0[0] OR osize[1] NE nx0[1] THEN $
        resize_image, ghfit, x0, nx0 ;resize the image back to original size

      C_p = cp_mse( image / sqrt(image_var), ghfit / sqrt(image_var), m ) ;estimate of the MSE
      error[k] = C_p

  endfor

;now use a spline to interpolate MSE to get best scale
  nscales2 = 200
  gammas2 = scale / 2.0 + dindgen(nscales2) / (nscales2 - 1.0) * scale
  error2 = interpol(error, gammas, gammas2, /spline)
  minerr = min(error2, minind)  ;find minimum error
;get shapelet coefficients for this scale
  coefs_k = shapelet_scale_transform( coefs, scale, gammas2[minind], n1, n2 )
  GH_k = construct_gh_matrix(nx, x0, gammas2[minind], nmax)

;find MSE in model with best scale
  ghfit = coefs_k ## GH_k
  ghfit = reform(ghfit, nx[0], nx[1])
  IF osize[0] NE nx0[0] OR osize[1] NE nx0[1] THEN $
    resize_image, ghfit, x0, nx0 ;resize the image back to original size

;use Mallow's C_p statistic to estimate the mean squared error (MSE)
  C_p = cp_mse( image / sqrt(image_var), ghfit / sqrt(image_var), m )

  image = 0
  image = temporary(image0)     ;resore inputs
  x0 = x00
  decomp = {coef:coefs_k, scale:gammas2[minind], n1:n1, n2:n2, C_p:C_p}

  status=0
  return
end

;Function to compute the shapelet coefficients by minimizing the
;chi-squared of the observed image with the image reconstructed from
;the first m shapelets. This is the ordinary least-squares method of
;computing the coefficients.
pro ols_decomp, image, nmax, x0, scale, decomp, image_var=image_var, $
                ghfit=ghfit, ghmatrix=gh, psf=psf

  nx = size(image, /dim)
  if n_elements(image_var) eq 0 then image_var = replicate(1.0, nx[0], nx[1])
                                ;construct the shapelet
                                ;(Gauss-Hermite) basis
  GH = construct_gh_matrix(nx, x0, scale, nmax, n1=n1, n2=n2)

  m = n_elements(n1)            ;number of coefficients to estimate

  image = (image / sqrt( image_var ))[*] ;normalize by the noise and make into a vector

  if n_elements(psf) gt 0 then begin
                                ;convolve shapelet basis with the PSF
                                ;and normalize by the noise. this
                                ;creates a new basis, E, which we
                                ;decompose the project the image onto
      E = dblarr(nx[0] * nx[1], m)
      phi = reform(GH[*,0], nx[0], nx[1]) ;get 2-d representation of first shapelet
      E[0,0] = (convolve(phi, psf, ft_psf = psf_fft))[*] ;do convolutions with FFT, faster
      for j = 1, m - 1 do begin
          phi = reform(GH[*,j], nx[0], nx[1])
          E[0,j] = (convolve(phi, psf, ft_psf = psf_fft))[*] ;zero-index here for speed,
                                                             ;this is really a vector
      endfor
      
  endif else E = GH     ;no input PSF, just use regular shapelet basis
  
  Estar = E / cmreplicate( (sqrt( image_var ))[*], m ) ;now normalize by the noise

                                ;solve via Cholesky Decomposition for
                                ;speed
  Etrans = transpose(Estar)
  Gram = Estar ## Etrans

  corr = image ## Etrans
  choldc, Gram, P
  coefs = cholsol( Gram, P, corr )

  if n_elements(psf) eq 0 then begin ;can optimize the shapelet scale if E = GH
                                ;now find best scale on a grid of
                                ;scale values from 0.5*scale to
                                ;1.5*scale
      nscales = 10
      gammas = scale / 2.0 + dindgen(nscales) / (nscales - 1.0) * scale
      
      error = fltarr(nscales)
      
      for k = 0, nscales - 1 do begin
                                ;transform the shapelet coefficients
                                ;to have scale gammas[k]
          coefs_k = shapelet_scale_transform( coefs, scale, gammas[k], n1, n2 )
                                ;construct new basis with this scale
          GH_k = construct_gh_matrix(nx, x0, gammas[k], nmax)
          GH_k = GH_k / cmreplicate( (sqrt( image_var ))[*], m )   
          
          ghfit = coefs_k ## GH_k
          
          C_p = cp_mse( image, ghfit, m )
          error[k] = C_p
          
      endfor
                                ;now use a spline to interpolate MSE to get best scale
      nscales2 = 200
      gammas2 = scale / 2.0 + dindgen(nscales2) / (nscales2 - 1.0) * scale
      error2 = interpol(error, gammas, gammas2, /spline)
      minerr = min(error2, minind) ;find minimum error
                                ;get shapelet coefficients for this scale
      coefs_k = shapelet_scale_transform( coefs, scale, gammas2[minind], n1, n2 )
      GH_k = construct_gh_matrix(nx, x0, gammas2[minind], nmax)
      
                                ;find MSE in model with best scale
      ghfit = coefs_k ## (GH_k / cmreplicate( (sqrt( image_var ))[*], m ))
      coefs = coefs_k
      
      scale = gammas2[minind]
      E = GH_k

  endif else ghfit = coefs ## Estar

;use Mallow's C_p statistic to estimate the mean squared error (MSE)
  C_p = cp_mse( image, ghfit, m )
  
  ghfit = coefs ## E ;get reconstructed image, not normalized by noise or corrected for PSF
  
  ghfit = reform(ghfit, nx[0], nx[1])

  image = reform(image, nx[0], nx[1]) * sqrt(image_var) ;restore input

  decomp = {coef:coefs, scale:scale, n1:n1, n2:n2, C_p:C_p}

  return
end

;Function to estimate the shapelet coefficients via principal
;component regression.
pro pcr_decomp, image, nmax, x0, scale, decomp, image_var=image_var, $
                ghfit=ghfit, ghmatrix=gh, svthresh=svthresh, psf=psf

  nx = size(image, /dim)
  if n_elements(image_var) eq 0 then image_var = replicate(1.0, nx[0], nx[1])
  if n_elements(svthresh) eq 0 then begin
      svthresh = 0.0 
      auto = 1
  endif else auto = 0
                                ;construct the shapelet
                                ;(Gauss-Hermite) basis
  GH = construct_gh_matrix(nx, x0, scale, nmax, n1=n1, n2=n2)

  m = n_elements(n1)            ;number of coefficients to estimate

  image = (image / sqrt( image_var ))[*] ;normalize by the noise and make into a vector

  if n_elements(psf) ne 0 then begin
                                ;convolve shapelet basis with the PSF
                                ;and normalize by the noise. this
                                ;creates a new basis, E, which we
                                ;decompose the project the image onto
      E = dblarr(nx[0] * nx[1], m)
      phi = reform(GH[*,0], nx[0], nx[1]) ;get 2-d representation of first shapelet
      E[0,0] = (convolve(phi, psf, ft_psf = psf_fft))[*] ;do convolutions with FFT, faster
      for j = 1, m - 1 do begin
          phi = reform(GH[*,j], nx[0], nx[1])
          E[0,j] = (convolve(phi, psf, ft_psf = psf_fft))[*] ;zero-index here for speed,
                                ;this is really a vector
      endfor
      
  endif else E = GH             ;no input PSF, just use regular shapelet basis

  Estar = E / cmreplicate( (sqrt( image_var ))[*], m ) 
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
  Z = Etrans ## (V / cmreplicate(d, m)) ;note that Z is same as U in SVD
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
          coefs_k = shapelet_scale_transform( coefs, scale, gammas[k], n1, n2 )
                                ;construct new basis with this scale
          GH_k = construct_gh_matrix(nx, x0, gammas[k], nmax)
          GH_k = GH_k / cmreplicate( (sqrt( image_var ))[*], m )   
          
          ghfit = coefs_k ## GH_k
          
          C_p = cp_mse( image, ghfit, p )
          error[k] = C_p
          
      endfor
                                ;now use a spline to interpolate MSE
                                ;to get best scale
      nscales2 = 200
      gammas2 = scale / 2.0 + dindgen(nscales2) / (nscales2 - 1.0) * scale
      error2 = interpol(error, gammas, gammas2, /spline)
      minerr = min(error2, minind) ;find minimum error
                                ;get shapelet coefficients for this scale
      coefs_k = shapelet_scale_transform( coefs, scale, gammas2[minind], n1, n2 )
      GH_k = construct_gh_matrix(nx, x0, gammas2[minind], nmax)
      
                                ;find MSE in model with best scale
      ghfit = coefs_k ## (GH_k / cmreplicate( (sqrt( image_var ))[*], m ))
      coefs = coefs_k

      scale = gammas2[minind]
      E = GH_k
      
  endif else ghfit = coefs ## Estar

;use Mallow's C_p statistic to estimate the mean squared error (MSE)
  C_p = cp_mse( image, ghfit, p )

  ghfit = coefs ## E ;get reconstructed image, not normalized to noise or corrected for PSF

  ghfit = reform(ghfit, nx[0], nx[1])

  image = reform(image, nx[0], nx[1]) * sqrt(image_var) ;restore input

  decomp = {coef:coefs, scale:scale, n1:n1, n2:n2, C_p:C_p, npc:p}

  return
end



;Function to compute the shapelet coefficients by minimizing
;chi-squared of the true image and the reconstructed one. This is done
;by finding the p < m nested subset of the shapelet basis
;that minimize Mallow's C_p. This is the most accurate decomposition
;technique but is also the slowest.
pro nss_decomp, image, nmax, x0, scale, decomp, image_var=image_var, $
                ghfit=ghfit, ghmatrix=gh

  nx = size(image, /dim)
                                ;construct the shapelet
                                ;(Gauss-Hermite) basis
  GH = construct_gh_matrix(nx, x0, scale, nmax, n1=n1, n2=n2)

  m = n_elements(n1)            ;number of coefficients to estimate

  image = (image / sqrt( image_var ))[*] ;normalize by the noise and make into a vector
  GH = GH / cmreplicate( (sqrt( image_var ))[*], m )

;get best nested subset via forward selection
  p = 0L
  ghfit = 0d
  set_IA = lindgen(m)           ;the inactive set
  minerr = 1d30
  coefs = dblarr(m)
  error = dblarr(m)

  repeat begin
      print, float(p) / m
      if p gt 0 then _shapelets_split_array, lindgen(m), set_A, set_IA ;get new inactive set
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

      error[p-1] = cp_mse( image, ghfit, p ) ;estimate of the MSE

      if error[p-1] lt minerr then begin
          coefs[set_A] = coefs_k
          minerr = error[p-1]
          nss = set_A           ;the best nested subset
      endif

  endrep until p eq m         ;stop when full least squares solution is reached

  if n_elements(nss) ne m then $
    _shapelets_split_array, lindgen(m), nss, set_IA ;get the inactive set for the best NSS

;now find best scale on a grid of scale values from 0.5*scale to
;1.5*scale
  nscales = 20
  gammas = scale / 2.0 + dindgen(nscales) / (nscales - 1.0) * scale

  m = n_elements(nss)  ;the number of shapelet states in the best nested subset
  error = fltarr(nscales)

  for k = 0, nscales - 1 do begin
                                ;transform the shapelet coefficients
                                ;to have scale gammas[k]
      coefs_k = shapelet_scale_transform( coefs, scale, gammas[k], n1, n2 )
                                ;construct new basis with this scale
      GH_k = construct_gh_matrix(nx, x0, gammas[k], nmax)
;    GH_k = GH_k / cmreplicate( (sqrt( image_var ))[*], m )

      ghfit = coefs_k[nss] ## GH_k[*,nss] ;only use those coefs in the best NSS
      ghfit = ghfit / (sqrt( image_var ))[*]

      C_p = cp_mse( image, ghfit, m )
      error[k] = C_p

  endfor

;now use a spline to interpolate MSE to get best scale
  nscales2 = 200
  gammas2 = scale / 2.0 + dindgen(nscales2) / (nscales2 - 1.0) * scale
  error2 = interpol(error, gammas, gammas2, /spline)
  minerr = min(error2, minind)  ;find minimum error
;get shapelet coefficients for this scale
  coefs_k = shapelet_scale_transform( coefs, scale, gammas2[minind], n1, n2 )
  coefs_k[set_IA] = 0d          ;set to zero those coefs not in the active set
  GH_k = construct_gh_matrix(nx, x0, gammas2[minind], nmax)

;find MSE in model with best scale
  ghfit = coefs_k[nss] ## GH_k[*,nss]
  ghfit = ghfit / (sqrt( image_var ))[*]

;use Mallow's C_p statistic to estimate the mean squared error (MSE)
  C_p = cp_mse( image, ghfit, m )

  ghfit = ghfit * (sqrt( image_var ))[*]
  ghfit = reform(ghfit, nx[0], nx[1])

  image = reform(image, nx[0], nx[1]) * sqrt(image_var) ;restore input
  decomp = {coef:coefs_k, scale:gammas2[minind], n1:n1, n2:n2, C_p:minerr}

  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;                               MAIN ROUTINE
;
; To Do: maybe keep uncorrected scale as well?
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro shapelet_decomp, nmax, decomp, gal, clr, $
                     dtype=dtype, $
                     image=image_in, $
                     recon=recon, ghmatrix=ghmatrix, $
                     psFieldInfo=psFieldInfo, max_size=max_size, $
                     status=status, $
                     svthresh=svthresh, $
                     psf_in=psf_in, nmax_psf=nmax_psf, $
                     correct_psf=correct_psf, $
                     trim_atlas=trim_atlas, $
                     trim_nsigma=trim_nsigma, $
                     trimmed_image=trimmed_image

  on_error, 2

  if n_params() lt 2 then begin
      print, 'Syntax- shapelet_decomp, nmax[, decomp, gal, clr, dtype=dtype, image=image, recon=recon, ghmatrix=, psFieldInfo=, max_size=, status=, svthresh=svthresh, psf_in=psf_in, nmax_psf=nmax_psf, /correct_psf, /trim_atlas]'
      return
  endif

  status = 1

  delvarx, decomp, recon, ghmatrix, psFieldInfo

  IF n_elements(trim_nsigma) NE 0 THEN BEGIN 
      IF NOT keyword_set(trim_atlas) THEN trim_atlas=1
  ENDIF 

  ;; Check input and set up defaults
  ngal = n_elements(gal)
  if ngal eq 0 then begin
      ngal = 1
      atlas = 0
  endif else atlas = 1

  npsf = n_elements(psf_in)

  if n_elements(dtype) eq 0 then dtype = 'ols'
  dtype = strlowcase(strtrim(dtype))
  if dtype ne 'ols' and dtype ne 'approx' and dtype ne 'pcr' then begin
      print, 'If input, DTYPE must be one of the following : '
      print, '   OLS, APPROX, PCR'
      message,'Halting'
  ENDIF
  IF n_elements(svthresh) GT 0 THEN BEGIN 
      IF (svthresh GT 1 OR svthresh LT 0) THEN BEGIN 
          print, 'SVTHRESH must be less than or equal to 1 ' + $
            'and greater or equal to 0.'
          message,'Halting'
      ENDIF 
  ENDIF

  IF atlas THEN BEGIN 
      colors = ['u','g','r','i','z']
      IF clr LT 0 OR clr GT 4 THEN BEGIN 
          print, 'CLR must be an integer in [0,4]'
          message,'Halting'
      ENDIF 
      band = colors[clr]
  ENDIF ELSE band=''

  ;; number of coeffs given nmax
  m = long((nmax + 1) * (nmax + 2) / 2)
  psf_nmax = n_elements(nmax_psf) eq 0 ? 10 : nmax_psf

  ;; psf saved info was scalar. E.S.S.
  mpsf = long((psf_nmax + 1) * (psf_nmax + 2) / 2)
  decomp_bprint = { nmax:nmax,                                $
                    scale:-1d,                                $
                    coef:dblarr(m),                           $
                    n1:intarr(m), n2:intarr(m),               $
                    flux:0d,                                  $
                    $
                    center: [0d, 0d],                         $ 
                    size: [0,0],                              $
                    $
                    C_p:9999d,                                $
                                                              $
                    psf_nmax:psf_nmax,                        $
                    psf_scale:-1d,                            $
                    psf_coef:dblarr(mpsf),                    $
                    psf_n1:intarr(mpsf), psf_n2:intarr(mpsf), $
                                                              $
                    dtype:dtype                               $
                  }

  IF dtype EQ 'approx' THEN BEGIN 
      decomp_bprint = create_struct(decomp_bprint,'ortho',-1d)
  ENDIF 
  IF dtype EQ 'pcr' THEN BEGIN 
      decomp_bprint = create_struct(decomp_bprint,'npc',0)
  ENDIF

  ;; Add in extra tags for atlas images. Save anough info that
  ;; we can just send the output struct back into this program

  IF atlas THEN BEGIN 
      addstruct = { $
                    run:0L, $
                    rerun:0, $
                    camcol:0, $
                    field:0, $
                    id:0L, $
                    rowc:fltarr(5), $
                    colc:fltarr(5), $
                    m_rr_cc:fltarr(5), $
                    m_rr_cc_psf:fltarr(5), $
                    clr:clr $
                  }

      ;; Not necessary, but copy if input
      IF tag_exist(gal, 'ra') AND tag_exist(gal, 'dec') THEN BEGIN 
          copyRADEC=1
          addStruct=create_struct(addStruct, 'ra', 0d, 'dec', 0d)
      ENDIF ELSE copyRADEC=0

      decomp_bprint = create_struct(addstruct,decomp_bprint)
  ENDIF 

  decomp = replicate(decomp_bprint, ngal)

;store useful info in the output structure
  if atlas then begin

      IF copyRADEC THEN BEGIN 
          decomp.ra = gal.ra
          decomp.dec = gal.dec
      ENDIF 

      decomp.run = gal.run
      decomp.rerun = gal.rerun
      decomp.camcol = gal.camcol
      decomp.field = gal.field
      decomp.id = gal.id

      decomp.rowc = gal.rowc[clr]
      decomp.colc = gal.colc[clr]
      
      decomp.m_rr_cc = gal.m_rr_cc[clr]
      decomp.m_rr_cc_psf = gal.m_rr_cc_psf[clr]

  endif

;decompose each galaxy one-at-a-time
  for i = 0, ngal - 1 do begin

      if atlas then begin       ;need to grap the atlas image
          
          ;; Use read_atlas instead of old get_atlas
          CASE clr OF
              0: read_atlas, gal[i], $
                row0=row0, col0=col0, imu=image, status=astatus, /silent
              1: read_atlas, gal[i], $
                row0=row0, col0=col0, img=image, status=astatus, /silent
              2: read_atlas, gal[i], $
                row0=row0, col0=col0, imr=image, status=astatus, /silent
              3: read_atlas, gal[i], $
                row0=row0, col0=col0, imi=image, status=astatus, /silent
              4: read_atlas, gal[i], $
                row0=row0, col0=col0, imz=image, status=astatus, /silent
          ENDCASE 

          ;; Was the atlas image found?
          IF astatus EQ 0 THEN BEGIN 
              
              atlas_dir, gal[i].run, gal[i].rerun, gal[i].camcol, atldir
              
;get the PSF data
              
              psp = sdss_read('psfield',$
                              gal[i].run,gal[i].camcol,fields=gal[i].field)
              psf_rec, psp, gal[i].rowc, gal[i].colc, psfp, clr ;get psf image
              psf = *psfp[0]
              npsf = n_elements(psf)

              skysig = (*psp[5]).skysig[clr] ;get the sky variance

                                ;get the galaxy center
              rowctr = gal[i].rowc[clr] - row0[clr]
              colctr = gal[i].colc[clr] - col0[clr]
              
              ;; Trim image of zeros if necessary.  refind the center
              IF keyword_set(trim_atlas) THEN BEGIN 

                  nskysig = 2
                  trimmed_image = $
                    shapelet_trim_image(image, 0.0, skysig, nskysig, $
                                        smoothing_window=4, slack=2, $
                                        minx=minx, miny=miny, status=trstatus)
                  IF trstatus NE 0 THEN BEGIN 
                      delvarx, recon, ghmatrix
                      return
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

                                ;initial guess for the galaxy scale
                                ;based adaptively-weighted moments
              scale = sqrt( gal[i].m_rr_cc[clr] / 2.0 )
              psf_scale = sqrt( gal[i].m_rr_cc_psf[clr] )

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

              ;; E.S.S.
              ;; This can be a significant memory leak!
              ptr_free, psp, psfp

;              imsize = size(image, /dim)
;              flux = total( image, /nan )
;              Ixx = total( total(image, 2) * (dindgen(imsize[0]) - x0[0])^2 ) / flux > 0
;              Iyy = total( total(image, 1) * (dindgen(imsize[1]) - x0[1])^2 ) / flux > 0
;              scale2 = sqrt((Ixx + Iyy) / 2.0)
;              print,'scale = ',scale
;              print,'scale2 unweighted = ',scale2

          ENDIF 

      endif else begin          ;just decompose input image

          image = image_in
          nx = size(image, /dim)
                                ;estimate sky noise from pixels along border
          sky = [image[*,0], image[*,nx[1]-1], reform(image[0,*]), reform(image[nx[0]-1,*])]
          iqr = iqr(sky, skysig) ;robust estimate of the sky noise
          smooth_image = median( image, 3 ) ;median-smooth the image to find center
          maxflux = max( smooth_image, center )
          rowctr = center / nx[0]
          colctr = center - rowctr * nx[0]
          x0 = [colctr, rowctr] ;the galaxy center
                                ;find the galaxy RMS width to use as
                                ;the scale
          flux = total( image, /nan )
          Ixx = total( total(image, 2) * (dindgen(nx[0]) - x0[0])^2 ) / flux > 0
          Iyy = total( total(image, 1) * (dindgen(nx[1]) - x0[1])^2 ) / flux > 0
          scale = sqrt((Ixx + Iyy) / 2.0)

          if scale eq 0 then scale = sqrt( (nx[0] + nx[1]) / 2.0 )

          if npsf gt 0 then begin ;PSF was input, so do same estimates, 
                                                ;but assume same sky noise
              psf = psf_in
              nxpsf = size(psf, /dim)

              smooth_image = median( psf, 3 ) ;median-smooth the psf to find center
              maxflux = max( smooth_image, center )
              rowctr = center / nxpsf[0]
              colctr = center - rowctr * nxpsf[0]
              psfx0 = [colctr, rowctr] ;the galaxy center
                                ;find the galaxy RMS width to use as
                                ;the scale
              flux = total( psf, /nan )
              Ixx = total( total(psf, 2) * (dindgen(nxpsf[0]) - psfx0[0])^2 ) / flux > 0
              Iyy = total( total(psf, 1) * (dindgen(nxpsf[1]) - psfx0[1])^2 ) / flux > 0
              psf_scale = sqrt((Ixx + Iyy) / 2.0)

              if psf_scale eq 0 then psf_scale = sqrt( (nxpsf[0] + nxpsf[1]) / 2.0 )

          endif

          astatus = 0

      endelse

      ;; Status != 0 when atlas image not found.
      IF astatus EQ 0 THEN BEGIN 

          ;; create matrix of image variances

          ;; estimate of the poisson variance
          muhat = $
            0.5 * sqrt( 4.0 * positive( image^2 - skysig^2 ) + 1.0 ) - 0.5

          image_var = muhat + skysig^2 ;estimate of the variance in the image

          if atlas then begin
              wpsf = where(psf EQ 0, nwpsf)

              if nwpsf ne 0 then $
                psf[wpsf] = psf[wpsf] + skysig * randomn(seed, nwpsf)
          endif
          
          ;; decompose image based on method
          case dtype of
              'approx' : BEGIN 
                  approx_decomp, image, nmax, x0, scale, decomp0, $
                    skysig=skysig, image_var=image_var, $
                    ghfit=recon, ghmatrix=ghmatrix, max_size=max_size, $
                    status=dstatus
                  IF dstatus NE 0 THEN BEGIN 
                      delvarx, recon, ghmatrix
                      return
                  ENDIF 
                  IF npsf NE 0 THEN BEGIN 
                      ols_decomp, psf, psf_nmax, psfx0, psf_scale, dpsf
                  ENDIF 
              END 
              'ols' : begin
                  if npsf gt 0 then begin
                                ;correct shapelet scale for PSF, must
                                ;be at least of size one pixel

                      IF keyword_set(correct_psf) THEN BEGIN 
                          scale = sqrt(positive(scale^2 - psf_scale^2)) > 1.0
                      ENDIF 
                      ols_decomp, psf, psf_nmax, psfx0, psf_scale, dpsf, $
                        image_var=psf_var
                  endif
                  if keyword_set(correct_psf) then $
                    ols_decomp, image, nmax, x0, scale, decomp0, $
                    image_var=image_var, ghfit=recon, ghmatrix=ghmatrix, $
                    psf=psf/total(psf) $
                  else $
                    ols_decomp, image, nmax, x0, scale, decomp0, $
                    image_var=image_var, ghfit=recon, ghmatrix=ghmatrix
              end
              'pcr' : begin
                  if npsf gt 0 then begin
                                ;correct shapelet scale for PSF, must
                                ;be at least of size one pixel

                      IF keyword_set(correct_psf) THEN BEGIN 
;                          print,'centroid = ',x0
;                          print,'scale original = ',scale
                          scale = sqrt(positive(scale^2 - psf_scale^2)) > 1.0
;                          print,'scale new = ',scale
                      ENDIF 
                      pcr_decomp, psf, psf_nmax, psfx0, psf_scale, dpsf, $
                        image_var=psf_var
                  endif

                  ;; An intrinsic idl proc stupidly crashes instead of
                  ;; exiting with an error status.  We can catch the error
                  ;; by using execute

                  command = $
                    'pcr_decomp, image, nmax, x0, scale, decomp0, '+$
                    '   image_var=image_var, ghfit=recon, ghmatrix=ghmatrix,'+$
                    '   svthresh=svthresh'
                  IF keyword_set(correct_psf) THEN $
                    command = command+', psf=psf/total(psf)'

                  IF NOT execute(command) THEN BEGIN 
                      delvarx, recon, ghmatrix
                      return
                  ENDIF 

;                  if keyword_set(correct_psf) then $
;                    pcr_decomp, image, nmax, x0, scale, decomp0, $
;                    image_var=image_var, ghfit=recon, ghmatrix=ghmatrix, $
;                    svthresh=svthresh, psf=psf/total(psf) $
;                  else $
;                    pcr_decomp, image, nmax, x0, scale, decomp0, $
;                    image_var=image_var, ghfit=recon, ghmatrix=ghmatrix, $
;                    svthresh=svthresh
              end
          endcase
                    
          ;; store shapelet decomposition data into output structure
          decomp[i].scale = decomp0.scale
          decomp[i].coef = decomp0.coef
          decomp[i].n1 = decomp0.n1
          decomp[i].n2 = decomp0.n2
          decomp[i].C_p = decomp0.C_p
          if npsf NE 0 then begin
              decomp[i].psf_scale = dpsf.scale
              decomp[i].psf_coef = dpsf.coef
              decomp[i].psf_n1 = dpsf.n1
              decomp[i].psf_n2 = dpsf.n2

              psf_flux = shapelet_flux( dpsf.coef, dpsf.scale, $
                                        dpsf.n1, dpsf.n2 )

              ;;normalize PSF to have flux of unity
              decomp[i].psf_coef = decomp[i].psf_coef / psf_flux 
          endif
          
          ;; normalize coefficients to have flux
          ;; of unity (in counts)
          flux = shapelet_flux( decomp0.coef, decomp0.scale, $
                                decomp0.n1, decomp0.n2 )
          decomp[i].coef = decomp[i].coef / flux
          decomp[i].flux = flux

          decomp[i].center = x0
          decomp[i].size = (size(image))[1:2]

          IF dtype EQ 'approx' THEN BEGIN 
              ortho = 1.0d - total( total(GHMatrix^2, 1) ) / m
              decomp[i].ortho = ortho 
          ENDIF 

          IF dtype EQ 'pcr' THEN BEGIN 
              decomp[i].npc = decomp0.npc
          ENDIF 

      ENDIF ;; Don't proceed if no atlas image found

  endfor

  status=0
  return
end
