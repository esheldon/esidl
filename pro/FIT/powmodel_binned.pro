PRO powmodel_binned, normrange, powrange, nnorm, npower, $
                     xvals, nbin, ind1, ind2, $
                     modelx, model, norm, power, weights=weights, $
                     modelerr=modelerr

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: '
      print,'powmodel_binned, normrange, powrange, nnorm, npower, $'
      print,'            xvals, nbin, ind1, ind2, $'
      print,'            modelx, model, norm, power, weights=weights'
      return
  ENDIF  

  ;; xvals must be sorted!!!

  n1=n_elements(ind1)
  n2=n_elements(ind2)
  IF (n1 NE nbin) OR (n2 NE nbin) THEN message,'# elements in index arrays must equal nbin'

  nx=n_elements(xvals)
  nw=n_elements(weights)
  IF n_elements(weights) EQ 0 THEN BEGIN
      weights = replicate(1.0, nx)
  ENDIF ELSE BEGIN 
      IF nx NE nw THEN message,'# in weights must be same as # in xvals'
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; find mean x values in each bin
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  modelx = fltarr(nbin)
  FOR ib=0L, nbin-1 DO BEGIN 

      tx = xvals[ ind1[ib]:ind2[ib] ]
      wmom, tx, 1.0, xmean, xsig, xerr, $
        inputweights = weights[ ind1[ib]:ind2[ib] ]

      modelx[ib] = xmean

  ENDFOR 

  norm = arrscl( findgen(nnorm), normrange[0], normrange[1] )
  power = arrscl( findgen(npower), powrange[0], powrange[1] )

  model = fltarr(npower, nnorm, nbin )
  modelerr = model

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Find mean of model in these bins
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  FOR inorm=0L, nnorm-1 DO BEGIN 
      FOR ipow=0L, npower-1 DO BEGIN 

          ;; Msolar/pc^2
          FOR ib=0L, nbin -1 DO BEGIN 

              tmod = norm[inorm]*( xvals[ ind1[ib]:ind2[ib] ] )^power[ipow]
              wmom, tmod, 1.0, ymean, ysig, yerr, $
                inputweights = weights[ ind1[ib]:ind2[ib] ]

              model[ipow, inorm, ib] = ymean
              modelerr[ipow, inorm, ib] = yerr

          ENDFOR 
      ENDFOR 
  ENDFOR 

  return
END 
