PRO combine_wtheta, files, outstruct
  
  ;; this combines sum files

  IF n_params() LT 1 THEN BEGIN
      print,'-Syntax: combine_shear, files, outstruct'
      return
  ENDIF 

  nf=n_elements(files)
  FOR i=0, nf-1 DO BEGIN

      struct1 = mrdfits(files[i], 1, hdr, /silent)

      IF i EQ 0 THEN BEGIN

          struct = struct1
          nbin = n_elements(struct1.rsum)
          nbinuse = nbin
          struct1 = 0

      ENDIF ELSE BEGIN 

          nbin2 =  n_elements(struct1.rsum)
          IF nbin2 NE nbin THEN BEGIN 
              message,'Files must contain same number of bins'
          ENDIF 
          IF (  (struct.binsize NE struct1.binsize) OR $
                (struct.rmin NE struct1.rmin) ) THEN BEGIN 
              message,'Either binsize, rmin, or h has changed'
          ENDIF 

          struct.rsum = struct.rsum + struct1.rsum
          struct.npair = struct.npair + struct1.npair

          struct.nlenses = struct.nlenses + struct1.nlenses
          struct.totpairs = struct.totpairs + struct1.totpairs

          FOR jj=0L, nbin-1 DO BEGIN
              struct.rmax_act[jj] = max( [struct.rmax_act[jj], struct1.rmax_act[jj]] )
          ENDFOR 

          struct1 = 0
      ENDELSE                   ; check complete bin, also check if already done
      IF (struct.rsum[nbin-1] EQ 0.) AND (nbinuse NE nbin-1) THEN BEGIN 
          ;print,'Last bin is incomplete. Cutting'
          nbinuse=nbin-1
      ENDIF 
  ENDFOR 
 
  arrval = fltarr(nbinuse)
  outstruct = wthetastruct(arrval)

  outstruct.nlenses  = struct.nlenses
  outstruct.totpairs = struct.totpairs
  outstruct.binsize  = struct.binsize
  outstruct.rmin     = struct.rmin
  outstruct.rmax     = struct.rmax
  outstruct.h        = struct.h

  npairold = 0
  FOR i=0L, nbinuse-1 DO BEGIN 

      ;; calculate area and density of background galaxies
      R1 = outstruct.rmin + i*outstruct.binsize
      R2 = outstruct.rmin + (i+1)*outstruct.binsize

      outstruct.area[i] = !pi*(R2^2 - R1^2)

      outstruct.npair[i] = struct.npair[i]
      outstruct.density[i] = outstruct.npair[i]/outstruct.area[i]/outstruct.nlenses

      outstruct.rmax_act[i] = struct.rmax_act[i]

      outstruct.meanr[i]    = struct.rsum[i]/struct.npair[i]

      ;; Now totals within each rmax_act[i]
      outstruct.tnpair[i] = total(struct.npair[0:i])

  ENDFOR 

  return
END 
