PRO zobjshear_setup_masks, stripe, flam, llam, $
                           maskdir=maskdir, etarangedir=etarangedir

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: zobjshear_setup_masks, stripe, FIRSTLAMBDA, LASTLAMBDA, maskdir=maskdir, etarangedir=etarangedir'
      return
  ENDIF 

  COMMON mask_block, seed, $
    FIRSTLAMBDA, LASTLAMBDA, ETAMIN, ETAMAX, $
    etarange, mask, use_etarange_edgecut, special

  nstripe = n_elements(stripe)

  FIRSTLAMBDA = flam
  LASTLAMBDA  = llam

  read_stripe_mask, stripe, mask, indir=maskdir

  primary_bound_multi, stripe, bound
  ETAMIN = min(bound.etamin)
  ETAMAX = max(bound.etamax)

  IF nstripe EQ 1 THEN BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; always cut relative to etarange for
      ;; individual stripes
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      use_etarange_edgecut=1
      read_etarange, stripe, etarange, indir=etarangedir

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; these are special stripes. Spectroscopy is not a rectangle
      ;; in lambda, eta so to generate random points need etarange
      ;; and generate uniformly in lambda rather than sin(clambda)
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      IF (stripe EQ 82 OR stripe EQ 10 OR $
          stripe EQ 86 OR stripe EQ 76) THEN BEGIN 
          special = 1
      ENDIF ELSE BEGIN 
          special = 0
      ENDELSE 
  ENDIF ELSE BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; multiple stripes we use edgecut relative to the
      ;; primary bounds.  No special: generate in sin(clambda)
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      use_etarange_edgecut = 0
      etarange = 0
      special = 0

  ENDELSE 

END 
