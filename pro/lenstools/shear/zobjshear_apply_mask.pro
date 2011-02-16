PRO zobjshear_edge_intersect, lambda, eta, angmax, edgeflag

  COMMON mask_block, seed, $
    FIRSTLAMBDA, LASTLAMBDA, ETAMIN, ETAMAX, $
    etarange, mask, use_etarange_edgecut, special

  ;; find the distance to the top and bottom edges, as well as the
  ;; distance to the ends.
  ;; Make an edge cut on the ends.  
  ;; for the edges, cut out objects that intersect both edges. 
  ;; Flag objects as hitting 
  ;; 0: no edges
  ;; 2b^0: FIRSTLAMBDA edge
  ;; 2b^1: LASTLAMBDA edge
  ;; 2b^2: the upper edge (maxeta)
  ;; 2b^3: the lower edge (mineta)

  nn = n_elements(lambda)
  edgeflag = bytarr(nn)

  IF use_etarange_edgecut THEN BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; this is for all single stripes
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      mineta = interpol(etarange.mineta, etarange.clambda, lambda)
      maxeta = interpol(etarange.maxeta, etarange.clambda, lambda)

  ENDIF ELSE BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; for multiple stripes use primary bound
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      mineta = replicate(ETAMIN, nn)
      maxeta = replicate(ETAMAX, nn)

  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; search bounds for each lens
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  lens_minlam = lambda - angmax
  lens_maxlam = lambda + angmax

  lens_mineta = eta-angmax
  lens_maxeta = eta+angmax

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; do we cross the bounds of the source catalog?
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  wminlam = where(lens_minlam LT FIRSTLAMBDA, nminlam)
  wmaxlam = where(lens_maxlam GT LASTLAMBDA, nmaxlam)

  wmineta = where(lens_mineta LT mineta, nmineta)
  wmaxeta = where(lens_maxeta GT maxeta, nmaxeta)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set the appropriate flags
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF nminlam NE 0 THEN edgeflag[wminlam] = edgeflag[wminlam] + !MINLAMFLAG
  IF nmaxlam NE 0 THEN edgeflag[wmaxlam] = edgeflag[wmaxlam] + !MAXLAMFLAG

  IF nmineta NE 0 THEN edgeflag[wmineta] = edgeflag[wmineta] + !MINETAFLAG
  IF nmaxeta NE 0 THEN edgeflag[wmaxeta] = edgeflag[wmaxeta] + !MAXETAFLAG

END 


PRO zobjshear_apply_edgecut, edgeflag, bad, good, etaor=etaor, $
                             noetacut=noetacut

  ;; always cut on lambda edges, but three different modes for
  ;; cutting on eta edge

  nn = n_elements(edgeflag)
  badcheck = bytarr(nn)
  IF bad[0] NE -1 THEN badcheck[bad] = 1b

  IF keyword_set(noetacut) THEN BEGIN

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; only cut on lambda edges
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      bad = where( ( $
                     ((edgeflag AND !MINLAMFLAG) NE 0b) OR $
                     ((edgeflag AND !MAXLAMFLAG) NE 0b) $
                   ) OR $
                   (badcheck EQ 1b), nbad, comp=good, ncomp=ngood)

  ENDIF ELSE IF keyword_set(etaor) THEN BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; cut if intersect one of the eta bounds or
      ;; one of the lambda bounds
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      bad = where( ((edgeflag AND !MINLAMFLAG) NE 0b) OR  $
                   ((edgeflag AND !MAXLAMFLAG) NE 0b) OR  $
                   ((edgeflag AND !MINETAFLAG) NE 0b) OR $
                   ((edgeflag AND !MAXETAFLAG) NE 0b) OR $
                   (badcheck EQ 1b), $
                   nbad, comp=good, ncomp=ngood)

  ENDIF ELSE BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; cut only if intersect BOTH of the eta bounds
      ;; or one of the lambda bounds
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;      help, where( ( $
;                     ((edgeflag AND !MINLAMFLAG) NE 0b) OR $
;                     ((edgeflag AND !MAXLAMFLAG) NE 0b)  $
;                   ) )
;      help, where( ( $
;                     ((edgeflag AND !MINETAFLAG) NE 0b) AND $
;                     ((edgeflag AND !MAXETAFLAG) NE 0b) $
;                   ) )

      bad = where( ((edgeflag AND !MINLAMFLAG) NE 0b) OR $
                   ((edgeflag AND !MAXLAMFLAG) NE 0b) OR $
                   ( $
                     ((edgeflag AND !MINETAFLAG) NE 0b) AND $
                     ((edgeflag AND !MAXETAFLAG) NE 0b) $
                   ) OR $
                   (badcheck EQ 1b), nbad, comp=good, ncomp=ngood)
                   
  ENDELSE 


END 
PRO zobjshear_apply_edgecut_old, lambda, eta, angmax, bad, good

  COMMON mask_block, seed, $
    FIRSTLAMBDA, LASTLAMBDA, ETAMIN, ETAMAX, $
    etarange, mask, use_etarange_edgecut, special

  ;; do nothing
  IF good[0] EQ -1 THEN return

  defval = 1000d

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; this makes indexing easier
  ;; Below, good and bad are absolute indices
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  nn = n_elements(lambda)
  maxeta = replicate(defval, nn)
  mineta = maxeta

  IF use_etarange_edgecut THEN BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; this is for all single stripes
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      maxeta[good] = $
        interpol(etarange.maxeta, etarange.clambda, lambda[good])
      mineta[good] = $
        interpol(etarange.mineta, etarange.clambda, lambda[good])

  ENDIF ELSE BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; for multiple stripes use primary bound
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      maxeta[good] = ETAMAX
      mineta[good] = ETAMIN

  ENDELSE 

  ;; too close to the edge or did not pass masks?
  bad=where( (eta+angmax GT maxeta) OR $
             (eta-angmax LT mineta) OR $
             (lambda+angmax GT LASTLAMBDA) OR $
             (lambda-angmax LT FIRSTLAMBDA) OR $
             (maxeta EQ defval), nbad, complement=good, ncomp=ngood)

END 

PRO zobjshear_apply_mask, lambda, eta, angmax, bad, good, edgeflag, $
                          etaor=etaor, noetacut=noetacut

  COMMON mask_block, seed, $
    FIRSTLAMBDA, LASTLAMBDA, ETAMIN, ETAMAX, $
    etarange, mask, use_etarange_edgecut, special

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Find intersections with the edges
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  zobjshear_edge_intersect, lambda, eta, angmax, edgeflag

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; first apply the mask
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  apply_mask2d, mask, lambda, eta, bad, good, edgecut=angmax
  IF good[0] EQ -1 THEN return

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; now check the edges
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  zobjshear_apply_edgecut, edgeflag, bad, good, etaor=etaor, noetacut=noetacut
;  zobjshear_apply_edgecut, lambda, eta, angmax, bad, good


END 
