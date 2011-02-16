;; copy input struct into a default par_struct
function correlate_pardef

  ;; psample: This is the parameters sample if you like.  This is the
  ;; version of the parameters used. For example, one might use the
  ;; same DPsample, RP, DS, RS but different max radius. Does not change
  ;; with output types.

  ;; Input samples:
  ;;
  ;; DPsample.  This is the input primary sample.  This only changes if
  ;;   the sample selection changes. e.g. different redshift range.
  ;;
  ;; RPsample.  Random primary sample.  Should match the geometry of the
  ;;   DPsample, but redshift dist will be sampled later.
  ;;
  ;; DSsample.  The secondaries. Only changes if selection criteria changes.
  ;; RSsample.  The random secondaries. Only changes if geometry changes.
  ;; 
  ;; Ksample: The k-correction matrix. 
  ;;
  ;; rzsample:  The way the redshift binning is done for histogram matching.
  ;;   Not exactly sure I we versioned this. 

  ;; numerator_output_type: Type of measurements to keep.  Currently 
  ;; supported values:
  ;;   'cl_clr': counts+lum binned by color-lum-rad
  ;;   'cl_r': counts+lum binned by rad
  ;;   'c_r': counts binned by rad
  
  pardef =           $
    {                $
      version:        '0.9', $
      psample:       'NONE', $
      DPsample:      'NONE', $
      RPsample: 'dr4plus01', $
      DSsample: 'dr4plus01', $
      RSsample: 'dr4plus01', $
      Ksample:    'morad01', $
      rzsample:      'rz01', $  ; 40 bins from 0.1 to 0.3
      numerator_output_type: 'cl_clr', $ ; Default clr. The user should set this.
      denominator_output_type: 'c_r', $ ; Default r. Should never change. 
      h:                1.0, $
      omega_m:         0.27, $
      nrad:              18, $
      rmin:            0.02, $
      rmax:            11.5, $
      nlum:              20, $
      lumband:            3, $ ; Band for lum limits/binning and output lum
      loglmin:          9.5, $
      loglmax:         11.7, $
      nkgmr:             20, $
      kgmrmin:          0.0, $
      kgmrmax:          2.0, $
      comoving:           0, $
      depth:              8  $
    }

  return,pardef

end 
