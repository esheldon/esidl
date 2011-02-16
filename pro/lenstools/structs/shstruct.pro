function shstruct, arrval_input

  if n_elements(arrval_input) eq 0 then begin 
      print,'-Syntax: shstruct = shstruct(arrval)'
      return, -1
  endif 
  
  arrval = double(arrval_input)

  tmparrval = arrval
  tmparrval[*] = 1.e7
  shst = $
    {$
      nlenses: 0L, $            ; total number of lenses used
      nlbin: long(arrval), $    ;number lenses used in each bin
      $
      totpairs: 0LL, $          ;total lens-source pairs
      $
      h: 0d, $                        
      omega_m: 0d, $
      interp_photoz: 0, $
      sigmacrit_style: 0, $
      logbin: 0, $
      nbin: 0, $
      binsize: 0d, $            ; binsize (linear bins)
      rmin: 0d, $
      rmax: 0d, $
      comoving: 0, $
      htm_depth: 0, $
      $
      rmax_act: arrval, $       ; actual max r for each bin
      rmin_act: tmparrval,$
      area_act: arrval, $       ; area from rmax,rmin_act
      $
      ssh: 0d, $                ; shear responsivity
      zmean: 0d, $              ; mean lens redshift
      meanang: arrval, $        ;mean in angular bin
      meanr: arrval, $          ;mean radius
      sigma: arrval, $          ; density contrast
      sigmaerr: arrval, $ 
      sigmaerr2: arrval, $      ;alternative error estimate
      orthosig: arrval, $       ; orthogonal density contrast
      orthosigerr: arrval, $
      orthosigerr2: arrval, $
      scritinv: arrval, $
      npair: long64(arrval), $  ; number of pairs used in each bin
      area: arrval, $           ; area of each bin
      density: arrval, $
      $
      tmeanr: arrval, $         ; These are averaged over interior 
      tsigma: arrval, $         ; bins
      tsigmaerr: arrval, $
      tsigmaerr2: arrval, $
      torthosig: arrval, $
      torthosigerr: arrval, $
      torthosigerr2: arrval, $
      tnpair: long64(arrval), $
      tarea: arrval,$
      $
      rsum: arrval, $
      angsum: arrval, $
      wsum: arrval, $
      wsum2: arrval, $
      wsum_mean: arrval, $
      wsum_err: arrval, $
      owsum: arrval, $
      wscritinvsum: arrval, $
      $
      lenswsum: 0d, $
      wsum_ssh: 0d,$            ;sums for combining later
      sshsum: 0d, $
      zsum: 0d, $
      sigerrsum: arrval, $
      orthosigerrsum: arrval $
    }

  return, shst




end 
