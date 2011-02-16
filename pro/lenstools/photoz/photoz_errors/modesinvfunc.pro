PRO modesinvfunc, zLarr, zSarr, sigzsarr, modes

  minzL = 0.02
  maxzL = 0.25
  NzL = 5
  zLarr = arrscl( findgen(NzL), minzL, maxzL ) 

  minzS = 0.2
  maxzS = 0.6
  NzS = 20
  zSarr = arrscl( findgen(Nzs), minzS, maxzS )

  nsigz = 1
  sigzsarr = 0.05

  modes = fltarr(nzL, nzS, nsigz)

  nrand = 50000L

  fac = 1.e4

  FOR iL = 0L, nzL-1 DO BEGIN 
      print,'iL = ',iL
      FOR iS = 0L, nzS-1 DO BEGIN 
          FOR isig=0L, nsigz-1 DO BEGIN 

              mode_sigmacritinv, zLarr[iL], zSarr[iS], sigzsarr[isig], nrand, $
                                 tmode
              modes[iL, iS, isig] = tmode
          ENDFOR 
      ENDFOR 
  ENDFOR 
  

END 
