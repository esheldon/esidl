PRO zobjshear_looplens_print_status, group, nstep, sigsum, sigwsum

  IF sigwsum GT 0.0 THEN BEGIN 
      mean_denscont = ntostr( sigsum/sigwsum )
      err_denscont = ntostr( sqrt(1./sigwsum) )
  ENDIF ELSE BEGIN 
      mean_denscont = '0'
      err_denscont = '0'
  ENDELSE 

  print,'Lens Group = ',ntostr(group+1)+'/'+ntostr(nstep),' Mean dens. cont. = ' + $
        mean_denscont+' '+!plusminus+' '+err_denscont

return

  IF wsum GT 0. THEN BEGIN 
      mean_shear_p = ntostr( etansum/wsum/2. )
      mean_shear_err_p = ntostr( sqrt( etanerrsum )/wsum/2. )
  ENDIF ELSE BEGIN 
      mean_shear_p = '0'
      mean_shear_err_p = '0'
  ENDELSE 
;  echo,['Lens Group = ',$
;        ntostr(group+1)+'/'+ntostr(nstep),$
;        '  Mean shear = ',$
;        mean_shear_p+' '+!plusminus+' '+mean_shear_err_p], $
;       color=['cyan','none','cyan','none'],$
;       bold=[1,0,1,0],$
;       nonewline=[1,1,1,0]

  print,'Lens Group = ',ntostr(group+1)+'/'+ntostr(nstep),' Mean shear = ' + $
        mean_shear_p+' '+!plusminus+' '+mean_shear_err_p

END 
