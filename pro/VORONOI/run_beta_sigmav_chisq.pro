PRO run_beta_sigmav_chisq

  FOR wclr=0,4 DO BEGIN 
      
      beta_sigmav_chisq,wclr
      beta_sigmav_chisq,wclr,/noweight
      beta_sigmav_chisq,wclr,/sigonly

      beta_sigmav_chisq_fixsig,wclr
      beta_sigmav_chisq_fixsig,wclr,/noweight
      beta_sigmav_chisq_fixsig,wclr,/sigonly

  ENDFOR 


END 
