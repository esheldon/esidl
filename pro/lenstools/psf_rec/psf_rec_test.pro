pro psf_rec_test,psp,row,col,psf,eigen1,eigen2,eigen3, ecoeff0, ecoeff1, ecoeff2

;; This tests by returning the eigenimages for r-band

;reconstruct the psf from the KL eigenimages
;stored in psp pointer which is read in with
;read_psf from the psField file
;row,col should be a (5) array such as from cat(1000).rowc

  if n_params() eq 0 then begin
      print,'-syntax psf_rec_test,psp,row,col,psf,eigen1,eigen2,eigen3'
      return
  endif

  rcs=.001
  
  band = 2
  nrow_b=((*psp[band]).nrow_b)[0]
  ncol_b=((*psp[band]).ncol_b)[0]
                                ;assumes they are the same for each eigen
                                ;so only use the 0 one
  rnrow=((*psp[band]).rnrow)[0]
  rncol=((*psp[band]).rncol)[0]
  
  nb=nrow_b*ncol_b
  coeffs=fltarr(nb)
  ecoeff=fltarr(3)
  cmat=(*psp[band]).c
  
  
  for i=0L, nb-1 do coeffs[i]=(row*rcs)^(i mod nrow_b) * (col*rcs)^(i/nrow_b)
  
  for j=0,2 do begin
      for i=0L, nb-1 do begin
          ecoeff[j]=ecoeff[j]+cmat(i/nrow_b,i mod nrow_b,j)*coeffs[i]
      endfor	
  ENDFOR
  
  eigen1 = reform( ((*psp[band]).rrows)[*,0]*ecoeff[0], rncol, rnrow )
  eigen2 = reform( ((*psp[band]).rrows)[*,1]*ecoeff[1], rncol, rnrow )
  eigen3 = reform( ((*psp[band]).rrows)[*,2]*ecoeff[2], rncol, rnrow )

  p=((*psp[band]).rrows)[*,0]*ecoeff[0]+((*psp[band]).rrows)[*,1]*ecoeff[1]+((*psp[band]).rrows)[*,2]*ecoeff[2]

  fac = 200000./total(eigen1)
  eigen1 = eigen1*fac
  eigen2 = eigen2*fac
  eigen3 = eigen3*fac

  fac = 200000./total(p)
  p = p*fac
  psf=reform(p,rncol,rnrow)

  ecoeff0=ecoeff[0]
  ecoeff1=ecoeff[1]
  ecoeff2=ecoeff[2]

  return
end



	
			
				
		


