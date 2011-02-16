pro psf_rec,psp,row,col,psf,bands
;reconstruct the psf from the KL eigenimages
;stored in psp pointer which is read in with
;read_psf from the psField file
;row,col should be a (5) array such as from cat(1000).rowc

  if n_params() eq 0 then begin
      print,'-syntax psf_rec,psp,row,col,psf,bands'
      return
  endif

  if n_elements(bands) eq 0 then bands=[0,1,2,3,4]
  nbands=n_elements(bands)
  psf=ptrarr(nbands)
  rcs=.001
  
  for k=0, nbands-1 do begin
      band=bands[k]
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
      
      
      for i=0L, nb-1 do coeffs[i]=(row[band]*rcs)^(i mod nrow_b) * (col[band]*rcs)^(i/nrow_b)
      
      for j=0,2 do begin
          for i=0L, nb-1 do begin
              ecoeff[j]=ecoeff[j]+cmat(i/nrow_b,i mod nrow_b,j)*coeffs[i]
          endfor	
      endfor
      p=((*psp[band]).rrows)[*,0]*ecoeff[0]+((*psp[band]).rrows)[*,1]*ecoeff[1]+((*psp[band]).rrows)[*,2]*ecoeff[2]

      p = p/total(p)*200000.

      psf[k]=ptr_new(reform(p,rncol,rnrow))
  endfor
  
  return
end



	
			
				
		


