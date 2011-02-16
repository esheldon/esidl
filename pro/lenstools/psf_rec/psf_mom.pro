pro psf_mom,psp,row,col,sig,ixx,iyy,ixy

;take the psp read in from the psField file
;and calculate the seconds moments 
;sig^2 of the reconstructed
;psf at each row,col point
;takes advantage of the linearity of moments
;so only has to calculate the moments of each
;template once per field ,not at every point
;row,col should be (5,num) arrays such as from cat.rowc

  ;; !! WILL BE SENDING EIGENIMAGES TO ADAPTIVE MOMENT CODE
  ;; MUST REWRITE. ALSO, ONLY NEED TO DO 3 BANDS

  if n_params() eq 0 then begin
      print,'-syntax psf_mom,psp,row,col,sig'
      return
  endif
 
;  setup_mystuff
;  bindir='/sdss3/usrdevel/philf/idl.lib/'
;  bindir = !mybindir
 
  ss=size(row)
  if ss[0] ne 2 then begin
      print,'row,col should be (5,num) arrays such as from cat.rowc'
      return
  endif

  bands=[0,1,2,3,4]             ;do all five bands
  nbands=5
  num=ss[2]                     ;number of objects
  rcs=.001	
;ns=11
  ns=51                         ;this will truncate the psf to size (ns,ns)
                                ;before measuring moments
  ns2=(ns-1)/2
  sig=fltarr(num,nbands)
  ixx=sig                       ;the sigma's that are output
  iyy=sig
  ixy=sig

  l=lindgen(ns*ns)              ;this part makes the weight * r^2 to 
  x=l mod ns                    ;integrate over the psf
  y=l /ns			;to change the weight then change this part
  x=x-ns2
  y=y-ns2
  x=double(reform(x,ns,ns))
  y=double(reform(y,ns,ns))
  xx=x^2
  yy=y^2
  xy=x*y
  rr2=xx+yy
;r2=rr2
  weightmom=exp(-.5*(rr2/18.0))
  wr2=rr2*weightmom
  wx2=xx*weightmom
  wy2=yy*weightmom
  wxy=xy*weightmom
  
  n_moms = 3
  moms=fltarr(3,n_moms)		;for each eigentemplate
  tmoms=fltarr(3)               ;to generalize could have mom=fltarr(3,n_moms)
  sums=fltarr(3)                             ;ie measure n_moms different moments
                                ;then just add a for loop inside the innermost for loop

 

  for k=0, nbands-1 do begin
      band=bands[k]
      nrow_b=((*psp[band]).nrow_b)[0]
      ncol_b=((*psp[band]).ncol_b)[0]
                                ;assumes they are the same for each eigen
                                ;so only use the 0 one , I believe this will always be true
      rnrow=((*psp[band]).rnrow)[0]
      rncol=((*psp[band]).rncol)[0] ;Size of eigentemplates. Size of rrows = rnrow*rncol
      numh=(rnrow-1)/2
      
      nb=nrow_b*ncol_b
      coeffs=fltarr(num,nb)
                                ;stores the row,col dependence of the coeeficients
      ecoeff=fltarr(num,3)
                                ;stores the three coeeficents of the three eigentemplates
      cmat=(*psp(band)).c
      

      for i=0L, nb-1 do coeffs[*,i]=(row[band,*]*rcs)^(i mod nrow_b) * (col[band,*]*rcs)^(i/nrow_b)
                                ;this is done quickly in IDL 	
      for j=0,2 do begin
          for i=0L, nb-1 do begin
              ecoeff[*,j]=ecoeff[*,j]+cmat[i/nrow_b,i mod nrow_b,j]*coeffs[*,i]
          endfor  
          eigtemp=  reform(  (*psp[band])[j].rrows , rnrow, rnrow)
                                ;read in the eigentemplate
                                ;assume square
          eigtemp=eigtemp[numh-ns2:numh+ns2, numh-ns2:numh+ns2]	
                                ;trim the size

          tmoms[j]=total(wr2*eigtemp)
          moms[0,j] = total(wx2*eigtemp)
          moms[1,j] = total(wy2*eigtemp)
          moms[2,j] = total(wxy*eigtemp)
          sums[j]=total(eigtemp*weightmom)	

          ;; default inputs for ad_momi
          sky=[0.]
          skysig = [(*psp[5]).skysig[band]]
          
;          shiftmax=10.
;          n=1L
;          nx=51L
;          ny=51L
;          mag=[0.]
;          cenx=[float(numh)]
;          ceny=[float(numh)]
;          sixx = [4.5]
;          siyy = [4.5]
;          sixy = [0.]
;          rho4=[0.]
;          uncert=[0.]
;          numiter=[0L]
;          wcenx=[0.]
;          wceny=[0.]
;          whyflag=[0L]
;          ad_momi2, eigtemp, cenx, ceny, shiftmax, sky, skysig, sixx, siyy, sixy, rho4, uncert, $
;            wcenx, wceny, numiter, whyflag, isum
;          ff=call_external(bindir+'ad_momi.so','ad_mom_', $
;                           cenx, ceny, $
;                           sixx, siyy, sixy, $
;                           n, shiftmax, $
;                           eigtemp, sky, nx, ny, uncert, skysig, mag,$
;                           numiter, wcenx, wceny, whyflag, rho4)
;          print,whyflag[0],numiter[0]
;          eixx[j] = sixx*isum
;          eiyy[j] = siyy*isum
;          eixy[j] = sixy*isum
;          eixxiyy[j] = eixx[j]+eiyy[j]
;          esums[j] = isum
                                ;this is the moment and sum for each eigentemplate
      endfor
      
                                ;now for the time saver part, we only need sums and moms (six numbers)
                                ;and ecoeff array to compute the moments anywhere
;      colprint,eixx,eiyy,eixy,eixxiyy,esums
      sig[*,band]=(ecoeff[*,0]*tmoms[0]+ecoeff[*,1]*tmoms[1]+ecoeff[*,2]*tmoms[2])  / $
        (ecoeff[*,0]*sums[0]+ecoeff[*,1]*sums[1]+ecoeff[*,2]*sums[2])
      ixx[*,band]=(ecoeff[*,0]*moms[0,0]+ecoeff[*,1]*moms[0,1]+ecoeff[*,2]*moms[0,2])  / $
        (ecoeff[*,0]*sums[0]+ecoeff[*,1]*sums[1]+ecoeff[*,2]*sums[2])
      iyy[*,band]=(ecoeff[*,0]*moms[1,0]+ecoeff[*,1]*moms[1,1]+ecoeff[*,2]*moms[1,2])  / $
        (ecoeff[*,0]*sums[0]+ecoeff[*,1]*sums[1]+ecoeff[*,2]*sums[2])
      ixy[*,band]=(ecoeff[*,0]*moms[2,0]+ecoeff[*,1]*moms[2,1]+ecoeff[*,2]*moms[2,2])  / $
        (ecoeff[*,0]*sums[0]+ecoeff[*,1]*sums[1]+ecoeff[*,2]*sums[2])
;      mysig[*,band] = (ecoeff[*,0]*eixxiyy[0]+ecoeff[*,1]*eixxiyy[1]+ecoeff[*,2]*eixxiyy[2])  / $
;        (ecoeff[*,0]*esums[0]+ecoeff[*,1]*esums[1]+ecoeff[*,2]*esums[2])
      
  endfor	
  
  sig=sqrt(sig > 0.0)
;  mysig = sqrt(mysig > 0.)
                                ;return the sigma not variance
  return
end















