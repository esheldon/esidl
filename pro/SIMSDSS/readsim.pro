PRO readsim,cat


  if N_params() eq 0 then begin
     print,'-Syntax: readsim, cat1,cat2....'
     print,''
     print,'Use doc_library,""  for more help.'  
     return
  endif

dir = '/sdss4/data1/esheldon/'
name = dir+'test_ell_3.fit'
fits_info,name,/silent,n_ext=n

j=0
FOR i=1,n DO BEGIN 
  j=j+1
  print,j
  tmp=0
  tmp=mrdfits(name,i,hdr,/silent,structyp='tteemmpp')
  IF i EQ 1 THEN cat=tmp ELSE cat=[cat,tmp]
ENDFOR

mwrfits,cat,'tmp3.fit'
cat=0

name = dir+'test_ell_4.fit'
fits_info,name,/silent,n_ext=n
FOR i=1,n DO BEGIN 
  j=j+1
  print,j
  tmp=0
  tmp=mrdfits(name,i,hdr,/silent,structyp='tteemmpp')
  IF i EQ 1 THEN cat=tmp ELSE cat=[cat,tmp]
ENDFOR
mwrfits,cat,'tmp4.fit'
cat=0


name = dir+'test_ell_5.fit'
fits_info,name,/silent,n_ext=n
FOR i=1,n DO BEGIN 
  j=j+1
  print,j
  tmp=0
  tmp=mrdfits(name,i,hdr,/silent,structyp='tteemmpp')
  IF i EQ 1 THEN cat=tmp ELSE cat=[cat,tmp]
ENDFOR
mwrfits,cat,'tmp5.fit'
cat=0


return
end































