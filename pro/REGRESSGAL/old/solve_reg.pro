pro solve_reg,matrixfile,evectfile,nbins,ata,btb,shear1,shear2,num
;reads in the matrix file and evectfile
;and computes the solution to the regression problem
;of galaxy-galaxy lensing

if n_params() eq 0 then begin
	print,'-syntax solve_reg,matrixfile,evectfile,nbins,ata,btb'
	print,'shear1,shear2,num'
	return
endif

if exist(matrixfile) eq 0 then begin
	print,'file: ',matrixfile,' not found'
	return
endif
if exist(evectfile) eq 0 then begin
	print,'file: ',evectfile,' not found'
	return
endif

a_inv=invert(ata)
b_inv=invert(btb)

if n_elements(num) eq 0 then begin
	num=numlines(matrixfile)
endif
	
get_lun,munit
openr,munit,matrixfile
get_lun,eunit
openr,eunit,evectfile

i=0L
j=0L
a=0.0
b=0.0

arow=fltarr(nbins)
brow=arow
ind=intarr(nbins)
shear1=fltarr(nbins)
shear2=shear1

n0=0
readf,eunit,ii,e1,e2

for n=0L, num-1 do begin
	if n mod 10000 eq 0 then print,n,i,ii,j
	
	readf,munit,i,j,a,b
	if i ne n0 then begin		
		;a new i has begun do
		;do the computation and start a new one
		
		shear1=shear1+e1*(a_inv##arow)
		shear2=shear2+e2*(b_inv##brow)	
		w=where(ind)		
		n0=i
		arow(w)=0.0
		brow(w)=0.0
		ind(w)=0
		readf,eunit,ii,e1,e2
		;now read the next ellipticity
	endif

	arow(j)=a
	brow(j)=b
	ind(j)=1
endfor	

close,munit
free_lun,munit
close,eunit
free_lun,eunit

return
end



 
