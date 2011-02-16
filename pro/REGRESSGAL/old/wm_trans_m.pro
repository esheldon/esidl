pro wm_trans_m,matrixfile,evectfile,ata,btb,nbins,num
;reads in the matrix file
;and computes the transpose(matrix)*matrix
;for the two matrices a and b
;a square matrix for the regression method
;of galaxy-galaxy lensing
; Now does weighting E.S.S.

vint = 0.32^2

if n_params() eq 0 then begin
	print,'-syntax  wm_trans_m,matrixfile,ata,btb,nbins,num'
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

if n_elements(num) eq 0 then begin
	num=numlines(matrixfile)
endif

get_lun,munit
openr,munit,matrixfile
openr,eunit,evectfile,/get_lun

print,'reading ',matrixfile
i=0L
j=0L
a=0.0
b=0.0

ata=fltarr(nbins,nbins)
btb=ata
arow=fltarr(nbins)
brow=arow
ind=intarr(nbins)


readf, eunit, n0, e1, e2, uncert
weight = 1./(vint + uncert^2)

for n=0L, num-1 do begin
	if n mod 10000 eq 0 then print,n
	readf,munit,i,j,a,b
	if i ne n0 then begin
		;a new i has begun do
		;do the computation and start a new one
		w=where(ind,nnz)
		for v=0,nnz-1 do begin
			for u=0,nnz-1 do begin
				ata(w(u),w(v))=ata(w(u),w(v))+$
				   arow(w(u))*arow(w(v))*weight
				
				btb(w(u),w(v))=btb(w(u),w(v))+$
				   brow(w(u))*brow(w(v))*weight	

			endfor
		endfor
                n0=i
                readf, eunit, ii, e1, e2, uncert
                weight = 1./(vint + uncert^2)
		arow(w)=0.0
		brow(w)=0.0
		ind(w)=0
	endif

	arow(j)=a
	brow(j)=b
	ind(j)=1
endfor	

close,munit
free_lun,munit

return
end



 
