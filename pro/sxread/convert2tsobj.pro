pro convert2tsobj,names,types,str,nnums,fail=fail
;convert sx tagnames to tsObj tagnames
;and make a structure

if n_params() eq 0 then begin
	print,'-syntax convert2tsobj,names,types,tsobj'
	return
endif

fail=0
s=sort(names)
	;asciibetical sort
names=names(s)
types=types(s)
n=n_elements(names)

n=n_elements(names)

oldbase=''
nt=0
nnums=intarr(n)
ntypes=strarr(n)
nnames=strarr(n)
nindex=intarr(n)

for i=0, n-1 do begin
	t=names(i)
	base=t	
	ss=str_sep(t,'_')
	nss=n_elements(ss)
	ext=-1
	if nss gt 1 and strlen(ss(nss-1)) eq 1 then begin
		last=ss(nss-1)
		blast=byte(last)
		blast=blast(0)
		if blast gt 47 and blast lt 58 then begin
			;we have a number between 0 and 9
			base=strmid(t,0,strlen(t)-2)
                        IF strupcase(base) EQ 'MODELCOUNTS' THEN base = 'COUNTS_MODEL'
			ext=fix(ss(nss-1))		
		endif
	endif
	if base ne oldbase then begin
		;we have a new one
                
		nnames(nt)=base	
		ntypes(nt)=types(i)
		nnums(nt)=1
		nindex(nt)=s(i)
		nt=nt+1
		oldbase=base
	endif else begin
		;same one just increment
		nnums(nt-1)=nnums(nt-1)+1
	endelse
endfor

nindex=nindex(0:nt-1)
us=sort(nindex)
nindex=nindex(us)
	;truncate and unsort them
nnums=nnums(0:nt-1)
nnums=nnums(us)
ntypes=ntypes(0:nt-1)
ntypes=ntypes(us)
nnames=nnames(0:nt-1)
nnames=nnames(us)

for i=0, nt-1 do begin
	case ntypes(i) OF
                'short':        val=0
		'int':          val=0L
                'long':         val=0LL
    		'float':        val=0.0
       	        'double':       val=0d
                'uint':         val=0UL
                'uint64':       val=0ULL
                'string':       val=''
        else: begin
                print,'illegal data type'
                print,types(i)
                print,nnames(i),ntypes(i)
                fail=1
		dewdewwe
               	return
        	end
        endcase

	;translate_tsobj,nnames(i),tsname
	;this function doesn't exist yet
	tsname=nnames(i)
	nnames(i)=tsname
	if i eq 0 then begin
		str=create_struct(nnames(i),val)
	endif else begin
		if nnums(i) gt 1 then val=replicate(val,nnums(i))
		str=create_struct(str,nnames(i),val)
	endelse
endfor

return
end


