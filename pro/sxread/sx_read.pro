pro sx_read,file,str,names=names,usenames=usenames
;read an sx_output which is in ascii .tbl format
;will read into an IDL structure with the right data types
;unless you do usenames=0, it will use the names from sx, tranlated in 
;a way like "model_counts_2" 
;it doesn't yet convert to tsObj names 
;if you do usename=0 it will use names tag0,tag1,tag2...
;if you give it "names" it will use those names
;this is a beta version but seems quite robust so far

if n_params() eq 0 then begin
	print,'-syntax sx_read,file,str'
	return
endif

if n_elements(usenames) eq 0 then usenames=1 

num=numlines(file)

print,'Reading file: ',file
print,ntostr(num),' lines'

line=''
get_lun,unit
openr,unit,file
i=0L
count=0L
perlast=-1

;; ESS
back = string(8B)
pstr = '00'
perlast = 0

while EOF(unit) ne 1 do begin
	readf,unit,line
	beg=strmid(line,0,1)
	if beg eq '#' then begin
		;print,line
		if i eq 1 then begin
			lin=strmid(line,1,5)
			nf=long(lin)
			print,nf,'  tags'
                        ;;print,pstr+'%',format='(a,$)'
		endif		
		if i eq 2 and usenames eq 1 then begin
			lin=strmid(line,2,strlen(line))
			tagnames=str_sep(lin,' ')
			for kk=0, nf-1 do begin
				;help,tagnames(kk)
				kkname=tagnames(kk)
				kks=str_sep(kkname,'.')
				nkk=n_elements(kks)
				if nkk gt 1 then begin
					kkname='"'+kks(nkk-1)
				endif
				pp=strpos(kkname,'[')
				if pp ne -1 then begin
					kkname=strmid(kkname,0,pp)+'_'+strmid(kkname,pp+1,1)+'"'
					tagnames(kk)=kkname
				endif
				pp=strpos(kkname,'(')
				if pp ne -1 then begin
					kkname=strmid(kkname,0,pp)+'"'
					tagnames(kk)=kkname
				endif
				kkname=strcompress(kkname,/rem)
				kkname=strmid(kkname,1,strlen(kkname)-2)
				tagnames(kk)=kkname
				;;print,tagnames(kk)
				;wait,.4
			endfor
			names=tagnames
		endif	
		if i eq 3 then begin
			lin=strmid(line,2,strlen(line))
			types=str_sep(lin,' ')
			j=0
			if j lt n_elements(names) then tname=names(j) else tname='tag'+strcompress(j,/rem)
			if types(j) eq 'int' then val=0
			if types(j) eq 'float' then val=0.0
			if types(j) eq 'double' then val=0d
			str=create_struct(tname,val)
			for j=1, nf-1 do begin
				if j lt n_elements(names) then tname=names(j) else tname='tag'+strcompress(j,/rem)
				if types(j) eq 'int' then val=0
				if types(j) eq 'float' then val=0.0
				if types(j) eq 'double' then val=0d
				if types(j) eq 'uint' then val=0L
				if types(j) eq 'string' then val=''
				str=create_struct(str,tname,val)
			endfor	
			str=replicate(str,num)
                        help,str,/str,output=strhelp
                        colprint,strhelp
                        print
		endif		

		
	endif else begin
		lin=strtrim(strcompress(line),2)
		vals=str_sep(lin,' ')
		if n_elements(vals) ne nf then begin
			print,'WARNING ',' line ',i
			print,"doesn't have ",nf," values"
			print,'skipping'
		endif else begin
			w=where(vals eq 'NoLink',wif)
			if wif gt 0 then vals(w)='0'	
			for k=0, nf-1 do str(count).(k)=vals(k)
			count=count+1
			
		ENDELSE

                per=fix(100.0*i/float(num))
                if per gt perlast then begin
                    perlast=per
                    pstr = strmid( strtrim(per,2), 0, 10)
                    IF strlen(pstr) EQ 1 THEN pstr='0'+pstr
                    print,back+back+back+pstr+'%',format='(a,$)'
;		PRINT, FORMAT = '(3x,1(i3.2,a3,3x))', per, '% '+STRING(BYTE(141))
                endif

	endelse

	i=i+1		
endwhile

print,back+back+back+'100'+'%'

str=str(0:count-1)
close,unit
free_lun,unit
return
end
