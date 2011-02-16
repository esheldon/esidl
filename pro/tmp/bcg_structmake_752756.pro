PRO bcg_structmake_752756

indir='/sdss5/data0/annis/clusters/'

files=findfile(indir+ntostr(1)+'/')
nf=n_elements(files)
ptrlist=ptrarr(nf)
numlist=lonarr(nf)
FOR i=0L,nf-1,1 DO BEGIN
yanny_read,indir+ntostr(1)+'/'+files[i],pdata
ptrlist[i]=ptr_new((*pdata),/no_copy)
numlist=n_elements((*pdata))
yanny_free,pdataa))
yanny_free,pdata
ENDFOR
tstruct=(*ptrlist[0])[0]
zero_struct,tstruct
comb_ptrstruct,tstruct,ptrlist,numlist,outstruct_old







; FOR j=2,6,1 DO BEGIN
;     files=findfile(indir+ntostr(j)+'/')
;     nf=n_elements(files)
;     ptrlist=ptrarr(nf)
;     numlist=lonarr(nf)
;     FOR i=0L,nf-1,1 DO BEGIN
;         yanny_read,indir+ntostr(j)+'/'+files[i],pdata
;         ptrlist[i]=ptr_new((*pdata),/no_copy)
;         numlist=n_elements((*pdata))
;         yanny_free,pdata
;     ENDFOR
;     tstruct=(*ptrlist[0])[0]
;     zero_struct,tstruct
;     comb_ptrstruct,tstruct,ptrlist,numlist,outstruct
;     concat_structs,outstruct_old,outstruct,outstruct_old
;     ptrlist=0
;     outstruct=0
;     numlist=0
; ENDFOR

; mwrfits,'/sdss5/data0/dgrin/wholelists/bcg_hamf_752756/'+'bcg_hamf_752756.fit',outstruct_old,/create,/destroy

return
end
