PRO testbug

color_index=2
run=756
rerun=1
camcol=3
field=67
id = 300
nf=40

fetch_dir, run, camcol, rerun, dir, atldir
read_tsobj, dir, pstruct, start=field, nf=nf
fetch_file_list,dir,files,field,nf,run, camcol, rerun, fnums

make_flag_struct,fs
fs.satur='N'
fs.satur_center='N'
fs.bright = 'N'
flag_select,pstruct,fs,color_index,_ss
IF (_ss[0] EQ -1) THEN BEGIN
    GOTO,jump                   ;Avoids embedded ifs
ENDIF 

make_flag_struct, fs
fs.blended='N'
fs.moved='N'
flag_select,pstruct[_ss],fs,color_index,_ss2,/objc
IF (_ss2[0] EQ -1) THEN BEGIN
    _ss = _ss2
    GOTO, jump                  ;Avoids embedded ifs
ENDIF 
_ss = _ss[_ss2]


FOR j=0L, nf-1 DO BEGIN
    w=where(pstruct[_ss].field EQ fnums[j], nw)
    fname=atldir+'fpAtlas-'
    atlas_name,fname,run,camcol,fnums[j]
    print,'Field ',fnums[j]
    FOR i=0L, nw-1 DO BEGIN 
        im=0
        ind = pstruct[ _ss[w[i]] ].id
        get_atlas3, ind, color_index, fname, im, s, row0, col2
    ENDFOR 
ENDFOR 

jump:
return
END 
