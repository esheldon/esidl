pro smooth_col, filename, ofilebeg, nframes=nframes, start=start, groupn=groupn, taglist=taglist

;+
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Applies the frame by frame shape smoothing function to an SDSS column....
;  
; Inputs:  filename: name of input fits file
;	   ofilebeg: output fits file (will overwrite existing file)
;          taglist:  A list of photo tags in all CAPS that the user wants
;                    in struct
;          nframes:  Optional parameter which tells how many frames to read
;                    from filename
;          start:    Beginning fram
;	   groupn: How many frames to take together
;
; Outputs: fits file containing smoothed col information
;
; Created: Tim McKay
; Date: 6/10/99
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Help message
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if n_params() eq 0 then begin
   print,'-syntax smooth_col, filename, ofilebeg, nframes=nframes, start=start, groupn=groupn, taglist=taglist'
   return
endif

t=systime(1)

print,' * Reading from file: ',filename
fits_info, filename, /silent, n_ext=framemax
framemax=long(framemax)
print,' * Total Number of Frames: ',framemax

IF NOT keyword_set(nframes) THEN nframes=framemax ELSE nframes=long(nframes)

IF keyword_set(start) THEN BEGIN
  start=long(start) > 1
  print,' * Beginning with frame ',start
  max = long(nframes) + start-1
  sta=start
ENDIF ELSE BEGIN 
  max=nframes
  sta=1
ENDELSE

max = max < framemax
n_ext = max-sta+1


if not keyword_set(groupn) then begin
	groupn = 30
endif


;Figure out how many groupn frame units there are, and what's left over....
nunits=fix(n_ext/groupn)
nextra=n_ext MOD groupn
print,' * Will read ',ntostr(n_ext),' fields.'
print,' * Will use ',ntostr(nunits),' groups ',$
  'and ',ntostr(nextra),' extra files'
print

;Get the taglist for these files....
if not keyword_set(taglist) then begin
	make_lenstags,tlist
endif else begin
	tlist=taglist
endelse

gfile=ofilebeg+'_smooth_gal.fit'
sfile=ofilebeg+'_smooth_star.fit'

;Start reading them in one group at a time
ntotal=0
for i=0L,nunits-1,1 do BEGIN
    print,' * Reading group ',ntostr(i+1),'/',ntostr(nunits)
    print
    l=0.
    read_photo_col,filename,l,start=sta+i*groupn,nframes=groupn,$
      taglist=tlist,maxf=max

    ; Must extract on a frame by frame basis because shapes change over time.
    sss=sta+i*groupn+10

    w=where(l.field EQ sss,nw)
    extract_corrstars,l[w],2,tmp1,/silent
    extract_lensgalaxies,l[w],2,tmp2,/silent, max_mag=22.0
    os=w[tmp1]
    gs=w[tmp2]
    FOR jj=sss+1, sss+groupn-1 DO BEGIN 
      w=where(l.field EQ jj,nw)
      extract_corrstars,l[w],2,tmp1,/silent
      extract_lensgalaxies,l[w],2,tmp2,/silent
      os=[os,w[tmp1]]
      gs=[gs,w[tmp2]]
    ENDFOR 

    print
    print,' * Smoothing Shapes by Field'
    shape_smooth,l[os],sms,/stars
    shape_smooth,l[gs],gms
    if (i eq 0) then begin
       	mwrfits,sms,sfile,/create
       	mwrfits,gms,gfile,/create
    endif else begin
	mwrfits,sms,sfile
	mwrfits,gms,gfile
    endelse
endfor

;Now read what's left over...
sss = sta+nunits*groupn > sta
if (nextra gt 5) then BEGIN
    print,' * Reading What is Left'
    l=0.
    read_photo_col,filename,l,start=sss,nframes=nextra,$
	taglist=tlist,maxf=max

    ; Must extract on a frame by frame basis because shapes change over time.
    sss=sss+10
    w=where(l.field EQ sss,nw)
    extract_corrstars,l[w],2,tmp1,/silent
    extract_lensgalaxies,l[w],2,tmp2,/silent, max_mag=22.0
    os=w[tmp1]
    gs=w[tmp2]
    FOR jj=sss+1, sss+nextra-1 DO BEGIN 
      w=where(l.field EQ jj,nw)
      extract_corrstars,l[w],2,tmp1,/silent
      extract_lensgalaxies,l[w],2,tmp2,/silent
      os=[os,w[tmp1]]
      gs=[gs,w[tmp2]]
    ENDFOR 
;    extract_corrstars,l,2,os
;    extract_lensgalaxies,l,2,gs,/silent
    shape_smooth,l[os],sms,/stars
    shape_smooth,l[gs],gms
    IF (nunits EQ 0) THEN BEGIN 
      mwrfits,sms,sfile,/create
      mwrfits,gms,gfile,/create
    ENDIF ELSE BEGIN 
      mwrfits,sms,sfile
      mwrfits,gms,gfile
    ENDELSE 
    nunits=nunits+1
endif

;Now put them all in ONE structure....
print
print,' * Creating One big structre.'
print
st='smooth2'
for i=0,nunits-1,1 do begin
    s = mrdfits(sfile,i+1,hdr,structyp=st, /silent)
    if (i eq 0) then begin
	sms=s
    endif else begin
	sms=[sms,s]
    endelse
endfor
mwrfits,sms,sfile,/create

for i=0,nunits-1,1 do begin
    g = mrdfits(gfile,i+1,hdr,structyp=st, /silent)
    if (i eq 0) then begin
	gms=g
    endif else begin
	gms=[gms,g]
    endelse
endfor
mwrfits,gms,gfile,/create

ptime,systime(1)-t
end




