;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; subroutine
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO runclose_match, str, clstr, ned_id,match_struct


  frac = 1.0
  x = [-frac,0,frac,0,-frac]
  y = [0,frac,0,-frac,0]
  usersym,x,y,/fill

  plot,str.ra,str.dec,psym=3,yrange=[min(str.dec),max(str.dec)]
  oplot,clstr.ra,clstr.dec,psym=8,symsize=2.0

  openu,1,'matchnum.dat',/append
  ndstr=create_struct('ned_id',0L,'ned_type', '', 'ned_name', '', 'ned_ra',0.0,'ned_dec',0.0, 'ned_z', 0.0, 'ned_mag', 0.0)
  tol = 40.0 ;; arcseconds
  tol = tol*2.8e-4 ;; degrees
  allow = 1  ;; allow one closest match in the photo structure

  close_match, clstr.ra, clstr.dec, str.ra, str.dec, mclstr, mstr,tol,allow
  IF (mstr[0] EQ -1) THEN nmatches=0 ELSE nmatches=n_elements(mstr)
  printf,1,'   Number of matches:              ',strtrim(string(nmatches),2)
  IF (mstr[0] NE -1) THEN BEGIN
    IF (n_elements(mstr) NE n_elements(mclstr)) THEN BEGIN
      print,''
      print,'SOMETHING IS TERRIBLY WRONG!!'
      help,mclstr
      help,mstr
      print,''
    ENDIF
    s=replicate(ndstr, nmatches)
    IF (nmatches EQ 1) THEN s.ned_id = (ned_id[mclstr])[0] $
    ELSE s.ned_id = ned_id[mclstr]
    s.ned_type = clstr[mclstr].type
    s.ned_name = clstr[mclstr].name
    s.ned_ra = clstr[mclstr].ra
    s.ned_dec = clstr[mclstr].dec
    s.ned_z = clstr[mclstr].z
    s.ned_mag = clstr[mclstr].mag
    combine_structs, s, str[mstr], tmpstr
    IF (n_elements(match_struct) EQ 0) THEN BEGIN 
      match_struct=tmpstr 
    ENDIF ELSE BEGIN
      concat_structs, match_struct, tmpstr, tmp
      match_struct=tmp
    ENDELSE 
    printf,1,'    Number in match_struct: ',n_elements(match_struct)
  ENDIF

close,1
return
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Main procedure
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


PRO nedclstr_match, match_struct, run, ned=ned, beg_col=beg_col,$
                    end_col=end_col

IF n_params() EQ 0 THEN BEGIN
  print,'-Syntax: nedclstr_match, match_struct, run=run, ned=ned, beg_col=beg_col, end_col=end_col'
  return
ENDIF

IF n_elements(beg_col) EQ 0 THEN beg_col=1
IF n_elements(end_col) EQ 0 THEN end_col=6

runstr=strtrim(string(run),2)
dir = run_dir(run)
get_ned, ned, clstr,'GClstr'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Match ned clusters to phils shape catalogs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


tol = 40.0 ;; arcseconds
tol = tol*2.8e-4 ;;degrees
allow = 1  ;; allow one closest match in the photo structure

camcolstr=['1','2','3','4','5','6']
make_corrected_tags, taglist

openw,1,'matchnum.dat'
printf,1,'Matching to run '+runstr
close,1
FOR i=beg_col-1,end_col-1 DO BEGIN
  fname = dir + 'adat'+camcolstr[i]+'c.fit'
  openu,1,'matchnum.dat',/append
  message = ' * Reading from file: '+fname
  print,message
  printf,1,message
  close,1
  fits_info, fname, /silent,n_ext=n_ext

  ;;;;;; step by hundreds through the file  ;;;;;;;
  nstep = n_ext/100
  left = n_ext MOD 100
  FOR jj=0, nstep-1 DO BEGIN     ;;; start with hdu=1 ;;;
    bg=strtrim(string(jj*100+1),2)
    ed=strtrim(string(jj*100+100),2)
    openu,1,'matchnum.dat',/append
    mes=' * Reading hdu '+bg+' to '+ed
    print,mes
    printf,1,mes
    close,1
    read_photo_col,fname,str,start=jj*100+1,nframes=100,$
                   struct_typ='corrected',taglist=taglist
    ;;; match to ned
    runclose_match, str, clstr, w, match_struct
  ENDFOR
  IF (left NE 0) THEN BEGIN
    bg=strtrim(string(nstep*100+1),2)
    ed=strtrim(string(nstep*100+left),2)
    openu,1,'matchnum.dat',/append
    mes = ' *Reading hdu '+bg+' to '+ed
    print,mes
    printf,1,mes
    close,1
    read_photo_col,fname,str,start=nstep*100+1,nframes=left,$
                   struct_typ='corrected',taglist=taglist
    runclose_match, str, clstr, w, match_struct
  ENDIF
  mwrfits,match_struct,dir+run'+runstr+$
                       '-'+camcolstr[i]+'-ned.fit'
;  IF (i EQ 5) THEN BEGIN
;    FOR jj=0,5 DO BEGIN 
;      m=mrdfits(dir+run'+runstr+'-'+$
;                camcolstr[jj]+'-ned.fit', structyp='blah')
;      IF (jj EQ 0) THEN mtot = m ELSE mtot = [mtot,m]
;    ENDFOR
;    mwrfits, mtot, dir+run'+runstr+'-ned.fit'
;  ENDIF
ENDFOR

close,1

spawn,'rm matchnum.dat'

return
END












