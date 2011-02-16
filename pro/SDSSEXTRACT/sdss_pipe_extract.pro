pro sdss_extract,im,cat, fwhm=fwhm, tmpdir=tmpdir, imfile=imfile, objfile=objfile, parfile=parfile, flagstr=flagstr, remove=remove

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;A wrapper for sdss_pipe to run on one image
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if n_params() eq 0 then begin
	print,'-syntax  sdss_extract, im, cat, fwhm=fwhm, tmpdir=tmpdir, imfile=imfile, objfile=objfile, parfile=parfile, flagstr=flagst'
	return
endif

IF (NOT keyword_set(imfile)) THEN imfile='test.fits'
IF (NOT keyword_set(objfile)) THEN objfile = 'test_sobj.fits'
IF (NOT keyword_set(tmpdir)) THEN tmpdir='/tmp/'
IF (NOT keyword_set(remove)) THEN remove=0 ELSE remove=1

writefits,tmpdir+imfile,im

if n_elements(parfile) eq 0 then begin
	parfile='~/idl.lib/sdss.par'
endif

IF (NOT keyword_set(fwhm)) THEN fwhm = 0

sdss_pipe, imfile, tmpdir, parfile, whyflag=whyflag, fwhm=fwhm
cat=mrdfits(tmpdir+objfile,1,hdr,/silent)

;;;; Remove the files.  If you want to check them then alter this

IF remove THEN BEGIN 
    print,'Removing ',tmpdir+imfile,' and ',tmpdir+objfile
    command = 'rm '+tmpdir+imfile+' '+tmpdir+objfile+' '+tmpdir+'check.fits'
    spawn,command
ENDIF 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  Check the flags
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

IF keyword_set(flagstr) THEN BEGIN 
  fst = create_struct('id',0, 'whyflag',0.0, 'flag','')
  w=where(whyflag NE 0, nw)
  IF w[0] NE -1 THEN BEGIN
    badflags = whyflag[w]
    flags = ['negsum','bigshift','negsum','negm','lowsum',$
             'negdegn','negweight','maxit','baddet']
    
    flagstr = replicate(fst, nw)
    flagstr.id = w
    flagstr.whyflag = whyflag[w]
    flagstr.flag = flags[ whyflag[w]-1 ]


    g=where(badflags NE 8, ng)  ;; 8 = Maxit Reached 
    print
    print,strtrim(string(nw-ng),2),' With Max Iterations Reached'

    IF ng NE 0 THEN BEGIN 
      FOR i=0, ng-1 DO $
        print,'Id  ',strtrim(string(w[g[i]]),2),'  ',flags[ badflags[g[i]]-1 ]
    ENDIF

  ENDIF
ENDIF


return
end



