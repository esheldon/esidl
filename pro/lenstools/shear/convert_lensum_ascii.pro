PRO convert_lensum_ascii,name,arrval,outstruct,$
                         SKIPLINE = skipline, NUMLINE = numline,DOUBLE=double, $
                         ENDSKIP = endskip

  if N_params() lt 2 then begin
     print,'Syntax - RDFLOATSTR, name, arrval, outstruct, '
     print,'                    /DOUBLE, SKIPLINE =, NUMLINE = ]'
     return
  endif

  nskip = 0

; Get number of lines in file

   get_lun,lun

   nlines = NUMLINES( name )
   if nlines LT 0 then return

   if not keyword_set( SKIPLINE ) then skipline = 0
   nlines = nlines - skipline
   IF n_elements(endskip) NE 0 THEN nlines = nlines-endskip
   if keyword_set( NUMLINE) then nlines = numline < nlines

;Read first line, and determine number of columns of data

   openr, lun, name
   temp = ''
   if skipline GT 0 then $
        for i=0,skipline-1 do readf, lun, temp
   readf,lun,temp
   colval = str_sep( strtrim( strcompress(temp),2),' ')
   ncol = N_elements(colval)

;Create big output array and read entire file into the array

   if keyword_set(DOUBLE) then bigarr = dblarr(ncol, nlines, /NOZERO) $ 
                          else bigarr = fltarr(ncol, nlines, /NOZERO) 

   close,lun
   openr, lun, name
   if skipline GT 0 then $
        for i=0,skipline-1 do readf, lun, temp

   ptime,t,/savetime
   print,'Reading ',ntostr(nlines),' lines from file ',name
   readf, lun, bigarr
   free_lun, lun

   print, ntostr(nlines) + ' lines of data read'
   print,'Copying into structure'
   print
   ;Nvector = (N_params()-1) < ncol

   nbin = n_elements(arrval)
   use_struct = zlensumstruct(arrval)
   use_struct = create_struct('z', 0., 'clambda',0d, 'ceta',0d, use_struct)

   outstruct = replicate(use_struct, nlines)

   outstruct.index =          reform( bigarr[0,*] )
   outstruct.zindex =         reform( bigarr[1,*] )

   outstruct.z =              reform( bigarr[2,*] )
   outstruct.clambda =        reform( bigarr[3,*] )
   outstruct.ceta =           reform( bigarr[4,*] )

   outstruct.pixelmaskflags = reform( bigarr[5,*] )
   outstruct.scritinv =       reform( bigarr[6,*] )
   outstruct.totpairs =       reform( bigarr[7,*] )
   outstruct.sshsum =         reform( bigarr[8,*] )
   outstruct.wsum_ssh =       reform( bigarr[9,*] )
   outstruct.weight =         reform( bigarr[10,*] )
   outstruct.ie =             reform( bigarr[11,*] )

   offset = 12
   FOR i=0L, nbin-1 DO BEGIN 
       outstruct.npair[i] =          reform( bigarr[0*nbin  + offset + i, *] )
       outstruct.rmax_act[i] =       reform( bigarr[1*nbin  + offset + i, *] )
       outstruct.rmin_act[i] =       reform( bigarr[2*nbin  + offset + i, *] )
       outstruct.rsum[i] =           reform( bigarr[3*nbin  + offset + i, *] )
       outstruct.sigma[i] =          reform( bigarr[4*nbin  + offset + i, *] )
       outstruct.sigmaerr[i] =       reform( bigarr[5*nbin  + offset + i, *] )
       outstruct.orthosig[i] =       reform( bigarr[6*nbin  + offset + i, *] )
       outstruct.orthosigerr[i] =    reform( bigarr[7*nbin  + offset + i, *] )
       outstruct.sigerrsum[i] =      reform( bigarr[8*nbin  + offset + i, *] )
       outstruct.orthosigerrsum[i] = reform( bigarr[9*nbin  + offset + i, *] )
       outstruct.wsum[i] =           reform( bigarr[10*nbin + offset + i, *] )
       outstruct.owsum[i] =          reform( bigarr[11*nbin + offset + i, *] )
   ENDFOR 


   return
END 
