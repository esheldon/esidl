pro plate_wrapper,platelist=platelist,overwrite=overwrite

   if not keyword_set(platelist) then begin
     spawn,'ls /data1/spectra/1d_15/tsObj/tsObj-????.fit',platelist
   endif
   nfiles=n_elements(platelist)

   time = systime(1)
   outdir = '/data1/spectra/1d_15/tsObj/lens/'

   for i=0,nfiles-1 do begin

       dirsep, platelist[i], indir, platefile
       outfile = outdir + repstr(platefile, ".fit", "_lens.fit")

       tmp = (str_sep(platelist[i], '-'))[1]
       platen = long(  ( str_sep(tmp, ".fit") )[0] )

       ;; don't redo plate unless /overwrite
       IF (NOT fexist(outfile)) AND (NOT keyword_set(overwrite)) THEN BEGIN 
           print
           print,'Reading plate file: ',platelist[i]
           
           collate_plate, pln, platef=platelist(i)

           print
           print,'Writing lens plate file: ',outfile
           
           mwrfits,pln,outdir+outfile,/create
       ENDIF ELSE BEGIN
           print,'Lens plate file: '+outfile+' already exists. Send /overwrite to redo'
       ENDELSE 
   endfor
   ptime,systime(1)-time

   return
   end
