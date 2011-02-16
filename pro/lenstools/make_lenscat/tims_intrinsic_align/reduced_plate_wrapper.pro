pro reduced_plate_wrapper,ptot,platelist=platelist,keepnomatch=keepnomatch

  if not keyword_set(platelist) then begin
      spawn,'ls /net/cheops1/data1/spectra/1d_15/tsObj/lens/tsObj-????_lens.fit',platelist
      IF platelist[0] EQ '' THEN message,'Error'
  endif
  nfiles=n_elements(platelist)
  
  ngood=0.
  nplates=0

  ntotal=0L
  ptrlist=ptrarr(nfiles)
  numlist=lonarr(nfiles)

  IF keyword_set(keepnomatch) THEN message,'Will keep objects that did not match',/inf

  struct_type = 'struct'+ntostr(long(systime(1)))
  for i=0,nfiles-1 do begin

      print
      print,'Reading lens plate '+ntostr(i+1)+'/'+ntostr(nfiles)+': '+platelist[i]
      
      pl=mrdfits(platelist(i),1,/silent)
      extract_from_plate,pl,pln,keepnomatch=keepnomatch
      
      IF n_elements(tstr) EQ 0 THEN BEGIN
          tstr = create_struct(name=struct_type, pln[0])
      ENDIF 

      if (pln(0).run ne 0) then BEGIN
          nn=n_elements(pln)

          ;; need to use named structures for comb_ptrstruct, all same
          ;; don't know if all pl are same, so do struct_assign
          tmp = replicate(tstr, nn)
          struct_assign, pln, tmp

          ntotal = ntotal + nn
          numlist[i] = nn
          ptrlist[i] = ptr_new(tmp,/no_copy)

          pl=0
          pln = 0
       endif
       
       print,'Ngood = '+ntostr(ntotal)
       
   endfor
 
   comb_ptrstruct, tstr, ptrlist, numlist, ptot
   

  
   return
end
