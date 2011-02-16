PRO combine_scrantoncat, clr

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: combine_scrantoncat, clr'
      return
  ENDIF 

  ;; extract petro[clr] < max_mag galaxies from big catalogs, with the tags we want

  CASE clr OF
      0: mag_max = 21.5
      ELSE: mag_max = 21.0
  ENDCASE 

  IF mag_max EQ 21.0 THEN mstr = '21' ELSE mstr = '21.5'

  indir='/sdss5/data0/wtheta/'
  clrstr = !colors[clr]

  outfile = indir+'wtheta-752-756-gal-allmag-'+clrstr+mstr+'.fit'
;  outfile = indir+'wtheta-752-756-gal-allmag-'+clrstr+mstr+'.fit'
  print
  print,'Outputting file: ',outfile
  print

  runs = [752,756]
;  runs=752
  runstr = ntostr(runs)

  nrun=n_elements(runs)
  nfpercol = 7
  ncol = 6
;  ncol=1

  galtype = 3

  keepst=create_struct('lambda',0d,$
                       'eta',0d,$
                       'petrocounts', fltarr(5),$
                       'gr',0.0)

  nf=ncol*nfpercol*nrun
  print
  print,'Number of files: ',nf
  print
  ptrlist = ptrarr(nf)
  numlist = lonarr(nf)

  ii=0L
  FOR i=0L,nrun-1 DO BEGIN 

      FOR icol=1,ncol DO BEGIN 

          cstr = ntostr(icol)
          FOR ifile=0, nfpercol-1 DO BEGIN 

              ;; read in one of the files
              ifstr = ntostr(ifile)
              file=indir+'wtheta-seg-'+runstr[i]+'-'+cstr+'-'+ifstr+'.fits'
              print,'Reading: ',file
              tmp=mrdfits_deja_vu(file,/silent, deja_vu=deja_vu)

              petrocounts = tmp.petrocounts-tmp.reddening
              gr=(tmp.counts_model[1]-tmp.reddening[1]) - $
                (tmp.counts_model[2]-tmp.reddening[2])
              w=where((petrocounts[clr,*] LE mag_max) AND $
                      (tmp.objc_type EQ galtype), nw)
              
              IF nw NE 0 THEN BEGIN 
                  keep = replicate(keepst, nw)
                  eq2survey, tmp[w].ra, tmp[w].dec, lambda, eta
                  keep.lambda = lambda
                  keep.eta = eta
                  keep.petrocounts = petrocounts[*,w]
                  keep.gr = gr[w]

                  ptrlist[ii] = ptr_new(keep, /no_copy)
                  numlist[ii] = nw

              ENDIF 
              
              ii = ii+1
          ENDFOR 

      ENDFOR 

  ENDFOR 

  comb_ptrstruct, keepst, ptrlist, numlist, outstruct

  print
  print,'Outputting file: ',outfile
  print,'Number of objects: '+ntostr(long(total(numlist)))
  print
  mwrfits2, outstruct, outfile, /create, /destroy

END 
