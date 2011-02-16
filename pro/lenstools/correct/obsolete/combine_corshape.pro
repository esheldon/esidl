PRO combine_corshape, runs, reruns, clr, struct, hirata=hirata

  ;; combine together corshape outputs for many runs

  IF keyword_set(hirata) THEN rsmearstr = '_h' ELSE rsmearstr=''

  nrun = n_elements(runs)
  ntot = nrun*6

  ;;binsize = 0.01
  nbins = 20
  nberbin = 180000L
;  min_e1 = -0.15
;  max_e1 =  0.10

;  min_e2 = -0.15
;  max_e2 =  0.10

  ndef = 18000000L
  galpsfe1 = replicate(-9999., ndef)
  galpsfe2 = galpsfe1
  gale1 = galpsfe1
  gale2 = galpsfe1
  gale1c = galpsfe1
  gale2c = galpsfe2

  beg = 0L
  ntot = 0L

  sdss_shapecorr_dir = sdssidl_config('shapecorr_dir')
  FOR ir=0L, nrun-1 DO BEGIN 
      rstr = ntostr(runs[ir])
      rstr2 = run2string(runs[ir])
      rrstr = ntostr(reruns[ir])
      FOR camcol=1,6 DO BEGIN 

          index = ir*6 + (camcol-1)

          cstr = ntostr(camcol)
          file = sdss_shapecorr_dir + 'corr'+rstr+'/'+rrstr+'/objcs/'+cstr+$
            '/'+'corshape_'+rstr2+'_'+cstr+'_'+$
            !colors[clr]+rsmearstr+'_N1.fit'

          IF fexist(file) THEN BEGIN
              print,'OK   ',file 
              
              t=mrdfits(file,1,/silent)
              ww=where(abs(t.gale1_corrected) LT 2 AND abs(t.gale2_corrected) LT 2, nw)
              ntot = ntot+nw

              binner, t.galpsfe1[ww], t.gale1[ww], binsize, e1xhist, e1yhist, e1ysig, rev_ind1, nume1psf, $
                      nbins=nbins
              binner, t.galpsfe1[ww], t.gale1_corrected[ww], binsize, e1xhistc, e1yhistc, e1ysigc, $
                      nbins=nbins
              binner, t.galpsfe2[ww], t.gale2[ww], binsize, e2xhist, e2yhist, e2ysig, rev_ind2, nume2psf, $
                      nbins=nbins
              binner, t.galpsfe2[ww], t.gale2_corrected[ww], binsize, e2xhistc, e2yhistc, e2ysigc,$
                      nbins=nbins

              IF (beg+nw-1) LE ndef THEN BEGIN 
                  galpsfe1[beg:beg+nw-1] = t.galpsfe1[ww]
                  galpsfe2[beg:beg+nw-1] = t.galpsfe2[ww]
                  gale1[beg:beg+nw-1] = t.gale1[ww]
                  gale2[beg:beg+nw-1] = t.gale2[ww]
                  gale1c[beg:beg+nw-1] = t.gale1_corrected[ww]
                  gale2c[beg:beg+nw-1] = t.gale2_corrected[ww]
              ENDIF ELSE BEGIN 
                  add_arrval, t.galpsfe1[ww], galpsfe1
                  add_arrval, t.gale1[ww], gale1
                  add_arrval, t.gale1_corrected[ww], gale1c

                  add_arrval, t.galpsfe2[ww], galpsfe2
                  add_arrval, t.gale2[ww], gale2
                  add_arrval, t.gale2_corrected[ww], gale2c
              ENDELSE 
              beg = beg+nw


              IF n_elements(mean_e1psf) EQ 0 THEN BEGIN 
                  ;;nhist1 = n_elements(e1xhist)-1
                  ;;nhist2 = n_elements(e2xhist)-1
                  nhist1 = nbins-1
                  nhist2 = nbins-1

                  ;; NOTE WE DUMP THE LAST BIN, since its always
                  ;; rather empty
                  mean_e1psf = fltarr(nrun*6, nhist1)
                  mean_e1gal = mean_e1psf
                  err_e1gal = mean_e1psf
                  mean_e1galc = mean_e1psf
                  err_e1galc = mean_e1psf

                  hist_e1psf = lonarr(nhist1)

                  mean_e2psf = fltarr(nrun*6, nhist2)
                  mean_e2gal = mean_e2psf
                  err_e2gal = mean_e2psf
                  mean_e2galc = mean_e2psf
                  err_e2galc = mean_e2psf

                  hist_e2psf = lonarr(nhist2)

              ENDIF 

              FOR binnum=0L, nhist1-1 DO BEGIN 

                  IF rev_ind1[binnum] NE rev_ind1[binnum+1] THEN BEGIN 
                      mean_e1psf[index, binnum] = e1xhist[binnum]
                      mean_e1gal[index, binnum] = e1yhist[binnum]
                      err_e1gal[index, binnum] = e1ysig[binnum]
                      mean_e1galc[index, binnum] = e1yhistc[binnum]
                      err_e1galc[index, binnum] = e1ysigc[binnum]

                      hist_e1psf[binnum] = hist_e1psf[binnum] + nume1psf[binnum]
                  ENDIF 
              ENDFOR 

              FOR binnum=0L, nhist2-1 DO BEGIN 

                  IF rev_ind2[binnum] NE rev_ind2[binnum+1] THEN BEGIN 
                      mean_e2psf[index, binnum] = e2xhist[binnum]
                      mean_e2gal[index, binnum] = e2yhist[binnum]
                      err_e2gal[index, binnum] = e2ysig[binnum]
                      mean_e2galc[index, binnum] = e2yhistc[binnum]
                      err_e2galc[index, binnum] = e2ysigc[binnum]

                      hist_e2psf[binnum] = hist_e2psf[binnum] + nume2psf[binnum]
                  ENDIF 
              ENDFOR 

          ENDIF ELSE print,'BAD   ',file
          

      ENDFOR 
  ENDFOR 

  print,'Ntot: ',ntot

  w=where(galpsfe1 NE -9999)
  galpsfe1 = galpsfe1[w]
  galpsfe2 = galpsfe2[w]
  gale1 = gale1[w]
  gale2 = gale2[w]
  gale1c = gale1c[w]
  gale2c = gale2c[w]

  struct = create_struct('galpsfe1', galpsfe1, $
                         'galpsfe2', galpsfe2, $
                         'gale1', gale1, $
                         'gale2', gale2, $
                         'gale1_corrected', gale1c, $
                         'gale2_corrected', gale2c, $
                         'mean_e1psf', mean_e1psf, $
                         'mean_e2psf', mean_e2psf, $
                         'mean_e1gal', mean_e1gal, $
                         'mean_e2gal', mean_e2gal, $
                         'mean_e1galc', mean_e1galc, $
                         'mean_e2galc', mean_e2galc, $
                         'err_e1gal', err_e1gal, $
                         'err_e2gal', err_e2gal, $
                         'err_e1galc', err_e1galc, $
                         'err_e2galc', err_e2galc, $
                         'hist_e1psf', hist_e1psf, $
                         'hist_e2psf', hist_e2psf)

  galpsfe1 = 0
  galpsfe2 = 0
  gale1 = 0
  gale2 = 0
  gale1c = 0
  gale2c = 0

  ;; construct the overall mean, variance
  ;; NOTE WE DUMP THE LAST BIN, since its always
  ;; rather empty

  tmean_e1psf = fltarr(nhist1)
  tmean_e1gal = tmean_e1psf
  terr_e1gal = tmean_e1psf
  tmean_e1galc = tmean_e1psf
  terr_e1galc = tmean_e1psf
  
  tmean_e2psf = fltarr(nhist2)
  tmean_e2gal = tmean_e2psf
  terr_e2gal = tmean_e2psf
  tmean_e2galc = tmean_e2psf
  terr_e2galc = tmean_e2psf

  FOR i=0L, nhist1-1 DO BEGIN 

      w=where(mean_e1psf[*,i] NE 0)
      tmean_e1psf[i] = mean(mean_e1psf[w,i])

      mom = moment(mean_e1gal[w, i])
      tmean_e1gal[i] = mom[0]
      terr_e1gal[i] = sqrt(mom[1]/ntot)

      mom = moment(mean_e1galc[w, i])
      tmean_e1galc[i] = mom[0]
      terr_e1galc[i] = sqrt(mom[1]/ntot)


  ENDFOR 

  FOR i=0L, nhist2-1 DO BEGIN 

      w=where(mean_e1psf[*,i] NE 0)
      tmean_e2psf[i] = mean(mean_e2psf[w,i])

      mom = moment(mean_e2gal[w, i])
      tmean_e2gal[i] = mom[0]
      terr_e2gal[i] = sqrt(mom[1]/ntot)

      mom = moment(mean_e2galc[w, i])
      tmean_e2galc[i] = mom[0]
      terr_e2galc[i] = sqrt(mom[1]/ntot)

  ENDFOR 

  struct = create_struct(struct, $
                         'tmean_e1psf',  tmean_e1psf, $
                         'tmean_e2psf',  tmean_e2psf, $
                         'tmean_e1gal',  tmean_e1gal, $
                         'tmean_e2gal',  tmean_e2gal, $
                         'tmean_e1galc', tmean_e1galc, $
                         'tmean_e2galc', tmean_e2galc, $
                         'terr_e1gal',   terr_e1gal, $
                         'terr_e2gal',   terr_e2gal, $
                         'terr_e1galc',  terr_e1galc, $
                         'terr_e2galc',  terr_e2galc)



return
END 
