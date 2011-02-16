PRO run_spectra_plot, dojack=dojack, type=type, lumclr=lumclr, start=start, comoving=comoving, logbin=logbin

  hirata=1
  recorr=1
  dops=1
  clr=[1,2,3]

  stripes = [9,10,11,12,13,14,15,$
             27,28,29,30,31,32,33,34,35,36,37,$
             76,86]

  combined_stripes = 1
  cstripe_strings = stripearr2string(stripes)

;  GOTO, jump

  spectra_plot, stripes, clr, dops=dops, $
    dojack=dojack, $
    combined_stripes=combined_stripes, $
    cstripe_strings=cstripe_strings, $
    recorr=recorr, logbin=logbin, $
    hirata=hirata

  start = 2
  spectra_plot, stripes, clr, dops=dops, $
    dojack=dojack, $
    combined_stripes=combined_stripes, $
    cstripe_strings=cstripe_strings, $
    recorr=recorr, logbin=logbin, $
    hirata=hirata, /rlrg_sources, start=start


  type = 'matchLRGNew'
  spectra_plot, stripes, clr, dops=dops, type=type, $
    dojack=dojack, $
    combined_stripes=combined_stripes, $
    cstripe_strings=cstripe_strings, $
    recorr=recorr, logbin=logbin, $
    hirata=hirata

  type = 'matchAll'
  start = 2
  spectra_plot, stripes, clr, dops=dops, type=type, $
    dojack=dojack, $
    combined_stripes=combined_stripes, $
    cstripe_strings=cstripe_strings, $
    recorr=recorr, logbin=logbin, $
    hirata=hirata, /rlrg_sources, start=start

;jump:

  FOR lumclr=2,2 DO BEGIN 

      FOR iclr=1,3 DO BEGIN 

          type = 'lum'+ntostr(iclr)+'threebinnum_matchLRGNew'
          spectra_plot, stripes, clr, type=type, $
            dops=dops, $
            dojack=dojack, $
            combined_stripes=combined_stripes, $
            cstripe_strings=cstripe_strings, $
            recorr=recorr, logbin=logbin, $
            lumclr=lumclr, hirata=hirata, $
            comoving=comoving

          start = 2
          type = 'lum'+ntostr(iclr)+'threebinnum_matchAll'
          spectra_plot, stripes, clr, type=type, start=start, $
            dops=dops, $
            dojack=dojack, $
            combined_stripes=combined_stripes, $
            cstripe_strings=cstripe_strings, $
            recorr=recorr, logbin=logbin, $
            lumclr=lumclr, hirata=hirata, $
            comoving=comoving, /rlrg_sources

      ENDFOR 
  ENDFOR 

  

return

  FOR lumclr=0,4 DO BEGIN 

      FOR iclr=1,4 DO BEGIN 

          type = 'lum'+ntostr(iclr)+'fourbinnum'
          spectra_plot, stripes, clr, type=type, $
            dops=dops, start=start, nf=nf, $
            dojack=dojack, $
            combined_stripes=combined_stripes, $
            cstripe_strings=cstripe_strings, $
            recorr=recorr, logbin=logbin, $
            lumclr=lumclr, hirata=hirata, $
            comoving=comoving, lrg_sources=lrg_sources
      ENDFOR 
  ENDFOR 


return


  stripes1 = [9,10,11,12,13,14,15]
  stripes2 = [27,28,29,30,31,32,33,34,35,36,37]
  stripes3 = 76
  stripes4 = 86
  cst1 = stripearr2string(stripes1)
  cst2 = stripearr2string(stripes2)
  cst3 = '76'
  cst4 = '86'

  combined_stripes=1

  cstripe_strings = [cst1, cst2]
  stripes = [stripes1, stripes2]

  spectra_plot, stripes1, clr, dops=dops, type=type, start=start, nf=nf, $
    dojack=dojack, indir=indir, $
    combined_stripes=combined_stripes, $
    cstripe_strings=cst1, $
    recorr=recorr, logbin=logbin, $
    lumclr=lumclr, hirata=hirata, comoving=comoving, lrg_sources=lrg_sources

  spectra_plot, stripes2, clr, dops=dops, type=type, start=start, nf=nf, $
    dojack=dojack, $
    combined_stripes=combined_stripes, $
    cstripe_strings=cst2, $
    recorr=recorr, logbin=logbin, $
    lumclr=lumclr, hirata=hirata, comoving=comoving, lrg_sources=lrg_sources

  spectra_plot, stripes, clr, dops=dops, type=type, start=start, nf=nf, $
    dojack=dojack, $
    combined_stripes=combined_stripes, $
    cstripe_strings=cstripe_strings, $
    recorr=recorr, logbin=logbin, $
    lumclr=lumclr, hirata=hirata, comoving=comoving, lrg_sources=lrg_sources

return

  FOR lumclr=0,4 DO BEGIN 
      spectra_plot, stripes, clr, dops=dops, type='lum1threebin', $
        start=start, nf=nf, $
        dojack=dojack, $
        combined_stripes=combined_stripes, $
        cstripe_strings=cstripe_strings, $
        recorr=recorr, logbin=logbin, $
        lumclr=lumclr, hirata=hirata, comoving=comoving, $
        lrg_sources=lrg_sources

      spectra_plot, stripes, clr, dops=dops, type='lum2threebin', $
        start=start, nf=nf, $
        dojack=dojack, $
        combined_stripes=combined_stripes, $
        cstripe_strings=cstripe_strings, $
        recorr=recorr, logbin=logbin, $
        lumclr=lumclr, hirata=hirata, comoving=comoving, $
        lrg_sources=lrg_sources
      
      spectra_plot, stripes, clr, dops=dops, type='lum3threebin', $
        start=start, nf=nf, $
        dojack=dojack, $
        combined_stripes=combined_stripes, $
        cstripe_strings=cstripe_strings, $
        recorr=recorr, logbin=logbin, $
        lumclr=lumclr, hirata=hirata, comoving=comoving, $
        lrg_sources=lrg_sources

  ENDFOR 

return

  FOR lumclr=0,4 DO BEGIN 
      spectra_plot, stripes, clr, dops=dops, type='lum1threebinnum', $
        start=start, nf=nf, $
        dojack=dojack, $
        combined_stripes=combined_stripes, $
        cstripe_strings=cstripe_strings, $
        recorr=recorr, logbin=logbin, $
        lumclr=lumclr, hirata=hirata, comoving=comoving, $
        lrg_sources=lrg_sources

      spectra_plot, stripes, clr, dops=dops, type='lum2threebinnum', $
        start=start, nf=nf, $
        dojack=dojack, $
        combined_stripes=combined_stripes, $
        cstripe_strings=cstripe_strings, $
        recorr=recorr, logbin=logbin, $
        lumclr=lumclr, hirata=hirata, comoving=comoving, $
        lrg_sources=lrg_sources

      spectra_plot, stripes, clr, dops=dops, type='lum3threebinnum', $
        start=start, nf=nf, $
        dojack=dojack, $
        combined_stripes=combined_stripes, $
        cstripe_strings=cstripe_strings, $
        recorr=recorr, logbin=logbin, $
        lumclr=lumclr, hirata=hirata, comoving=comoving, $
        lrg_sources=lrg_sources

  ENDFOR 



  FOR lumclr=0,4 DO BEGIN 

      spectra_plot, stripes, clr, dops=dops, type='redlum1twobin', $
        start=start, nf=nf, $
        dojack=dojack, $
        combined_stripes=combined_stripes, $
        cstripe_strings=cstripe_strings, $
        recorr=recorr, logbin=logbin, $
        lumclr=lumclr, hirata=hirata, comoving=comoving, $
        lrg_sources=lrg_sources

      spectra_plot, stripes, clr, dops=dops, type='redlum2twobin', $
        start=start, nf=nf, $
        dojack=dojack, $
        combined_stripes=combined_stripes, $
        cstripe_strings=cstripe_strings, $
        recorr=recorr, logbin=logbin, $
        lumclr=lumclr, hirata=hirata, comoving=comoving, $
        lrg_sources=lrg_sources
      
  ENDFOR 

  FOR lumclr=0,4 DO BEGIN 

      spectra_plot, stripes, clr, dops=dops, type='redlum1threebin', $
        start=start, nf=nf, $
        dojack=dojack, $
        combined_stripes=combined_stripes, $
        cstripe_strings=cstripe_strings, $
        recorr=recorr, logbin=logbin, $
        lumclr=lumclr, hirata=hirata, comoving=comoving, $
        lrg_sources=lrg_sources

      spectra_plot, stripes, clr, dops=dops, type='redlum2threebin', $
        start=start, nf=nf, $
        dojack=dojack, $
        combined_stripes=combined_stripes, $
        cstripe_strings=cstripe_strings, $
        recorr=recorr, logbin=logbin, $
        lumclr=lumclr, hirata=hirata, comoving=comoving, $
        lrg_sources=lrg_sources

      spectra_plot, stripes, clr, dops=dops, type='redlum3threebin', $
        start=start, nf=nf, $
        dojack=dojack, $
        combined_stripes=combined_stripes, $
        cstripe_strings=cstripe_strings, $
        recorr=recorr, logbin=logbin, $
        lumclr=lumclr, hirata=hirata, comoving=comoving, $
        lrg_sources=lrg_sources

  ENDFOR 

  return
END 
