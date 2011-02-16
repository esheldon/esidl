PRO spectra_sub_meane, stripe, type, $
                       ext=ext, $
                       basedir=basedir, indir=indir, $
                       nooverwrite=nooverwrite, $
                       recorr=recorr, ri=ri, gri=gri, hirata=hirata,$
                       lrg_sources=lrg_sources
  
  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: spectra_sub_meane, stripe, type, ext=ext, basedir=basedir, indir=indir, nooverwrite=nooverwrite, lrg_sources=lrg_sources'
      print
      print,'Allowed types: main, mainseeing, lrg, lumtwobin, lumfourbin'
      return
  ENDIF 

  ;; BASEDIR is base, under which there are stripe directories, etc.
  ;; INDIR is the actual directory, ignoring basedir

  eclasscut = '-0.06'
  ;;eclasscut = '-0.02'
  
  gmrcut = '0.7'
  gmrmin = '0.1'
  gmrmax = '1.1'

  minz   = '0.02'
  maxz   = '0.3'


  IF keyword_set(nooverwrite) THEN overwrite=0 ELSE overwrite=1
  IF keyword_set(recorr) THEN recorrstr='_recorr' ELSE recorrstr=''

  IF keyword_set(ri) AND keyword_set(gri) THEN message,'Do not send both /ri and /gri'
  IF NOT keyword_set(ri) AND NOT keyword_set(gri) THEN gri=1
  IF keyword_set(hirata) THEN hirstr = '_h' ELSE hirstr = ''
  IF keyword_set(lrg_sources) THEN lrgstr = '_lrg' ELSE lrgstr=''

  msun = [6.39,5.07,4.62,4.52,4.48]
  mstar = [-18.34,-20.04,-20.83,-21.26,-21.55]
  lstar = 10.0^((mstar-msun)/(-2.5))

  IF n_elements(ext) EQ 0 THEN ext = 'N1.fit'

  stripestr = 'stripe'+stripearr2string(stripe)+'_'

  IF n_elements(basedir) EQ 0 THEN basedir='/net/cheops2/home/esheldon/lensout/'
  
  IF n_elements(indir) EQ 0 THEN BEGIN
      indir = basedir+'stripe'+stripearr2string(stripe)+'/'
  ENDIF 

  fri = indir + 'zgal_gal_'+stripestr+'ri'+recorrstr+hirstr+'_lensum_'+ext
  fgri = indir + 'zgal_gal_'+stripestr+'gri'+recorrstr+hirstr+'_lensum_'+ext

  rfri = indir + 'zrand_'+stripestr+'ri'+recorrstr+hirstr+'_lensum_'+ext
  rfgri = indir + 'zrand_'+stripestr+'gri'+recorrstr+hirstr+'_lensum_'+ext

  IF keyword_set(ri) THEN BEGIN 
      file = fri
      rfile = rfri
  ENDIF 
  IF keyword_set(gri) THEN BEGIN 
      file = fgri
      rfile = rfgri
  ENDIF 

  CASE type OF
      'zsub': BEGIN 

          wstr = '(lensum.z le '+maxz+') and (lensum.z ge '+minz+')'

          addstr = 'zsub'
          print
          print,'Doing Sub Sample: ',addstr

          zsub_sample, file, rfile, wstr, $
                       addstr=addstr, $
                       overwrite=overwrite
          return
      END 
      'main': BEGIN 

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Just main galaxy sample (obsolete since we cut
          ;; out LRG's before hand)
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          wstr = '(lensum.z le '+maxz+') and (lensum.z ge '+minz+')'
          wstr = wstr + ' and ( (lensum.primtarget and 2L^6) ne 0)'
          
          addstr = 'main'
          print
          print,'Doing Sub Sample: ',addstr

          zsub_sample, file, rfile, wstr, $
                       addstr=addstr, $
                       overwrite=overwrite
          return
      END 
      'mainseeing': BEGIN 
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Just main galaxy sample with seeing cut
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          wstr = '(lensum.z LE '+maxz+') AND (lensum.z GE '+minz+')'
          wstr = wstr + ' AND ( (lensum.primtarget AND 2L^6) NE 0)'
          wstr = wstr + ' AND ( lensum.seeing[2] LT 1.75 )'
          wstr = wstr + ' AND ( lensum.seeing[2] GT 0.0  )'
          
          addstr = 'mainseeing'
          print
          print,'Doing Sub Sample: ',addstr
          zsub_sample, file, rfile, wstr, $
                       addstr=addstr, $
                       overwrite=overwrite
          return
      END 
      'vlim': BEGIN 
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Just main galaxy sample but volume limited
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          minz = '0.02'
          maxz = '0.1'
          maxmag = '17.6'

          wstr = '(lensum.z LE '+maxz+') AND (lensum.z GE '+minz+')'
          wstr = wstr + ' AND ( (lensum.primtarget AND 2L^6) NE 0)'
          wstr = wstr + ' AND ( lensum.abscounts[2] LT '+maxmag+' )'
          
          addstr = 'vlim'
          print
          print,'Doing Sub Sample: ',addstr
          zsub_sample, file, rfile, wstr, $
                       addstr=addstr, $
                       overwrite=overwrite
      END 
      'vlim1': BEGIN 
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Now volume-absmag limited samples
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          minz = '0.027'
          maxz = '0.051'         
          minmag = '-20.0'
          maxmag = '-18.5'

          wstr = '(lensum.z LE '+maxz+') AND (lensum.z GE '+minz+')'
          wstr = wstr + ' AND ( (lensum.primtarget AND 2L^6) NE 0)'
          wstr = wstr + ' AND ( lensum.absmag[2] LE '+maxmag+' )'
          wstr = wstr + ' AND ( lensum.absmag[2] GT '+minmag+' )'

          addstr = 'vlim1'
          print
          print,'Doing Sub Sample: ',addstr
          zsub_sample, file, rfile, wstr, $
                       addstr=addstr, $
                       overwrite=overwrite
      END 
      'vlim2': BEGIN 
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Now volume-absmag limited samples
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          minz = '0.052'
          maxz = '0.097'
          minmag = '-21.5'
          maxmag = '-20.0'
          

          wstr = '(lensum.z LE '+maxz+') AND (lensum.z GE '+minz+')'
          wstr = wstr + ' AND ( (lensum.primtarget AND 2L^6) NE 0)'
          wstr = wstr + ' AND ( lensum.absmag[2] LE '+maxmag+' )'
          wstr = wstr + ' AND ( lensum.absmag[2] GT '+minmag+' )'

          addstr = 'vlim2'
          print
          print,'Doing Sub Sample: ',addstr
          zsub_sample, file, rfile, wstr, $
                       addstr=addstr, $
                       overwrite=overwrite
      END 
      'vlim3': BEGIN 
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Now volume-absmag limited samples
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          minz = '0.1'
          maxz = '0.174'          
          minmag = '-23.0'
          maxmag = '-21.5'

          wstr = '(lensum.z LE '+maxz+') AND (lensum.z GE '+minz+')'
          wstr = wstr + ' AND ( (lensum.primtarget AND 2L^6) NE 0)'
          wstr = wstr + ' AND ( lensum.absmag[2] LE '+maxmag+' )'
          wstr = wstr + ' AND ( lensum.absmag[2] GT '+minmag+' )'

          addstr = 'vlim3'
          print
          print,'Doing Sub Sample: ',addstr
          zsub_sample, file, rfile, wstr, $
                       addstr=addstr, $
                       overwrite=overwrite
      END 
      'lrg': BEGIN 
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; LRG
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          maxz = '0.6'
          minz = '0.15'

          wstr = '(lensum.z le '+maxz+') and (lensum.z ge '+minz+')'
          wstr = wstr + ' and ( (lensum.primtarget and 2L^5) ne 0)'
  
          addstr = 'lrg'
          print
          print,'Doing Sub Sample: ',addstr
          zsub_sample, file, rfile, wstr, $
                       addstr=addstr, $
                       overwrite=overwrite

      END 
      'lumtwobin': BEGIN 

          print
          print,'Reading: ',file
          len = mrdfits(file, 1)
          print,'Reading: ',rfile
          rlen = mrdfits(rfile, 1)

          FOR iclr = 0,4 DO BEGIN 
              outdir = indir + 'sublum/'+!colors[iclr]+'/'

              ;;!!!! BINBYNUM HAS BUILT-IN Z CUTS!!!! 0.6,0.02
              binbynum_2bin, len, iclr, w1, w2

              addstr = 'lum1twobin'
              print,'Sub sample: ',addstr
              zsub_sample, file, rfile, $
                           addstr=addstr,uselens=len,userand=rlen, $
                           overwrite=overwrite, indices=w1, $
                           outdir=outdir
                  
              addstr = 'lum2twobin'
              print,'Sub sample: ',addstr
              zsub_sample, file, rfile, $
                           addstr=addstr,uselens=len,userand=rlen, $
                           overwrite=overwrite, indices=w2, $
                           outdir=outdir


          ENDFOR 

          delvarx, len, rlen

      END 
      'lumthreebin': BEGIN 

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; only bin by r-mag for now
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          print
          print,'Reading: ',file
          len = mrdfits(file, 1)
          print,'Reading: ',rfile
          rlen = mrdfits(rfile, 1)

          FOR clr=0,4 DO BEGIN 
              clrstr = ntostr(clr)
              print
              print,'Doing '+!colors[clr]+'-band sublum'
              print

              outdir = indir + 'sublum/'+!colors[clr]+'/'

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; first bin
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              ;; zehavi: -20.0 < M_r    < -18.5
              ;;          0.83 < M_r-M* < 2.33
              ;;    0.46558609 > L/L*   > 0.11694994
              ;; absmag_min = '-20.0'
              ;; absmag_max = '-16.5'
              mdiff_min = '0.83'
              mdiff_max = '2.33'
              
              wstr = $
                '(lensum.z LE '+maxz+') AND (lensum.z GE '+minz+')'
              wstr = wstr + $
                ' AND ( (lensum.primtarget and 2L^6) NE 0)'
              wstr = wstr + $
                ' AND ( lensum.absmag['+clrstr+'] ne -9999. )'
              wstr = wstr + $
                ' AND ( lensum.absmag['+clrstr+'] - !mstar['+clrstr+'] GT '+mdiff_min+' )'
              wstr = wstr + $
                ' AND ( lensum.absmag['+clrstr+'] - !mstar['+clrstr+'] LE '+mdiff_max+' )'
              
              addstr = 'lum1threebin'
              print,'Sub sample: ',addstr
              zsub_sample, file, rfile, wstr, $
                           addstr=addstr,uselens=len,userand=rlen, $
                           overwrite=overwrite, outdir=outdir
              
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; Second bin
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              ;; zehavi: -21.5 < M_r      < -20.0
              ;;         -0.67 < M_r - M* < 0.83
              ;;     1.8535316 > L/L*     > 0.46558609
              ;; absmag_min = '-21.5'
              ;; absmag_max = '-20.0'
              mdiff_min = '-0.67'
              mdiff_max = '0.83'
              
              wstr = $
                '(lensum.z LE '+maxz+') AND (lensum.z GE '+minz+')'
              wstr = wstr + $
                ' AND ( (lensum.primtarget and 2L^6) NE 0)'
              wstr = wstr + $
                ' AND ( lensum.absmag['+clrstr+'] ne -9999. )'
              wstr = wstr + $
                ' AND ( lensum.absmag['+clrstr+'] - !mstar['+clrstr+'] GT '+mdiff_min+' )'
              wstr = wstr + $
                ' AND ( lensum.absmag['+clrstr+'] - !mstar['+clrstr+'] LE '+mdiff_max+' )'
              
              addstr = 'lum2threebin'
              print,'Sub sample: ',addstr
              zsub_sample, file, rfile, wstr, $
                           addstr=addstr,uselens=len,userand=rlen, $
                           overwrite=overwrite, outdir=outdir
              
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; Third bin
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              ;; zehavi: -23.0 < M_r      < -21.5
              ;;         -2.17 < M_r - M* < -0.67
              ;;     7.3790423 > L/L*     > 1.8535316
              ;; absmag_min = '-23.0'
              ;; absmag_max = '-21.5'
              mdiff_min = '-2.17'
              mdiff_max = '-0.67'
              
              wstr = $
                '(lensum.z LE '+maxz+') AND (lensum.z GE '+minz+')'
              wstr = wstr + $
                ' AND ( (lensum.primtarget and 2L^6) NE 0)'
              wstr = wstr + $
                ' AND ( lensum.absmag['+clrstr+'] ne -9999. )'
              wstr = wstr + $
                ' AND ( lensum.absmag['+clrstr+'] - !mstar['+clrstr+'] GT '+mdiff_min+' )'
              wstr = wstr + $
                ' AND ( lensum.absmag['+clrstr+'] - !mstar['+clrstr+'] LE '+mdiff_max+' )'
              
              addstr = 'lum3threebin'
              print,'Sub sample: ',addstr
              zsub_sample, file, rfile, wstr, $
                           addstr=addstr,uselens=len,userand=rlen, $
                           overwrite=overwrite, outdir=outdir
          ENDFOR 
          delvarx, len, rlen

      END 
      'lumthreebinnum': BEGIN 

          read_cuts, range_struct
          
          FOR clr=0, 4 DO BEGIN 

              print
              print,'** Binning in '+!colors[clr]+'-mag'
              print
              cstr = ntostr(clr)
              outdir = indir + 'sublum/'+!colors[clr]+'/'

              lowcut = ntostr(range_struct.main_minmag_threebin[clr,0])
              highcut = ntostr(range_struct.main_maxmag_threebin[clr,0])
              wstr = $
                'lensum.absmag['+cstr+'] le '+highcut+' AND '+$
                'lensum.absmag['+cstr+'] gt '+lowcut+ ' AND '+$
                'lensum.z gt '+minz+' AND '+$
                'lensum.z lt '+maxz
              
              addstr = 'lum1threebinnum'
              print,'Sub sample: ',addstr
              print,lowcut+' < M['+!colors[clr]+'] < '+highcut
              zsub_sample, file, rfile, wstr, $
                           addstr=addstr,uselens=len,userand=rlen, $
                           overwrite=overwrite, indices=w1, $
                           outdir=outdir
              
              lowcut = ntostr(range_struct.main_minmag_threebin[clr,1])
              highcut = ntostr(range_struct.main_maxmag_threebin[clr,1])
              wstr = $
                'lensum.absmag['+cstr+'] le '+highcut+' AND '+$
                'lensum.absmag['+cstr+'] gt '+lowcut+ ' AND '+$
                'lensum.z gt '+minz+' AND '+$
                'lensum.z lt '+maxz
              
              addstr = 'lum2threebinnum'
              print,'Sub sample: ',addstr
              print,lowcut+' < M['+!colors[clr]+'] < '+highcut
              zsub_sample, file, rfile, wstr, $
                           addstr=addstr,uselens=len,userand=rlen, $
                           overwrite=overwrite, indices=w2, $
                           outdir=outdir

              lowcut = ntostr(range_struct.main_minmag_threebin[clr,2])
              highcut = ntostr(range_struct.main_maxmag_threebin[clr,2])
              wstr = $
                'lensum.absmag['+cstr+'] le '+highcut+' AND '+$
                'lensum.absmag['+cstr+'] gt '+lowcut+ ' AND '+$
                'lensum.z gt '+minz+' AND '+$
                'lensum.z lt '+maxz
              
              addstr = 'lum3threebinnum'
              print,'Sub sample: ',addstr
              print,lowcut+' < M['+!colors[clr]+'] < '+highcut
              zsub_sample, file, rfile, wstr, $
                           addstr=addstr,uselens=len,userand=rlen, $
                           overwrite=overwrite, indices=w2, $
                           outdir=outdir


          ENDFOR 
          delvarx, len, rlen

      END 
      'lumfourbinnum': BEGIN 

          print
          print,'Reading: ',file
          len = mrdfits(file, 1)
          print,'Reading: ',rfile
          rlen = mrdfits(rfile, 1)

          FOR iclr = 0,4 DO BEGIN 
              outdir = indir + 'sublum/'+!colors[iclr]+'/'

              ;;!!!! BINBYNUM HAS BUILT-IN Z CUTS!!!! 0.6,0.02
              binbynum_4bin, len, iclr, w1, w2, w3, w4

              addstr = 'lum1fourbinnum'
              print,'Sub sample: ',addstr
              zsub_sample, file, rfile, $
                           addstr=addstr,uselens=len,userand=rlen, $
                           overwrite=overwrite, indices=w1, $
                           outdir=outdir
              
              addstr = 'lum2fourbinnum'
              print,'Sub sample: ',addstr
              zsub_sample, file, rfile, $
                           addstr=addstr,uselens=len,userand=rlen, $
                           overwrite=overwrite, indices=w2, $
                           outdir=outdir
              
              addstr = 'lum3fourbinnum'
              print,'Sub sample: ',addstr
              zsub_sample, file, rfile, $
                           addstr=addstr,uselens=len,userand=rlen, $
                           overwrite=overwrite, indices=w3, $
                           outdir=outdir
              
              addstr = 'lum4fourbinnum'
              print,'Sub sample: ',addstr
              zsub_sample, file, rfile, $
                           addstr=addstr,uselens=len,userand=rlen, $
                           overwrite=overwrite, indices=w4, $
                           outdir=outdir

          ENDFOR 

          delvarx, len, rlen

      END 
      'mgtr': BEGIN 
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; three samples: v > vmin
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          wstr = $
            '(lensum.z LE '+maxz+') AND (lensum.z GE '+minz+')'
          wstr = wstr + $
            ' AND ( (lensum.primtarget and 2L^6) NE 0)'
          wstr = wstr + $
            ' AND ( lensum.absmag[2] lt -18 )'
          wstr = wstr + $
            ' and ( lensum.absmag[2] gt -24 )'

          addstr = 'mgtr18'
          print
          print,'Doing Sub Sample: ',addstr
          zsub_sample, file, rfile, wstr, $
                       addstr=addstr, $
                       overwrite=overwrite

          wstr = $
            '(lensum.z LE '+maxz+') AND (lensum.z GE '+minz+')'
          wstr = wstr + $
            ' AND ( (lensum.primtarget and 2L^6) NE 0)'
          wstr = wstr + $
            ' AND ( lensum.absmag[2] lt -19 )'
          wstr = wstr + $
            ' and ( lensum.absmag[2] gt -24 )'

          addstr = 'mgtr19'
          print
          print,'Doing Sub Sample: ',addstr
          zsub_sample, file, rfile, wstr, $
                       addstr=addstr, $
                       overwrite=overwrite

          wstr = $
            '(lensum.z LE '+maxz+') AND (lensum.z GE '+minz+')'
          wstr = wstr + $
            ' AND ( (lensum.primtarget and 2L^6) NE 0)'
          wstr = wstr + $
            ' AND ( lensum.absmag[2] lt -19.8 )'
          wstr = wstr + $
            ' and ( lensum.absmag[2] gt -24 )'

          addstr = 'mgtr198'
          print
          print,'Doing Sub Sample: ',addstr
          zsub_sample, file, rfile, wstr, $
                       addstr=addstr, $
                       overwrite=overwrite


      END 
      'gmr': BEGIN 
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Just main galaxy sample with a cut in g-r
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          print
          print,'Reading: ',file
          len = mrdfits(file, 1)
          print,'Reading: ',rfile
          rlen = mrdfits(rfile, 1)

          wstr = $
            '(lensum.z LE '+maxz+') AND (lensum.z GE '+minz+')'
          wstr = wstr + $
            ' AND ( (lensum.primtarget and 2L^6) NE 0)'
          wstr = wstr + $
            ' AND ( (lensum.abscounts[1]-lensum.abscounts[2]) LT '+gmrmax+' )'
          wstr = wstr + $
            ' AND ( (lensum.abscounts[1]-lensum.abscounts[2]) GE '+gmrmin+' )'
          wstr = wstr + $
            ' AND ( (lensum.abscounts[1]-lensum.abscounts[2]) LT '+gmrcut+' )'

          addstr = 'gmr1'
          print
          print,'Doing Sub Sample: ',addstr
          zsub_sample, file, rfile, wstr, $
                       addstr=addstr,uselens=len,userand=rlen, $
                       overwrite=overwrite

          wstr = $
            '(lensum.z LE '+maxz+') AND (lensum.z GE '+minz+')'
          wstr = wstr + $
            ' AND ( (lensum.primtarget and 2L^6) NE 0)'
          wstr = wstr + $
            ' AND ( (lensum.abscounts[1]-lensum.abscounts[2]) LT '+gmrmax+' )'
          wstr = wstr + $
            ' AND ( (lensum.abscounts[1]-lensum.abscounts[2]) GE '+gmrmin+' )'
          wstr = wstr + $
            ' AND ( (lensum.abscounts[1]-lensum.abscounts[2]) GE '+gmrcut+' )'
          
          addstr = 'gmr2'
          print
          print,'Doing Sub Sample: ',addstr
          zsub_sample, file, rfile, wstr, $
                       addstr=addstr,uselens=len,userand=rlen, $
                       overwrite=overwrite
          delvarx, len, rlen

      END 
      'redlumtwobin': BEGIN 

          read_cuts, range_struct
          
          FOR clr=0, 4 DO BEGIN 

              cstr = ntostr(clr)
              outdir = indir + 'sublum/'+!colors[clr]+'/'

              lowcut = ntostr(range_struct.red_minmag_twobin[clr,0])
              highcut = ntostr(range_struct.red_maxmag_twobin[clr,0])
              wstr = $
                '(lensum.abscounts[1]-lensum.abscounts[2]) LT '+gmrmax+' AND '+$
                '(lensum.abscounts[1]-lensum.abscounts[2]) GE '+gmrmin+' AND '+$
                '(lensum.abscounts[1]-lensum.abscounts[2]) GE '+gmrcut+' AND '+$
                'lensum.absmag['+cstr+'] le '+highcut+' AND '+$
                'lensum.absmag['+cstr+'] gt '+lowcut+ ' AND '+$
                'lensum.z gt '+minz+' AND '+$
                'lensum.z lt '+maxz
              
              addstr = 'redlum1twobin'
              print,'Sub sample: ',addstr
              print,lowcut+' < M['+!colors[clr]+'] < '+highcut
              zsub_sample, file, rfile, wstr, $
                           addstr=addstr,uselens=len,userand=rlen, $
                           overwrite=overwrite, indices=w1, $
                           outdir=outdir
              
              lowcut = ntostr(range_struct.red_minmag_twobin[clr,1])
              highcut = ntostr(range_struct.red_maxmag_twobin[clr,1])
              wstr = $
                '(lensum.abscounts[1]-lensum.abscounts[2]) LT '+gmrmax+' AND '+$
                '(lensum.abscounts[1]-lensum.abscounts[2]) GE '+gmrmin+' AND '+$
                '(lensum.abscounts[1]-lensum.abscounts[2]) GE '+gmrcut+' AND '+$
                'lensum.absmag['+cstr+'] le '+highcut+' AND '+$
                'lensum.absmag['+cstr+'] gt '+lowcut+ ' AND '+$
                'lensum.z gt '+minz+' AND '+$
                'lensum.z lt '+maxz
              
              addstr = 'redlum2twobin'
              print,'Sub sample: ',addstr
              print,lowcut+' < M['+!colors[clr]+'] < '+highcut
              zsub_sample, file, rfile, wstr, $
                           addstr=addstr,uselens=len,userand=rlen, $
                           overwrite=overwrite, indices=w2, $
                           outdir=outdir

          ENDFOR 
          delvarx, len, rlen

      END 
      'redlumthreebin': BEGIN 

          read_cuts, range_struct
          
          FOR clr=0, 4 DO BEGIN 

              cstr = ntostr(clr)
              outdir = indir + 'sublum/'+!colors[clr]+'/'

              lowcut = ntostr(range_struct.red_minmag_threebin[clr,0])
              highcut = ntostr(range_struct.red_maxmag_threebin[clr,0])
              wstr = $
                '(lensum.abscounts[1]-lensum.abscounts[2]) LT '+gmrmax+' AND '+$
                '(lensum.abscounts[1]-lensum.abscounts[2]) GE '+gmrmin+' AND '+$
                '(lensum.abscounts[1]-lensum.abscounts[2]) GE '+gmrcut+' AND '+$
                'lensum.absmag['+cstr+'] le '+highcut+' AND '+$
                'lensum.absmag['+cstr+'] gt '+lowcut+ ' AND '+$
                'lensum.z gt '+minz+' AND '+$
                'lensum.z lt '+maxz
              
              addstr = 'redlum1threebin'
              print,'Sub sample: ',addstr
              print,lowcut+' < M['+!colors[clr]+'] < '+highcut
              zsub_sample, file, rfile, wstr, $
                           addstr=addstr,uselens=len,userand=rlen, $
                           overwrite=overwrite, indices=w1, $
                           outdir=outdir

              lowcut = ntostr(range_struct.red_minmag_threebin[clr,1])
              highcut = ntostr(range_struct.red_maxmag_threebin[clr,1])
              wstr = $
                '(lensum.abscounts[1]-lensum.abscounts[2]) LT '+gmrmax+' AND '+$
                '(lensum.abscounts[1]-lensum.abscounts[2]) GE '+gmrmin+' AND '+$
                '(lensum.abscounts[1]-lensum.abscounts[2]) GE '+gmrcut+' AND '+$
                'lensum.absmag['+cstr+'] le '+highcut+' AND '+$
                'lensum.absmag['+cstr+'] gt '+lowcut+ ' AND '+$
                'lensum.z gt '+minz+' AND '+$
                'lensum.z lt '+maxz
                            
              addstr = 'redlum2threebin'
              print,'Sub sample: ',addstr
              print,lowcut+' < M['+!colors[clr]+'] < '+highcut
              zsub_sample, file, rfile, wstr, $
                           addstr=addstr,uselens=len,userand=rlen, $
                           overwrite=overwrite, indices=w2, $
                           outdir=outdir

              lowcut = ntostr(range_struct.red_minmag_threebin[clr,2])
              highcut = ntostr(range_struct.red_maxmag_threebin[clr,2])
              wstr = $
                '(lensum.abscounts[1]-lensum.abscounts[2]) LT '+gmrmax+' AND '+$
                '(lensum.abscounts[1]-lensum.abscounts[2]) GE '+gmrmin+' AND '+$
                '(lensum.abscounts[1]-lensum.abscounts[2]) GE '+gmrcut+' AND '+$
                'lensum.absmag['+cstr+'] le '+highcut+' AND '+$
                'lensum.absmag['+cstr+'] gt '+lowcut+ ' AND '+$
                'lensum.z gt '+minz+' AND '+$
                'lensum.z lt '+maxz
              
              addstr = 'redlum3threebin'
              print,'Sub sample: ',addstr
              print,lowcut+' < M['+!colors[clr]+'] < '+highcut
              zsub_sample, file, rfile, wstr, $
                           addstr=addstr,uselens=len,userand=rlen, $
                           overwrite=overwrite, indices=w2, $
                           outdir=outdir

          ENDFOR 
          delvarx, len, rlen

      END 
      'eclass': BEGIN 
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Just main galaxy sample with a cut in g-r
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          print
          print,'Reading: ',file
          len = mrdfits(file, 1)
          print,'Reading: ',rfile
          rlen = mrdfits(rfile, 1)

          wstr = $
            '(lensum.z LE '+maxz+') AND (lensum.z GE '+minz+')'
          wstr = wstr + $
            ' AND ( (lensum.primtarget and 2L^6) NE 0)'
          wstr = wstr + $
            ' AND ( lensum.eclass lt '+eclasscut+' )'

          addstr = 'eclass1'
          print
          print,'Doing Sub Sample: ',addstr
          zsub_sample, file, rfile, wstr, $
                       addstr=addstr,uselens=len,userand=rlen, $
                       overwrite=overwrite
          wstr = $
            '(lensum.z LE '+maxz+') AND (lensum.z GE '+minz+')'
          wstr = wstr + $
            ' AND ( (lensum.primtarget and 2L^6) NE 0)'
          wstr = wstr + $
            ' AND ( lensum.eclass ge '+eclasscut+' )'

          addstr = 'eclass2'
          print
          print,'Doing Sub Sample: ',addstr
          zsub_sample, file, rfile, wstr, $
                       addstr=addstr,uselens=len,userand=rlen, $
                       overwrite=overwrite
          delvarx, len, rlen

      END 
      'earlylumtwobin': BEGIN 

          read_cuts, range_struct
          
          FOR clr=0, 4 DO BEGIN 

              cstr = ntostr(clr)
              outdir = indir + 'sublum/'+!colors[clr]+'/'

              lowcut = ntostr(range_struct.early_minmag_twobin[clr,0])
              highcut = ntostr(range_struct.early_maxmag_twobin[clr,0])
              wstr = $
                'lensum.eclass lt '+eclasscut+' AND '+$
                'lensum.absmag['+cstr+'] le '+highcut+' AND '+$
                'lensum.absmag['+cstr+'] gt '+lowcut+ ' AND '+$
                'lensum.z gt '+minz+' AND '+$
                'lensum.z lt '+maxz
              
              addstr = 'earlylum1twobin'
              print,'Sub sample: ',addstr
              print,lowcut+' < M['+!colors[clr]+'] < '+highcut
              zsub_sample, file, rfile, wstr, $
                           addstr=addstr,uselens=len,userand=rlen, $
                           overwrite=overwrite, indices=w1, $
                           outdir=outdir
              
              lowcut = ntostr(range_struct.early_minmag_twobin[clr,1])
              highcut = ntostr(range_struct.early_maxmag_twobin[clr,1])
              wstr = $
                'lensum.eclass lt '+eclasscut+' AND '+$
                'lensum.absmag['+cstr+'] le '+highcut+' AND '+$
                'lensum.absmag['+cstr+'] gt '+lowcut+ ' AND '+$
                'lensum.z gt '+minz+' AND '+$
                'lensum.z lt '+maxz
              
              addstr = 'earlylum2twobin'
              print,'Sub sample: ',addstr
              print,lowcut+' < M['+!colors[clr]+'] < '+highcut
              zsub_sample, file, rfile, wstr, $
                           addstr=addstr,uselens=len,userand=rlen, $
                           overwrite=overwrite, indices=w2, $
                           outdir=outdir

          ENDFOR 
          delvarx, len, rlen

      END 
      'earlylumthreebin': BEGIN 

          read_cuts, range_struct
          
          FOR clr=0, 4 DO BEGIN 

              cstr = ntostr(clr)
              outdir = indir + 'sublum/'+!colors[clr]+'/'

              lowcut = ntostr(range_struct.early_minmag_threebin[clr,0])
              highcut = ntostr(range_struct.early_maxmag_threebin[clr,0])
              wstr = $
                'lensum.eclass lt '+eclasscut+' AND '+$
                'lensum.absmag['+cstr+'] le '+highcut+' AND '+$
                'lensum.absmag['+cstr+'] gt '+lowcut+ ' AND '+$
                'lensum.z gt '+minz+' AND '+$
                'lensum.z lt '+maxz
              
              addstr = 'earlylum1threebin'
              print,'Sub sample: ',addstr
              print,lowcut+' < M['+!colors[clr]+'] < '+highcut
              zsub_sample, file, rfile, wstr, $
                           addstr=addstr,uselens=len,userand=rlen, $
                           overwrite=overwrite, indices=w1, $
                           outdir=outdir
              
              lowcut = ntostr(range_struct.early_minmag_threebin[clr,1])
              highcut = ntostr(range_struct.early_maxmag_threebin[clr,1])
              wstr = $
                'lensum.eclass lt '+eclasscut+' AND '+$
                'lensum.absmag['+cstr+'] le '+highcut+' AND '+$
                'lensum.absmag['+cstr+'] gt '+lowcut+ ' AND '+$
                'lensum.z gt '+minz+' AND '+$
                'lensum.z lt '+maxz
              
              addstr = 'earlylum2threebin'
              print,'Sub sample: ',addstr
              print,lowcut+' < M['+!colors[clr]+'] < '+highcut
              zsub_sample, file, rfile, wstr, $
                           addstr=addstr,uselens=len,userand=rlen, $
                           overwrite=overwrite, indices=w2, $
                           outdir=outdir

              lowcut = ntostr(range_struct.early_minmag_threebin[clr,2])
              highcut = ntostr(range_struct.early_maxmag_threebin[clr,2])
              wstr = $
                'lensum.eclass lt '+eclasscut+' AND '+$
                'lensum.absmag['+cstr+'] le '+highcut+' AND '+$
                'lensum.absmag['+cstr+'] gt '+lowcut+ ' AND '+$
                'lensum.z gt '+minz+' AND '+$
                'lensum.z lt '+maxz
              
              addstr = 'earlylum3threebin'
              print,'Sub sample: ',addstr
              print,lowcut+' < M['+!colors[clr]+'] < '+highcut
              zsub_sample, file, rfile, wstr, $
                           addstr=addstr,uselens=len,userand=rlen, $
                           overwrite=overwrite, indices=w2, $
                           outdir=outdir

          ENDFOR 
          delvarx, len, rlen

      END 
      'vdis': BEGIN 
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Just main galaxy sample with a cut in g-r
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;          print
;          print,'Reading: ',file
;          len = mrdfits(file, 1)
;          print,'Reading: ',rfile
;          rlen = mrdfits(rfile, 1)

;          vcut = '232'
          
          vmax = '400'
          vmin = '50'
          vcut = '182'
          wstr = $
            '(lensum.z LE '+maxz+') AND (lensum.z GE '+minz+')'
          wstr = wstr + $
            ' AND ( lensum.vel_dis gt '+vmin+' )'
          wstr = wstr + $
            ' AND ( lensum.vel_dis lt '+vmax+' )'
          wstr = wstr + $
            ' AND ( lensum.vel_dis lt '+vcut+' )'

          addstr = 'vdis1'
          print
          print,'Doing Sub Sample: ',addstr
          zsub_sample, file, rfile, wstr, $
                       addstr=addstr,uselens=len,userand=rlen, $
                       overwrite=overwrite

          wstr = $
            '(lensum.z LE '+maxz+') AND (lensum.z GE '+minz+')'
          wstr = wstr + $
            ' AND ( lensum.vel_dis gt '+vmin+' )'
          wstr = wstr + $
            ' AND ( lensum.vel_dis lt '+vmax+' )'
          wstr = wstr + $
            ' AND ( lensum.vel_dis ge '+vcut+' )'

          addstr = 'vdis2'
          print
          print,'Doing Sub Sample: ',addstr
          zsub_sample, file, rfile, wstr, $
                       addstr=addstr,uselens=len,userand=rlen, $
                       overwrite=overwrite

          delvarx, len, rlen

      END 
      'test': BEGIN 
          dir = '~/lensout/stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_76_86/'
          lensumfile = dir + $
            'zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_76_86_gri_recorr_h_lensum_N1.fit'
          randsumfiles = dir + $
            'zrand'+ntostr(lindgen(10))+$
            '_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_76_86_gri_recorr_h_lensum_N1.fit'


          print,lensumfile
          print,randsumfiles

          wstr = '(lensum.z gt 0.15)'
          addstr = 'testing'
          ;;outdir = '~/tmp/'

          print
          print,'Doing Sub Sample: ',addstr
          zsub_sample_multirand, lensumfile, randsumfiles, wstr, $
                       addstr=addstr,uselens=len,userand=rlen, $
                       overwrite=overwrite, outdir=outdir

      END 
      ELSE: message,'Unknown type: '+ntostr(type)
  ENDCASE
 

END 
