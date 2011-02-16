FUNCTION run_shapelet_decomp_addextra, decomp, tstr, vagc_index

  new_decomp = $
    create_struct(decomp, $
                  'stripe', tstr.stripe, $
                  'ra',tstr.ra,$
                  'dec',tstr.dec,$
                  'z', tstr.z, $
                  'cmodel_counts', tstr.cmodel_counts, $
                  'm_e1', tstr.m_e1,$
                  'm_e2', tstr.m_e2, $
                  'm_e1_psf',tstr.m_e1_psf,$
                  'm_e2_psf',tstr.m_e2_psf,$
                  'vagc_index', vagc_index)

  return,new_decomp

END 



PRO run_shapelet_decomp, stripes, dtype, nmax, clr, $
                         doplot=doplot, png=png, zbuff=zbuff, prompt=prompt,$
                         status=status

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: run_shapelet_decomp, stripes, dtype, nmax, clr, /doplot, /png, /zbuff, /prompt, status='
      return
  ENDIF 

  on_error, 2

  shapelet = obj_new('sdss_shapelet')

  status = 1

  tmtot = systime(1)
  
  minmag = 15.0
  maxmag = 18.0
  maxz_vlim = 0.1
  max_abs_mag = -19.7

  outDir = shapelet->dir()

  dtype = strlowcase(dtype)

  colors = ['u','g','r','i','z']
  cstr = colors[clr]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set up some plotting variables
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(doplot) THEN BEGIN 
      title = 'Method: '+strupcase(dtype)+'   Nmax: '+ntostr(nmax)
      !p.multi=[0,2,2]
      
      siglevels = [3,5,10,20,50]
      c_lines=[0,0,0,0,0]
      
      IF keyword_set(zbuff) THEN BEGIN 
          setupplot,'Z'
          device, set_resolution=[770,890]
          
          loadct, 0
          
          c_colors = [!p.color,!p.color,!p.color,100,100]
      ENDIF ELSE BEGIN 
          setupplot,'X'
          c_colors = [!white, !green, !red, !blue, !darkGreen]
      ENDELSE 
  ENDIF 

  nStripe = n_elements(stripes)



  FOR si=0L, nStripe-1 DO BEGIN 

      stripe = stripes[si]
      sstr = stripe2string(stripe)

      outfile = shapelet->output_file(stripe, clr, dtype, nmax)

      query = $
        'select run,rerun,camcol,field,id,stripe,'+$
                'ra, dec, rowc, colc, m_rr_cc, m_rr_cc_psf, m_e1, m_e2,'+$
                'm_e1_psf, m_e2_psf, z, cmodel_counts, absmodelmag, '+$
                'match_rerun, match_id, flags, primtarget '+$
        'from specgal_main where stripe = '+ntostr(stripe)
      print,query


      str = pgsql_query(query)



      wg = where(str.match_rerun NE -1 AND $
                 str.m_rr_cc[clr] GT 0.0 AND $
                 str.cmodel_counts[2] LT maxmag AND $
                 str.cmodel_counts[2] GT minmag AND $
                 str.z LE maxz_vlim, nwg)

      make_flag_struct,fs
      fs.satur='N'
      fs.satur_center = 'N'
      fs.bright='N'
      fs.edge = 'N'
      flag_select, str, fs, 2, w, input_index=wg



      ;; Remove duplicates in the match space. Will also sort by id
      str.rerun = str.match_rerun
      str.id = str.match_id
      phid = sdss_photoid(str[w])
      rmd = rem_dup(phid)
      
      w = w[rmd]

      nprocess = n_elements(w)

      ;; Only do random?
      IF 0 THEN BEGIN 
          ;; first randomize
          rs = sort(randomu(seed, nprocess))
          nprocess = long(nprocess*0.154) > 1

          w = w[rs[0:nprocess-1]]

      ENDIF 

      print
      print,'Processing '+ntostr(nprocess)+' galaxies'
      print
      print,'Outfile: ',outfile

      tm=systime(1)

      first = 1
      FOR i=0L, nprocess-1 DO BEGIN 
          ii = w[i]
          
          tstr = str[ii]

          ;;tstr.rerun = ntostr(tstr.match_rerun)
          ;;tstr.id = tstr.match_id
          
          
          decomp = shapelet->atlas_decomp(tstr, dtype, nmax, clr, $
                                          /trim_atlas, $
                                          psFieldInfo=psFieldInfo, $
                                          recon=recon, $
                                          trimmed_image=trimmed_image, $
                                          status=dstatus)

;          shapelet_decomp,nmax,decomp,tstr,clr,dtype=dtype,$
;            recon=recon,psFieldInfo=psFieldInfo, $
;            trimmed_image=trimmed_image, /trim_atlas, $
;            status=status


          IF dstatus EQ 0 THEN BEGIN 
              decomp = run_shapelet_decomp_addextra(decomp, tstr, ii)

              IF first THEN BEGIN 
                  write_idlstruct, decomp, outFile, info_struct=info_struct
                  first = 0
              ENDIF ELSE BEGIN 
                  write_idlstruct, decomp, outFile, info_struct=info_struct, $
                    /append
              ENDELSE 
          ENDIF 




          IF keyword_set(doplot) AND dstatus EQ 0 THEN BEGIN 

              IF keyword_set(png) THEN BEGIN 
                  pngFile = shapelet->png_file(decomp,clr,/createdir)
                  print,'   Stripe: '+sstr+' '+pngFile
              ENDIF ELSE BEGIN 
                  idstring = shapelet->idstring(tstr)
                  print,'   Stripe: '+sstr+' '+idstring
              ENDELSE 




              IF size(psFieldInfo,/tname) EQ 'STRUCT' THEN BEGIN 
                  skysig = psFieldInfo.skysig
              ENDIF ELSE BEGIN 
                  skysig=5.0
              ENDELSE 
              
              
              view_atlas, tstr, /silent, clr=clr, imtot=image, /hideradec
              
              ;; Just plot blank if decomposition was not successful
              IF n_elements(recon) EQ 0 THEN BEGIN 
                  image_sz=size(image,/dim)
                  recon=intarr(image_sz[0], image_sz[1])
              ENDIF ELSE BEGIN 
                  nrecon=n_elements(recon)
                  recon[*]=recon[*] + skysig*randomu(seed,nrecon,/normal)
              ENDELSE 
              tvasinh, recon, title=title, sky=0.0
              
              IF dtype EQ 'approx' THEN BEGIN 
                  orthostr = $
                    'ortho = '+$
                    ntostr(decomp.ortho/1.e-7,5,/round) + $
                    !csym.times+'10!U'+!csym.minus+'7!N'
                  orthostr = $
                    'ortho = '+$
                    ntostr(decomp.ortho/1.e-7) + $
                    !csym.times+'10!U'+!csym.minus+'7!N'
                  legend,orthostr,box=0,charsize=1
              ENDIF 
              ;; plot contours at various sigma
              
              levels = siglevels*skysig
              ;;image_contour, image, levels=levels, c_colors=c_colors
              image_contour, trimmed_image, $
                levels=levels, c_colors=c_colors, sky=0.0
              image_contour, recon, $
                levels=levels, c_colors=c_colors, sky=0.0
              
              legend,ntostr(siglevels)+' '+!csym.sigma,$
                /right,box=0,charsize=0.7, lines=c_lines, colors=c_colors
              
              legend,'C_p = '+ntostr(decomp.c_p),/bottom,/left,box=0,charsize=1

              IF keyword_set(png) THEN BEGIN 
                  IF keyword_set(zbuff) THEN BEGIN 
                      write_png, pngFile, tvrd(), rmap, gmap, bmap
                  ENDIF ELSE BEGIN 
                      write_png, pngFile, tvrd(/true)
                  ENDELSE 
              ENDIF 

          ENDIF 
          

          IF keyword_set(prompt) THEN BEGIN 
              key=prompt_kbrd("Hit a key")
              IF key EQ 'q' THEN BEGIN 
                  obj_destroy, shapelet
                  return
              ENDIF 
          ENDIF 

      ENDFOR 

      print,'Finished: ',outfile
      print,'For stripe: '+sstr
      ptime,systime(1)-tm

  ENDFOR ;; Loop over stripes

  IF keyword_set(zbuff) THEN setupplot,'X'

  IF nStripe GT 1 THEN BEGIN 
      print
      print,'Total time'
      ptime,systime(1)-tmtot
  ENDIF 

  !p.multi=0

  status = 0
  obj_destroy, shapelet

END 
