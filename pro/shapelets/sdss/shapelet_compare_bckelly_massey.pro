PRO shapelet_compare_bckelly_massey_compare, original, recon_bckelly, recon_massey, skysig, idstring, dtype

  !p.multi=[0,2,2]

  nrk=n_elements(recon_bckelly)
  recon_bckelly[*] = recon_bckelly[*] + skysig*randomu(seed,nrk,/normal)

  nrm=n_elements(recon_massey)
  recon_massey[*] = recon_massey[*] + skysig*randomu(seed,nrm,/normal)


  siglevels = [3,5,10,20,50]
  c_lines=[0,0,0,0,0]

  IF (!d.name EQ 'Z') OR (!d.name EQ 'PS') THEN BEGIN 

      IF !d.name EQ 'Z' THEN BEGIN 
          setupplot,'Z'
          device, set_resolution=[770,890]
      ENDIF 
      loadct, 0
      c_colors = [!p.color,!p.color,!p.color,100,100]

  ENDIF ELSE IF !d.name EQ 'X' THEN BEGIN 
      setupplot,'X'
      c_colors = [!white, !green, !red, !blue, !darkGreen]
  ENDIF 

  levels = siglevels*skysig

  image_contour, original, levels=levels, c_colors=c_colors, $
    title=idstring, sky=0.0

  image_contour, recon_bckelly, levels=levels, c_colors=c_colors, $
    title='Kelly Reconstruction "'+dtype+'"', sky=0.0
  image_contour, recon_massey, levels=levels, c_colors=c_colors, $
    title='Massey Reconstruction "ls"', sky=0.0

  legend,ntostr(siglevels)+' '+!csym.sigma,$
    /right,box=0,charsize=0.7, lines=c_lines, colors=c_colors

  

END 


PRO shapelet_compare_bckelly_massey, str, psFieldInfo=psFieldInfo

  IF n_elements(str) EQ 0 THEN BEGIN 
      file = '~/shapelet_outputs/stripe16-vagc-shapelets-approx-nmax15.fit'
      str = mrdfits(file,1)
  ENDIF 

  nmax = 15
;  dtype = 'approx'
  dtype = 'pcr'

  window,0,xsize=800,ysize=890
;  window,1,xsize=800,ysize=512

  outDir = '~/shapelet_outputs/compare_bckelly_massey/'
  bcKellyDir = outDir + 'bckelly/'
  masseyDir = outDir + 'massey/'
  compareDir = outDir + 'compare/'

  outFileFront = 'stripe16-nmax' + ntostr(nmax) + '-'

  bcKellyTime = 0d
  MasseyTime = 0d

  nstr=n_elements(str)
  

  FOR nn=0L, nstr-1 DO BEGIN 

      !p.multi=[0,2,2]
;      wset,0

      idString = shapelet_idstring(str[nn])

      BCKellyFile = $
        bcKellyDir + outFileFront + idString+'-bckelly-'+dtype+'.png'
      MasseyFile  = $
        masseyDir  + outFileFront + idString+'-massey-ls.png'
      compareFile = $
        compareDir + outFileFront + idString+'-compare-ls-'+dtype+'.png'

      read_atlas, str[nn], imr=imr, col0=col0, row0=row0

      print
      print,'Calling Brandons code'
      print,BCKellyFile

      tm = systime(1)
      shapelet_decomp, nmax, decomp, str[nn], 2, recon=recon, $
        psFieldInfo=psFieldInfo, dtype=dtype
      tm = systime(1)-tm
      print,'BCKelly time'
      ptime,tm
      bcKellyTime = bcKellyTime + tm

      shapelet_view_recon, imr, recon, psFieldInfo, /addnoise
      write_png, BCKellyFile, tvrd(/true)
      
      w=where(imr EQ 0, nw)

      imrcopy = imr
      imrcopy[w] = psFieldInfo.skysig*randomu(seed, nw, /normal)
      beta = sqrt(str[nn].m_rr_cc[2]/2.0)
      x0 = [str[nn].colc[2],str[nn].rowc[2]] - float([col0[2], row0[2]])

      print
      print,'Calling Masseys code'
      print,MasseyFile
      tm = systime(1)
      shapelets_decomp, imrcopy, beta, nmax, decomp2, recomp2, $
        noise=psFieldInfo.skysig, x0=x0
      tm = systime(1)-tm
      print,'Massey time'
      ptime,tm
      MasseyTime = MasseyTime + tm

      shapelet_view_recon, imr, recomp2, psFieldInfo, /addnoise
      write_png, MasseyFile, tvrd(/true)

;      wset,1
      shapelet_compare_bckelly_massey_compare, $
        imr, recon, recomp2, psFieldInfo.skysig, idstring, dtype
      write_png, compareFile, tvrd(/true)

      key=get_kbrd(1)

  ENDFOR 
  print,'Total BCKelly time'
  ptime,BCKellyTime
  print,'Total Massey time'
  ptime,MasseyTime

END 
