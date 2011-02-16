PRO make_lumdist

  ext='N1.fit'
  massclr = 2

  nclr = 5
  basedir =  '/sdss5/data0/lensout/'
  outdir = basedir+'luminosities/'
  outfile = outdir + 'lumdist_data_'+!colors[massclr]+'lensing.fit'
  print
  print,'Output file will be here: ',outfile
  print

  stripes = [10, 36, 37, 42, 43, 82]
  stripestr = 'stripe'+ntostr(stripes)
  nlum=4

  nst=n_elements(stripes)

  FOR i=0L, nst-1 DO BEGIN 
      
      FOR clr=0L, nclr-1 DO BEGIN 
          indir = basedir+stripestr[i]+'/sublum/'+!colors[clr]+'/'
          FOR lumi=1, nlum DO BEGIN 
              
              file=indir+'lum'+ntostr(lumi)+'_zgal_gal_'+stripestr[i]+'_'+$
                !colors[massclr]+'_lensum_N1.fit'
              print,'File: ',file
              tlum=mrdfits(file,1)

              val = !colors[clr]+'lum'+ntostr(lumi)
              command = 'add_arrval, tlum.lum[clr], '+val
              IF NOT execute(command) THEN message,'Error'

              val = !colors[clr]+'absmag'+ntostr(lumi)
              command = 'add_arrval, tlum.absmag[clr], '+val
              IF NOT execute(command) THEN message,'Error'

              tlum=0
          ENDFOR 
      ENDFOR 
      
  ENDFOR 

  lums=create_struct('ulum1',ulum1,$
                   'ulum2',ulum2,$
                   'ulum3',ulum3,$
                   'ulum4',ulum4,$
                   'glum1',glum1,$
                   'glum2',glum2,$
                   'glum3',glum3,$
                   'glum4',glum4,$
                   'rlum1',rlum1,$
                   'rlum2',rlum2,$
                   'rlum3',rlum3,$
                   'rlum4',rlum4,$
                   'ilum1',ilum1,$
                   'ilum2',ilum2,$
                   'ilum3',ilum3,$
                   'ilum4',ilum4,$
                   'zlum1',zlum1,$
                   'zlum2',zlum2,$
                   'zlum3',zlum3,$
                   'zlum4',zlum4 $
                   )
  abss=create_struct('uabsmag1',uabsmag1,$
                     'uabsmag2',uabsmag2,$
                     'uabsmag3',uabsmag3,$
                     'uabsmag4',uabsmag4,$
                     'gabsmag1',gabsmag1,$
                     'gabsmag2',gabsmag2,$
                     'gabsmag3',gabsmag3,$
                     'gabsmag4',gabsmag4,$
                     'rabsmag1',rabsmag1,$
                     'rabsmag2',rabsmag2,$
                     'rabsmag3',rabsmag3,$
                     'rabsmag4',rabsmag4,$
                     'iabsmag1',iabsmag1,$
                     'iabsmag2',iabsmag2,$
                     'iabsmag3',iabsmag3,$
                     'iabsmag4',iabsmag4,$
                     'zabsmag1',zabsmag1,$
                     'zabsmag2',zabsmag2,$
                     'zabsmag3',zabsmag3,$
                     'zabsmag4',zabsmag4 $
                    )
  
  ts = create_struct(lums, abss)

  print
  print,'Outputting file: ',outfile
  print
  mwrfits, ts, outfile, /create

END 
