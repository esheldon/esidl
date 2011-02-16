PRO zobjshear_names, prename, datfile, sumfile, zfile, lensumfile, groupfile, $
                     indfile, psfile, asciifile, extno=extno

  IF arg_present(asciifile) THEN doascii=1 ELSE doascii=0

  IF n_elements(extno) NE 0 THEN BEGIN 

      ext_str = 'N'+ntostr(long(extno))

      datfile = prename + '_'+ext_str+'.fit'
      sumfile = prename + '_sum_'+ext_str+'.fit'
      zfile   = prename + '_z_'+ext_str+'.fit'
      lensumfile = prename + '_lensum_'+ext_str+'.fit'
      groupfile = prename +'_groupstats_'+ext_str+'.fit'
      indfile = prename + '_ind_'+ext_str+'.dat'
      psfile = prename + '_'+ext_str+'.ps'
      asciifile = prename + '_lensum_'+ext_str+'.dat'

  ENDIF ELSE BEGIN 
  
      datfile = prename + '_N1.fit'
      sumfile = prename + '_sum_N1.fit'
      zfile   = prename + '_z_N1.fit'
      lensumfile = prename + '_lensum_N1.fit'
      groupfile = prename +'_groupstats_N1.fit'
      indfile = prename + '_ind_N1.dat'
      psfile = prename + '_N1.ps'
      asciifile = prename + '_lensum_N1.dat'

      WHILE fexist(datfile) OR fexist(psfile) DO BEGIN
          datfile = newname(datfile)
          sumfile = newname(sumfile)
          zfile = newname(zfile)
          lensumfile = newname(lensumfile)
          groupfile = newname(groupfile)
          indfile = newname(indfile)
          psfile = newname(psfile)
          asciifile = newname(asciifile)
      ENDWHILE 

  ENDELSE 
  ;;  print
  ;;  print,'Dat file: ',datfile
  ;;  print,'Sum file: ',sumfile
  ;;  print,'Lens sum file: ',lensumfile
  ;;  print,'Redshift file: ',zfile
  ;;  print,'Group stats file: ',groupfile
  ;;  print,'Index file: ',indfile
  ;;  print
      
  print
  echo,['Dat file: ',datfile],color=['green','none'],bold=[1,0],no=[1,0]
  ;;  echo,['Sum file: ',sumfile],color=['green','none'],bold=[1,0],no=[1,0]
  echo,['Lens sum file: ',lensumfile],color=['green','none'],bold=[1,0],no=[1,0]
  echo,['Redshift file: ',zfile],color=['green','none'],bold=[1,0],no=[1,0]
  echo,['Group stats file: ',groupfile],color=['green','none'],bold=[1,0],no=[1,0]
  echo,['psfile: ',psfile],color=['green','none'],bold=[1,0],no=[1,0]
  IF doascii THEN echo,['asciifile: ',asciifile],color=['green','none'],bold=[1,0],no=[1,0]
  ;;  echo,['Index file: ',indfile],color=['green','none'],bold=[1,0],no=[1,0]
  print

END 
