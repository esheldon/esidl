PRO mean_error_legend, names, means, errorlow, errorhigh, $
                       nkeep=nkeep,_extra=_extra

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax:  mean_error_legend, names, means, errorlow [, errorhigh, nkeep=nkeep, _extra=_extra]'
      return
  ENDIF 

  nn=n_elements(means)

  IF n_params() EQ 4 THEN BEGIN 

      FOR i=0L, nn-1 DO BEGIN 
          
          IF n_elements(nkeep) NE 0 THEN BEGIN               
              tms=ntostr(means[i],nkeep[i],/round)
              tel=ntostr(errorlow[i],nkeep[i],/round)
              teh=ntostr(errorhigh[i],nkeep[i],/round)
          ENDIF ELSE BEGIN 
              tms=ntostr(means[i])
              tel=ntostr(errorlow[i])
              teh=ntostr(errorhigh[i])
          ENDELSE 
          
          add_arrval, names[i] +': '+ tms +'!S!U+'+teh+'!R!D'+!csym.minus+tel,$
                      mess

          IF i LT nn-1 THEN add_arrval, '', mess
      ENDFOR 

  ENDIF ELSE IF n_params() EQ 3 THEN BEGIN 

      FOR i=0L, nn-1 DO BEGIN 
          
          IF n_elements(nkeep) NE 0 THEN BEGIN 
              add_arrval, ntostr(means[i],nkeep[i],/round), meanstring
              add_arrval, ntostr(errorlow[i],nkeep[i],/round), errorlowstring
          ENDIF ELSE BEGIN 
              add_arrval, ntostr(means[i]), meanstring
              add_arrval, ntostr(errorlow[i]), errorlowstring
          ENDELSE 
      ENDFOR 

      mess = names +': '+ meanstring + !csym.plusminus + errorlowstring+'!N'

  ENDIF ELSE BEGIN 
      print,'-Syntax:  mean_error_legend, means, errorlow [, errorhigh, digit=digit, nkeep=nkeep, _extra=_extra]'
      return
  ENDELSE 

  legend, mess, _extra=_extra

END 
