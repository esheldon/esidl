FUNCTION daycnv_str, mjd

  months = ['Jan','Feb','Mar','Apr','May','Jun','Jul',$
            'Aug','Sep','Oct','Nov','Dec']

  n_mjd = n_elements(mjd)

  IF n_mjd EQ 1 THEN date_string = '' ELSE date_string=strarr(n_mjd)

  FOR i=0L, n_mjd-1 DO BEGIN 

      daycnv,(mjd[i])+(2400000.5d), yr, mn, day

      date_string[i] = $
        strn(day,len=2,padchar='0') + '-'+$
        months[mn-1]+'-'+strn(yr)
  ENDFOR 

  return,date_string

END 
