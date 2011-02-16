FUNCTION mysql_struct2coldefs, struct, max_string_length=max_string_length,$
                               sdss=sdss, $
                               outtags=outtags

  ;; Generate a mysql create table statement from the input
  ;; structure. If /sdss, then 5-elements arrays are assumede to
  ;; be u,g,r,i,z object measurements. In this case, the mysql 
  ;; columns are appended _u, _g, _r, _i, _z rather than the default
  ;; of _1, _2, _3, ...

  ;; 

  IF n_elements(max_string_length) EQ 0 THEN max_string_length = 255

  tags = strlowcase( tag_names(struct) )
  ntags = n_elements(tags)

  sdss_bands = ['u','g','r','i','z']

  FOR i=0L, ntags-1 DO BEGIN 

      tmpvar = struct[0].(i)
      nelem = n_elements(tmpvar)

      tn = size(tmpvar,/tname)

      IF tn EQ 'STRING' THEN BEGIN 

          maxlen = max( strlen(struct.(i)) ) > max_string_length 
          mysql_type = mysql_idl2mysql('STRING',length=maxlen)

      ENDIF ELSE BEGIN 
          mysql_type = mysql_idl2mysql(tn)
      ENDELSE 

      print,tags[i],' '+tn

      IF nelem GT 1 THEN BEGIN 

          IF nelem EQ 5 AND keyword_set(sdss) THEN BEGIN 
              FOR j=0L, nelem-1 DO BEGIN 
                  
                  tag = tags[i]+'_'+sdss_bands[j]
                  coldef = $
                    tag +' '+mysql_type+' NOT NULL'
                  print,' |--> '+coldef
                  
                  add_arrval, coldef, coldefs
                  add_arrval, tag, outtags
              ENDFOR 
          ENDIF ELSE BEGIN 

              FOR j=0L, nelem-1 DO BEGIN 

                  tag = tags[i]+'_'+strtrim(j,2)
                  coldef = $
                    tag+' '+mysql_type+' NOT NULL'
                  print,' |--> '+coldef

                  add_arrval, coldef, coldefs
                  add_arrval, tag, outtags
              ENDFOR 

          ENDELSE 
      ENDIF ELSE BEGIN 
          coldef = $
            tags[i]+' '+mysql_type+' NOT NULL'
          
          print,' |--> '+coldef

          add_arrval, coldef, coldefs
          add_arrval, tags[i], outtags
      ENDELSE 

  ENDFOR 

  return,coldefs

END 
