PRO setzero, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, $
             v15, v16

  np=n_params()
  IF np EQ 0 THEN BEGIN 
      print,'-Syntax: setzero, v1, v2, ....'
      print,'sets variables to 0L.  Up to 16 arguments'
      return
  ENDIF 


  CASE np OF 
      1: v1=0L
      2:  BEGIN & v1=0L & v2=0L & END 
      3:  BEGIN & v1=0L & v2=0L & v3=0L & END 
      4:  BEGIN & v1=0L & v2=0L & v3=0L & v4=0L & END 
      5:  BEGIN & v1=0L & v2=0L & v3=0L & v4=0L & v5=0L & END 
      6:  BEGIN & v1=0L & v2=0L & v3=0L & v4=0L & v5=0L & v6=0L & END 
      7:  BEGIN & v1=0L & v2=0L & v3=0L & v4=0L & v5=0L & v6=0L & v7=0L & END 
      8:  BEGIN & v1=0L & v2=0L & v3=0L & v4=0L & v5=0L & v6=0L & v7=0L & v8=0L & END 
      9:  BEGIN & v1=0L & v2=0L & v3=0L & v4=0L & v5=0L & v6=0L & v7=0L & v8=0L & v9=0L & END
      10: BEGIN & v1=0L & v2=0L & v3=0L & v4=0L & v5=0L & v6=0L & v7=0L & v8=0L & v9=0L & v10=0L & END
      11: BEGIN & v1=0L & v2=0L & v3=0L & v4=0L & v5=0L & v6=0L & v7=0L & v8=0L & v9=0L & v10=0L & v11=0L & END
      12: BEGIN & v1=0L & v2=0L & v3=0L & v4=0L & v5=0L & v6=0L & v7=0L & v8=0L & v9=0L & v10=0L & v11=0L & v12=0L & END
      13: BEGIN & v1=0L & v2=0L & v3=0L & v4=0L & v5=0L & v6=0L & v7=0L & v8=0L & v9=0L & v10=0L & v11=0L & v12=0L & v13=0L & END
      14: BEGIN & v1=0L & v2=0L & v3=0L & v4=0L & v5=0L & v6=0L & v7=0L & v8=0L & v9=0L & v10=0L & v11=0L & v12=0L & v13=0L & v14=0L & END
      15: BEGIN & v1=0L & v2=0L & v3=0L & v4=0L & v5=0L & v6=0L & v7=0L & v8=0L & v9=0L & v10=0L & v11=0L & v12=0L & v13=0L & v14=0L & v15=0L & END
      16: BEGIN & v1=0L & v2=0L & v3=0L & v4=0L & v5=0L & v6=0L & v7=0L & v8=0L & v9=0L & v10=0L & v11=0L & v12=0L & v13=0L & v14=0L & v15=0L & v16=0L & END
      ELSE: message,'Too many parameters' 
  ENDCASE 
END 
