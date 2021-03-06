PRO setarrzero, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, $
                v15, v16

  np=n_params()
  IF np EQ 0 THEN BEGIN 
      print,'-Syntax: setarrzero, v1, v2, ....'
      print,'sets array elements to 0   Up to 16 arguments'
      print,'Must be number arrays'
      return
  ENDIF 


  CASE np OF 
      1:  v1[*]=0
      2:  BEGIN & v1[*]=0 & v2[*]=0 & END 
      3:  BEGIN & v1[*]=0 & v2[*]=0 & v3[*]=0 & END 
      4:  BEGIN & v1[*]=0 & v2[*]=0 & v3[*]=0 & v4[*]=0 & END 
      5:  BEGIN & v1[*]=0 & v2[*]=0 & v3[*]=0 & v4[*]=0 & v5[*]=0 & END 
      6:  BEGIN & v1[*]=0 & v2[*]=0 & v3[*]=0 & v4[*]=0 & v5[*]=0 & v6[*]=0 & END 
      7:  BEGIN & v1[*]=0 & v2[*]=0 & v3[*]=0 & v4[*]=0 & v5[*]=0 & v6[*]=0 & v7[*]=0 & END 
      8:  BEGIN & v1[*]=0 & v2[*]=0 & v3[*]=0 & v4[*]=0 & v5[*]=0 & v6[*]=0 & v7[*]=0 & v8[*]=0 & END 
      9:  BEGIN & v1[*]=0 & v2[*]=0 & v3[*]=0 & v4[*]=0 & v5[*]=0 & v6[*]=0 & v7[*]=0 & v8[*]=0 & v9[*]=0 & END
      10: BEGIN & v1[*]=0 & v2[*]=0 & v3[*]=0 & v4[*]=0 & v5[*]=0 & v6[*]=0 & v7[*]=0 & v8[*]=0 & v9[*]=0 & v10[*]=0 & END
      11: BEGIN & v1[*]=0 & v2[*]=0 & v3[*]=0 & v4[*]=0 & v5[*]=0 & v6[*]=0 & v7[*]=0 & v8[*]=0 & v9[*]=0 & v10[*]=0 & v11[*]=0 & END
      12: BEGIN & v1[*]=0 & v2[*]=0 & v3[*]=0 & v4[*]=0 & v5[*]=0 & v6[*]=0 & v7[*]=0 & v8[*]=0 & v9[*]=0 & v10[*]=0 & v11[*]=0 & v12[*]=0 & END
      13: BEGIN & v1[*]=0 & v2[*]=0 & v3[*]=0 & v4[*]=0 & v5[*]=0 & v6[*]=0 & v7[*]=0 & v8[*]=0 & v9[*]=0 & v10[*]=0 & v11[*]=0 & v12[*]=0 & v13[*]=0 & END
      14: BEGIN & v1[*]=0 & v2[*]=0 & v3[*]=0 & v4[*]=0 & v5[*]=0 & v6[*]=0 & v7[*]=0 & v8[*]=0 & v9[*]=0 & v10[*]=0 & v11[*]=0 & v12[*]=0 & v13[*]=0 & v14[*]=0 & END
      15: BEGIN & v1[*]=0 & v2[*]=0 & v3[*]=0 & v4[*]=0 & v5[*]=0 & v6[*]=0 & v7[*]=0 & v8[*]=0 & v9[*]=0 & v10[*]=0 & v11[*]=0 & v12[*]=0 & v13[*]=0 & v14[*]=0 & v15[*]=0 & END
      16: BEGIN & v1[*]=0 & v2[*]=0 & v3[*]=0 & v4[*]=0 & v5[*]=0 & v6[*]=0 & v7[*]=0 & v8[*]=0 & v9[*]=0 & v10[*]=0 & v11[*]=0 & v12[*]=0 & v13[*]=0 & v14[*]=0 & v15[*]=0 & v16[*]=0 & END
      ELSE: message,'Too many parameters' 
  ENDCASE 
END 
