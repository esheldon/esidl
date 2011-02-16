PRO print_matrix, m, lun=lun, file=file, format=format

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: print_matrix, matrix, lun=lun, file=file, format=format'
      print
      print,'Max in x-direction is 18 for now (colprint)'
      return
  ENDIF 

  ndim = size(m, /n_dim)
  IF ndim NE 2 THEN message,'Matrix must be 2 dimensional'

  sz = size(m)
  nx = sz[1]
  ny = sz[2]

  CASE nx OF
      2:  colprint,m[0,*], lun=lun, file=file, format=format
      3:  colprint,m[0,*],m[1,*], lun=lun, file=file, format=format
      4:  colprint,m[0,*],m[1,*],m[2,*], lun=lun, file=file, format=format
      5:  colprint,m[0,*],m[1,*],m[2,*],m[3,*], lun=lun, file=file, format=format
      6:  colprint,m[0,*],m[1,*],m[2,*],m[3,*],m[4,*], lun=lun, file=file, format=format
      7:  colprint,m[0,*],m[1,*],m[2,*],m[3,*],m[4,*],m[5,*], lun=lun, file=file, format=format
      8:  colprint,m[0,*],m[1,*],m[2,*],m[3,*],m[4,*],m[5,*],m[6,*], lun=lun, file=file, format=format
      9:  colprint,m[0,*],m[1,*],m[2,*],m[3,*],m[4,*],m[5,*],m[6,*],m[7,*],m[8,*], lun=lun, file=file, format=format
      10:  colprint,m[0,*],m[1,*],m[2,*],m[3,*],m[4,*],m[5,*],m[6,*],m[7,*],m[8,*],m[9,*], lun=lun, file=file, format=format
      11: colprint,m[0,*],m[1,*],m[2,*],m[3,*],m[4,*],m[5,*],m[6,*],m[7,*],m[8,*],m[9,*],m[10,*], lun=lun, file=file, format=format
      12: colprint,m[0,*],m[1,*],m[2,*],m[3,*],m[4,*],m[5,*],m[6,*],m[7,*],m[8,*],m[9,*],m[10,*],m[11,*], lun=lun, file=file, format=format
      13: colprint,m[0,*],m[1,*],m[2,*],m[3,*],m[4,*],m[5,*],m[6,*],m[7,*],m[8,*],m[9,*],m[10,*],m[11,*],m[12,*], lun=lun, file=file, format=format
      14: colprint,m[0,*],m[1,*],m[2,*],m[3,*],m[4,*],m[5,*],m[6,*],m[7,*],m[8,*],m[9,*],m[10,*],m[11,*],m[12,*],m[13,*], lun=lun, file=file, format=format
      15: colprint,m[0,*],m[1,*],m[2,*],m[3,*],m[4,*],m[5,*],m[6,*],m[7,*],m[8,*],m[9,*],m[10,*],m[11,*],m[12,*],m[13,*],m[14,*], lun=lun, file=file, format=format
      16: colprint,m[0,*],m[1,*],m[2,*],m[3,*],m[4,*],m[5,*],m[6,*],m[7,*],m[8,*],m[9,*],m[10,*],m[11,*],m[12,*],m[13,*],m[14,*],m[15,*], lun=lun, file=file, format=format
      17: colprint,m[0,*],m[1,*],m[2,*],m[3,*],m[4,*],m[5,*],m[6,*],m[7,*],m[8,*],m[9,*],m[10,*],m[11,*],m[12,*],m[13,*],m[14,*],m[15,*],m[16,*], lun=lun, file=file, format=format
      18: colprint,m[0,*],m[1,*],m[2,*],m[3,*],m[4,*],m[5,*],m[6,*],m[7,*],m[8,*],m[9,*],m[10,*],m[11,*],m[12,*],m[13,*],m[14,*],m[15,*],m[16,*],m[17,*], lun=lun, file=file, format=format
      ELSE: message,'Max allowed in x-direction is 18'
  ENDCASE 

END 
