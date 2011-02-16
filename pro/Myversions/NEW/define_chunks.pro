PRO define_chunks

  chunknums = [3,4,5,7,8,9]

  defsysv, '!CHUNKNUMS', chunknums

  defsysv, '!CHUNK3', [259,273]
  defsysv, '!CHUNK4', [752,756]
  defsysv, '!CHUNK5', [1140,1231]
  defsysv, '!CHUNK7', [1336,1339, 1356,1359]
  defsysv, '!CHUNK8', [94, 125, 1033, 1056]
  defsysv, '!CHUNK9', [1035, 1043]

  message,'CHUNK system variables have been added',/inf


  return
END 
