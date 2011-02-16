PRO getseq,word

word=''
tmp=get_kbrd(1)
i=0
WHILE 1 DO BEGIN 
    key=get_kbrd(0)
    IF key EQ '' THEN return 
    word=word+key
ENDWHILE 

return
END 
