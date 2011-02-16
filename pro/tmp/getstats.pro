PRO getstats, str, int, dex, wis, cha, ntrial=ntrial

COMMON seed,seed

IF n_elements(ntrial) EQ 0 THEN ntrial = 1

str = intarr(3)
int = intarr(3)
dex = intarr(3)
wis = intarr(3)
cha = intarr(3)

FOR i=0, 2 DO BEGIN

    str[i] = max(round(arrscl(randomu(seed,ntrial),1,6,arrmin=0.,arrmax=1.)))
    int[i] = max(round(arrscl(randomu(seed,ntrial),1,6,arrmin=0.,arrmax=1.)))
    dex[i] = max(round(arrscl(randomu(seed,ntrial),1,6,arrmin=0.,arrmax=1.)))
    wis[i] = max(round(arrscl(randomu(seed,ntrial),1,6,arrmin=0.,arrmax=1.)))
    cha[i] = max(round(arrscl(randomu(seed,ntrial),1,6,arrmin=0.,arrmax=1.)))

ENDFOR

str = fix( total(str) )
int = fix( total(int) )
dex = fix( total(dex) )
wis = fix( total(wis) )
cha = fix( total(cha) )

print,' str: ',str
print,' int: ',int
print,' dex: ',dex
print,' wis: ',wis
print,' cha: ',cha

return
END  
