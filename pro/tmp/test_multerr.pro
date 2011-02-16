PRO test_multerr

COMMON seed,seed
n=50
x=findgen(n)
a = replicate( 20., n)
b = a
aerr = replicate(1., n)
berr = aerr

a[*] = a[*] + randomn(seed, n)
b[*] = b[*] + randomn(seed, n)

aberr = sqrt( a^2*aerr^2 + b^2*berr^2 )

ploterr, x, a, aerr, psym=1
key=get_kbrd(1)
ploterr, x, b, berr, psym=1

key=get_kbrd(1)
ploterr, x, a*b, aberr, psym=1

wmom, a*b, aberr, wmean, wsig, wmerr

print, 'Mean = ',ntostr(wmean),' +/- ',ntostr(wmerr)

return
END  
