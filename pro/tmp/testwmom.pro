PRO testwmom

COMMON seed, seed
val = 10.0
n = 10000L
ntest = 20

index=lindgen(ntest*n)
t=dblarr(ntest*n)
err=t
v=t

wmean=dblarr(ntest)
wsig=wmean
werr=wmean

FOR i=0, ntest-1 DO BEGIN 

    ii = index[i*n:(i+1)*n -1 ]
    err[ii] = randomn(seed, n)
    t[ii] = replicate(val, n)+err[ii]
    v = ( moment(t[ii]) )[1]
    
    wmom, t[ii], sqrt(v + err[ii]^2), wm, ws, we
    wmean[i] = wm
    wsig[i] = ws
    werr[i] = we

ENDFOR 
ploterr, wmean, werr, yrange=[9.8, 10.2],psym=7

v = ( moment(t) )[1] 
print,v
wmom, t, sqrt( v + err^2), wm, ws, we
print,wm, ws, we

wmom, wmean, sqrt( wsig^2 + werr^2 ), wm, ws, we

print
print,wm, ws, we, median( wsig )

;; What I learned:
;; Use median of sig from bins to find total sig, but can use error
;; from weighted scheme as above.

return
END 
