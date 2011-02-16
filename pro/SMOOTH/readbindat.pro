PRO readbindat, file, nlines, binarea, e1err

readcol, file, v1, v2, v3, v4, v5,/silent

binarea = v1[0:nlines-1]

e1err    = v2[0:nlines-1]
e1errerr = v3[0:nlines-1]
e2err    = v4[0:nlines-1]
e2errerr = v5[0:nlines-1]

e1sig    = v2[nlines:2*nlines-1]
e1sigerr = v3[nlines:2*nlines-1]
e2sig    = v4[nlines:2*nlines-1]
e2sigerr = v5[nlines:2*nlines-1]

e1       = v2[2*nlines:3*nlines-1]
e1uncert = v3[2*nlines:3*nlines-1]
e2       = v4[2*nlines:3*nlines-1]
e2uncert = v5[2*nlines:3*nlines-1]

nobj     = v2[3*nlines:4*nlines-1]
theo     = v3[3*nlines:4*nlines-1]
ntot     = v4[3*nlines:4*nlines-1]

xtitle = 'Square Degrees'
ytitle = ''

ploterr, binarea, e1err, e1errerr, /xlog, /ylog, psym=1, xtitle=xtitle
oploterr, binarea, e1sig, e1sigerr, psym=3
oplot, binarea, theo*.43/.32

legend, ['<e1err>', 'bin2bin spread','.43/sqrt(N)'],psym=[1,3,0],$
          position=[.6,.06]


return
END 


