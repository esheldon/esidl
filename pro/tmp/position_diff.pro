PRO position_diff, cdiff, cdev, rdiff, rdev

colmin = 1
colmax = 6

cdiff = fltarr(6, 5)
rdiff = fltarr(6, 5)
cdev = fltarr(6, 5)
rdev = fltarr(6, 5)

FOR col=colmin, colmax DO BEGIN 

    str = 0
    fetch_dir, 756, col, 1, dir
    read_tsobj, dir, str, start=200, nf=50

    nn = n_elements(str)

    ucdiff = str.colc[2] - str.colc[0]
    gcdiff = str.colc[2] - str.colc[1]
    icdiff = str.colc[2] - str.colc[3]
    zcdiff = str.colc[2] - str.colc[4]

    urdiff = str.rowc[2] - str.rowc[0]
    grdiff = str.rowc[2] - str.rowc[1]
    irdiff = str.rowc[2] - str.rowc[3]
    zrdiff = str.rowc[2] - str.rowc[4]

    cdiff[col-1, 0] = median(ucdiff)
    cdiff[col-1, 1] = median(gcdiff)
    cdiff[col-1, 3] = median(icdiff)
    cdiff[col-1, 4] = median(zcdiff)

    rdiff[col-1, 0] = median(urdiff)
    rdiff[col-1, 1] = median(grdiff)
    rdiff[col-1, 3] = median(irdiff)
    rdiff[col-1, 4] = median(zrdiff)

    cdev[col-1, 0] = sdev(ucdiff)
    cdev[col-1, 1] = sdev(gcdiff)
    cdev[col-1, 3] = sdev(icdiff)
    cdev[col-1, 4] = sdev(zcdiff)

    rdev[col-1, 0] = sdev(urdiff)
    rdev[col-1, 1] = sdev(grdiff)
    rdev[col-1, 3] = sdev(irdiff)
    rdev[col-1, 4] = sdev(zrdiff)

    print,'--------------------------------------'
    print,'Column ',ntostr(col)
    print
    print,'u col diff = ',cdiff[col-1, 0],' +/- ',cdev[col-1, 0]
    print,'g col diff = ',cdiff[col-1, 1],' +/- ',cdev[col-1, 1]
    print,'i col diff = ',cdiff[col-1, 3],' +/- ',cdev[col-1, 3]
    print,'z col diff = ',cdiff[col-1, 4],' +/- ',cdev[col-1, 4]
    print
    print,'u row diff = ',rdiff[col-1, 0],' +/- ',rdev[col-1, 0]
    print,'g row diff = ',rdiff[col-1, 1],' +/- ',rdev[col-1, 1]
    print,'i row diff = ',rdiff[col-1, 3],' +/- ',rdev[col-1, 3]
    print,'z row diff = ',rdiff[col-1, 4],' +/- ',rdev[col-1, 4]
    print
    print,'--------------------------------------'

ENDFOR 
return
END 
