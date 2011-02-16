pro plt_wrap_moms2n


read_dir = '/sdss4/data1/esheldon/'


aratio = '95'

mname = read_dir + 'mids2nexptot'+aratio+'.fits'
mtot = mrdfits(mname,1,hdr)
w = where(mtot.s2n le 80.0)
plt_all_moms2n, mtot(w), aratio

aratio = '9'

mname = read_dir + 'mids2nexptot'+aratio+'.fits'
mtot = mrdfits(mname,1,hdr)
w = where(mtot.s2n le 80.0)
plt_all_moms2n, mtot(w), aratio

aratio = '8'

mname = read_dir + 'mids2nexptot'+aratio+'.fits'
mtot = mrdfits(mname,1,hdr)
w = where(mtot.s2n le 80.0)
plt_all_moms2n, mtot(w), aratio

aratio = '6'

mname = read_dir + 'mids2nexptot'+aratio+'.fits'
mtot = mrdfits(mname,1,hdr)
w = where(mtot.s2n le 80.0)
plt_all_moms2n, mtot(w), aratio

aratio = '4'

mname = read_dir + 'mids2nexptot'+aratio+'.fits'
mtot = mrdfits(mname,1,hdr)
w = where(mtot.s2n le 80.0)
plt_all_moms2n, mtot(w), aratio

aratio = '3'

mname = read_dir + 'mids2nexptot'+aratio+'.fits'
mtot = mrdfits(mname,1,hdr)
w = where(mtot.s2n le 80.0)
plt_all_moms2n, mtot(w), aratio

return
end



