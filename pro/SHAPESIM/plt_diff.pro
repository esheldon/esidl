pro plt_diff,lup_med,ad_med


read_dir = '/sdss4/data1/esheldon/'
print_dir = '/sdss4/data1/esheldon/PS/Exp/DIFF/'

ad_diff = fltarr(6)
lup_diff = fltarr(6)
e1_all = fltarr(6)

aratio='95'
find_e1,.95,e1
fname = read_dir + 'mids2nexptot'+aratio+'.fits'
tot = mrdfits(fname,1,hdr)
w = where(tot.s2n le 35 and tot.s2n ge 30)
lup_med = median(tot[w].e1_lupton)
ad_med = median(tot[w].e1_ad)
e1_all[0] = e1
lup_diff[0] = lup_med - e1
ad_diff[0] = ad_med-e1


aratio='9'
find_e1,.9,e1
fname = read_dir + 'mids2nexptot'+aratio+'.fits'
tot = mrdfits(fname,1,hdr)
w = where(tot.s2n le 35 and tot.s2n ge 30)
lup_med = median(tot[w].e1_lupton)
ad_med = median(tot[w].e1_ad)
e1_all[1] = e1
lup_diff[1] = lup_med - e1
ad_diff[1] = ad_med-e1


aratio='8'
find_e1,.8,e1
fname = read_dir + 'mids2nexptot'+aratio+'.fits'
tot = mrdfits(fname,1,hdr)
w = where(tot.s2n le 35 and tot.s2n ge 30)
lup_med = median(tot[w].e1_lupton)
ad_med = median(tot[w].e1_ad)
e1_all[2] = e1
lup_diff[2] = lup_med - e1
ad_diff[2] = ad_med-e1


aratio='6'
find_e1,.6,e1
fname = read_dir + 'mids2nexptot'+aratio+'.fits'
tot = mrdfits(fname,1,hdr)
w = where(tot.s2n le 35 and tot.s2n ge 30)
lup_med = median(tot[w].e1_lupton)
ad_med = median(tot[w].e1_ad)
e1_all[3] = e1
lup_diff[3] = lup_med - e1
ad_diff[3] = ad_med-e1


aratio='4'
find_e1,.4,e1
fname = read_dir + 'mids2nexptot'+aratio+'.fits'
tot = mrdfits(fname,1,hdr)
w = where(tot.s2n le 35 and tot.s2n ge 30)
lup_med = median(tot[w].e1_lupton)
ad_med = median(tot[w].e1_ad)
e1_all[4] = e1
lup_diff[4] = lup_med - e1
ad_diff[4] = ad_med-e1

aratio='3'
find_e1,.3,e1
fname = read_dir + 'mids2nexptot'+aratio+'.fits'
tot = mrdfits(fname,1,hdr)
w = where(tot.s2n le 35 and tot.s2n ge 30)
lup_med = median(tot[w].e1_lupton)
ad_med = median(tot[w].e1_ad)
e1_all[5] = e1
lup_diff[5] = lup_med - e1
ad_diff[5] = ad_med-e1

print_name = print_dir + 'ad_lup_diff.ps'
title = 'Difference from input ellipticity'
xtitle = 'input e1'
ytitle = 'Deviation from input'
begplot, name = print_name, /invbw,/landscape
plot, e1_all, lup_diff,psym=1,yrange=[-.4,.2],title=title,xtitle=xtitle,$
	ytitle=ytitle
oplot, e1_all,ad_diff,psym=4
legend,['Adaptive','R^-2'],psym=[4,1]
endplot,/noprint

return
end









