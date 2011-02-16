pro plt_moms2n,tot,w_type,aratio

if n_params() eq 0 then begin
	print,'-syntax  plt_moms2n,tot,w_type,aratio   aratio is a string'
	return
endif

w = where(tot[*].s2n le 400 and tot[*].e1_ad ge -2)
tot2 = tot[w]
s= sort(tot2[*].s2n)

;making bigger bins for this calculation
n = ((n_elements(tot2))/100)*100
print, 'HERE IS N:  ',n
print, 'HERE IS N/100: ',n/100

s2n = tot2[s].s2n

e1_ad = tot2[s].e1_ad
if w_type eq 'lupton' then begin
	e1 = tot2[s].e1_lupton*2
endif
if w_type eq 'gaussian' then begin
	e1 = tot2[s].e1_gaussian
endif
if w_type eq 'cowboy' then begin
	e1 = tot2[s].e1_cowboy
endif
if w_type eq 'none' then begin
	e1 = tot2[s].e1
endif

s2n = s2n[1:n]
e1 = e1[1:n+1]
s2n = rebin(s2n,n/100)

ad_s2n = fltarr(n/100)
e1_s2n = fltarr(n/100)
for i=0L,((n/100) -1) do begin
	var = moment( e1[i*100:(i+1)*100] )
	var_ad = moment(e1_ad[i*100:(i+1)*100])
	e1_s2n[i] = var[0]/sqrt(var[1])
	ad_s2n[i] = var_ad[0]/sqrt(var_ad[1])
endfor

e1_s2n = smooth(e1_s2n,3)
ad_s2n = smooth(ad_s2n,3)

print_ratio = aratio
if (print_ratio ne '1') then begin
	print_ratio = '.'+print_ratio
endif

if (w_type eq 'lupton') then begin
	w_type = 'R^-2'
endif
!p.multi = [0,1,2]
title = 'e1 Measurement for exponential Galaxies  '+ 'aratio='+print_ratio
;title = 'e1 Measurement for deVaucouleurs Galaxies  '+ 'aratio='+print_ratio
plot,s2n,e1_s2n,psym=1,xtitle='Flux S/N',ytitle='e1 S/N',title=title, $
	yrange=[0.0, max(ad_s2n)]
oplot,s2n,ad_s2n,psym=4
legend,['Adaptive', w_type],psym=[4,1]
devid = ad_s2n/e1_s2n
plot,s2n,devid, psym=4, title = 'Ratio of two e1 S/N curves', $
	xtitle='Flux S/N', ytitle='e1 S/N ratio',yrange=[0,max(devid)]

!p.multi = [0,1,1]
return
end


