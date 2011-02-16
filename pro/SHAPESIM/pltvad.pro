pro pltvad,tot,actual,w_type,aratio,delta

if n_params() eq 0 then begin
	print,'-syntax  pltvad,tot,actual,w_type,aratio,delta  aratio is a string'
	return
endif


w = where(tot[*].s2n le 400 and tot[*].e1_ad ge -2)
tot2 = tot[w]
n = ((n_elements(tot2))/10)*10

s= sort(tot2[*].s2n)

s2n = tot2[s].s2n
e1 = tot2[s].e1
e1_lupton = tot2[s].e1_lupton*2
e1_ad = tot2[s].e1_ad
e1_gaussian = tot2[s].e1_gaussian
e1_cowboy = tot2[s].e1_cowboy

s2n = s2n[1:n]
e1 = e1[1:n]
e1_lupton = e1_lupton[1:n]
e1_ad = e1_ad[1:n]
e1_gaussian = e1_gaussian[1:n]
e1_cowboy = e1_cowboy[1:n]

s2n = rebin(s2n,n/10)
e1 = rebin(e1,n/10)
e1_lupton = rebin(e1_lupton,n/10)
e1_ad = rebin(e1_ad,n/10)
e1_gaussian = rebin(e1_gaussian,n/10)
e1_cowboy = rebin(e1_cowboy,n/10)

e1 = smooth(e1,5)
e1_lupton = smooth(e1_lupton,5)
e1_ad = smooth(e1_ad,5)
e1_gaussian = smooth(e1_gaussian,5)
e1_cowboy = smooth(e1_cowboy,5)

ymin=actual-delta
ymax=actual+delta
if w_type eq 'none' then begin
	tit = 'adaptive and unweighted  aratio = '+aratio
	plot,s2n,e1,psym=1,yrange=[ymin,ymax],xtitle='s2n',ytitle='e1',$
	title = tit
endif 
if w_type eq 'lupton' then begin
	tit = 'adaptive and lupton  aratio = '+aratio
	plot,s2n,e1_lupton,psym=1,yrange=[ymin,ymax],xtitle='s2n',ytitle='e1',$
	title = tit
endif
if w_type eq 'gaussian' then begin
	tit = 'adaptive and gaussian aratio = '+aratio
	plot,s2n,e1_gaussian,psym=1,yrange=[ymin,ymax],xtitle='s2n',ytitle='e1',$
	title = tit
endif
if w_type eq 'cowboy' then begin
	tit = 'adaptive and cowboy aratio = '+aratio
	plot,s2n,e1_cowboy,psym=1,yrange=[ymin,ymax],xtitle='s2n',ytitle='e1',$
	title = tit
endif

oplot,s2n,e1_ad,psym=4
oplot,[0,400],[actual,actual]

return
end




