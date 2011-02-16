pro plot_e1_s2n,read

if n_params() eq 0 then begin
	print,'-syntax plot_e1_s2n,read'
	return
endif

print_dir = "/sdss4/data1/esheldon/PS/Exp/S2N/"
read_dir = "/sdss4/data1/esheldon/"
xmin = 0.0
xmax = 100.0


ymin=-.5
ymax=.5


if read eq 1 then begin
	fname = read_dir + 'mids2nexptot1.fits'
	tot1 = mrdfits(fname,1,hdr)
endif


;get_stars,tot1,w
find_e1,1,e1

begplot,name=print_dir+"admomVs2n_ar1.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_ad,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="adaptive weighting, aratio=1"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"luptonVs2n_ar1.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_lupton,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="lupton weighting, aratio=1"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"noweightVs2n_ar1.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="no weighting, aratio=1"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"cowboyVs2n_ar1.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_cowboy,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="cowboy hat weighting, aratio=1"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"gaussianVs2n_ar1.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_gaussian,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="gaussian weighting, aratio=1"
oplot,[0,400],[e1,e1]
endplot,/noprint

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ymin=-.6
ymax=1.6


if read eq 1 then begin
	fname = read_dir + 'mids2nexptot9.fits'
	tot1 = mrdfits(fname,1,hdr)
endif


;get_stars,tot1,w
find_e1,.9,e1

begplot,name=print_dir+"admomVs2n_ar9.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_ad,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="adaptive weighting, aratio=.9"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"luptonVs2n_ar9.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_lupton,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="lupton weighting, aratio=.9"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"noweightVs2n_ar9.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="no weighting, aratio=.9"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"cowboyVs2n_ar9.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_cowboy,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="cowboy hat weighting, aratio=.9"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"gaussianVs2n_ar9.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_gaussian,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="gaussian weighting, aratio=.9"
oplot,[0,400],[e1,e1]
endplot,/noprint

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ymin=-.4
ymax=.6

if read eq 1 then begin
	fname = read_dir + 'mids2nexptot8.fits'
	tot1 = mrdfits(fname,1,hdr)
endif
;get_stars,tot1,w
find_e1,.8,e1

begplot,name=print_dir+"luptonVs2n_ar8.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_lupton,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="lupton weighting, aratio=.8"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"admomVs2n_ar8.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_ad,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="adaptive weighting, aratio=.8"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"noweightVs2n_ar8.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="no weighting, aratio=.8"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"cowboyVs2n_ar8.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_cowboy,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="cowboy hat weighting, aratio=.8"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"gaussianVs2n_ar8.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_gaussian,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="gaussian weighting, aratio=.8"
oplot,[0,400],[e1,e1]
endplot,/noprint

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ymin=0.0
ymax=.6

if read eq 1 then begin
	fname = read_dir + 'mids2nexptot6.fits'
	tot1 = mrdfits(fname,1,hdr)
endif
;get_stars,tot1,w
find_e1,.6,e1

begplot,name=print_dir+"luptonVs2n_ar6.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_lupton,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="lupton weighting, aratio=.6"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"admomVs2n_ar6.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_ad,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="adaptive weighting, aratio=.6"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"noweightVs2n_ar6.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="no weighting, aratio=.6"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"cowboyVs2n_ar6.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_cowboy,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="cowboy hat weighting, aratio=.6"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"gaussianVs2n_ar6.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_gaussian,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="gaussian weighting, aratio=.6"
oplot,[0,400],[e1,e1]
endplot,/noprint

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ymin=-1.0
ymax=1.0

if read eq 1 then begin
	fname = read_dir + 'mids2nexptot95.fits'
	tot1 = mrdfits(fname,1,hdr)
endif
;get_stars,tot1,w
find_e1,.95,e1

begplot,name=print_dir+"luptonVs2n_ar95.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_lupton,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="lupton weighting, aratio=.95"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"admomVs2n_ar95.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_ad,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="adaptive weighting, aratio=.95"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"noweightVs2n_ar95.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="no weighting, aratio=.95"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"cowboyVs2n_ar95.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_cowboy,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="cowboy hat weighting, aratio=.95"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"gaussianVs2n_ar95.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_gaussian,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="gaussian weighting, aratio=.95"
oplot,[0,400],[e1,e1]
endplot,/noprint

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

xmin = 0.0
xmax = 140.0

ymin=.2
ymax=.8

if read eq 1 then begin
	fname = read_dir + 'mids2nexptot4.fits'
	tot1 = mrdfits(fname,1,hdr)
endif
;get_stars,tot1,w
find_e1,.4,e1

begplot,name=print_dir+"luptonVs2n_ar4.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_lupton,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="lupton weighting, aratio=.4"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"admomVs2n_ar4.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_ad,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="adaptive weighting, aratio=.4"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"noweightVs2n_ar4.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="no weighting, aratio=.4"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"cowboyVs2n_ar4.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_cowboy,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="cowboy hat weighting, aratio=.4"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"gaussianVs2n_ar4.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_gaussian,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="gaussian weighting, aratio=.4"
oplot,[0,400],[e1,e1]
endplot,/noprint

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ymin=.5
ymax=1.0

if read eq 1 then begin
	fname = read_dir + 'mids2nexptot3.fits'
	tot1 = mrdfits(fname,1,hdr)
endif
;get_stars,tot1,w
find_e1,.3,e1

begplot,name=print_dir+"luptonVs2n_ar3.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_lupton,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="lupton weighting, aratio=.3"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"admomVs2n_ar3.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_ad,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="adaptive weighting, aratio=.3"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"noweightVs2n_ar3.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="no weighting, aratio=.3"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"cowboyVs2n_ar3.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_cowboy,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="cowboy hat weighting, aratio=.3"
oplot,[0,400],[e1,e1]
endplot,/noprint

begplot,name=print_dir+"gaussianVs2n_ar3.ps",/invbw,/landscape
plot,tot1.s2n,tot1.e1_gaussian,psym=3,yrange=[ymin,ymax],xrange=[xmin,xmax],ytitle="e1",xtitle = "s2n",title="gaussian weighting, aratio=.3"
oplot,[0,400],[e1,e1]
endplot,/noprint

return
END







