pro plt_all,tot,actual,aratio,delta

if n_params() eq 0 then begin
	print,'-syntax  plt_all,tot,actual,aratio,delta  aratio is a string'
	return
endif

path = '/sdss/data1/esheldon/PS/NEWWAY/fwhm10/'

name = path + 'ad_lup' + aratio + '.ps'
type = 'lupton'
begplot,name=name,/invbw,/landscape
pltvad,tot,actual,type,aratio,delta
endplot,/noprint

name = path + 'ad_gaussian' + aratio + '.ps'
type = 'gaussian'
begplot,name=name,/invbw,/landscape
pltvad,tot,actual,type,aratio,delta
endplot,/noprint

name = path + 'ad_unweighted' + aratio + '.ps'
type = 'none'
begplot,name=name,/invbw,/landscape
pltvad,tot,actual,type,aratio,delta
endplot,/noprint

name = path + 'ad_cowboy' + aratio + '.ps'
type = 'cowboy'
begplot,name=name,/invbw,/landscape
pltvad,tot,actual,type,aratio,delta
endplot,/noprint



return
end






