pro plt_all_moms2n,tot,aratio

if n_params() eq 0 then begin
	print,'-syntax  plt_all_moms2n,tot,aratio    aratio is a string used'
	print,'in the filename and plot title.  dont give the decimal point'
	return
endif


print,'---------------------------------------------'

path = '/sdss4/data1/esheldon/PS/Exp/MOMS2N/'
tail = aratio+'.ps'

name = path + 'moms2n_lup'+tail
print,name
type = 'lupton'
begplot,name=name,/invbw,/landscape
plt_moms2n,tot,type,aratio
endplot,/noprint

name = path + 'moms2n_unweighted'+tail
print,name
type = 'none'
begplot,name=name,/invbw,/landscape
plt_moms2n,tot,type,aratio
endplot,/noprint

name = path + 'moms2n_gaussian'+tail
print,name
type = 'gaussian'
begplot,name=name,/invbw,/landscape
plt_moms2n,tot,type,aratio
endplot,/noprint

name = path + 'moms2n_cowboy'+tail
print,name
type = 'cowboy'
begplot,name=name,/invbw,/landscape
plt_moms2n,tot,type,aratio
endplot,/noprint

return
end