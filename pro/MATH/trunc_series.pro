pro trunc_series, theta, order, terms, cumul, total

if n_params() eq 0 then begin
	print,"-syntax trunc_series, theta, order, terms, cumul, total"
	return
endif

d = .01

cumul=findgen(order)
k=findgen(order)
my_fact, 2*k, k2fact
my_fact, k, kfact
numerator = 1.0*k2fact*k*(sin(theta))^(2*k+1)
denominator = (1.0*kfact*(2*k+1)*2^(2*k))

terms = numerator/denominator
for i=0, order-1 do begin
	cumul[i] = total(terms[0:i])
endfor
total=total(terms)


return 
end