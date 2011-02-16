function xi_int_kernal,x
return,(1d +x^2)^(!xi_int_kernal_exp)
end

function xi_pl_int,R0,R1,r,b
;perform the numerical integral 
; \frac{1}{\pi} \int_R0^R1 \frac{dR -W'(R)}{\sqrt{R^2-r_j^2} 
;when -W'(R) is a powerlaw R^(-b)
;transforms variables first

x1=sqrt((R1/r)^2-1)
x0=sqrt((R0/r)^2-1)

defsysv,"!xi_int_kernal_exp",-(b+1.0)/2.0
func="xi_int_kernal"
return,qromb(func,x0,x1,/double,eps=1e-8)*r^(-b)/!pi
end
