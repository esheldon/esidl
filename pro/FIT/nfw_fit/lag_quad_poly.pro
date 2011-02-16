;+
; NAME: lag_quad_poly
;
; PURPOSE:
; Computes the uniq quadratic poly through three x,y points
;
; CALLING SEQUENCE:
;    p=lag_quad_poly(x,y)
;
; INPUTS:
; x,y -  set of three points
;
; OUTPUTS:
;     polynomial -  three numbers , usual IDL notation
;
; METHOD:
;     Uses Lagrange formula
;-

function lag_quad_poly,x,y

n=n_elements(x)

if n ne 3 then begin
    print,'ERROR - x,y must have exactly 3 points each' 
    return ,-1
endif

q=dblarr(3)
p=q

q[0]=y[0]/((x[0]-x[1])*(x[0]-x[2]))
q[1]=y[1]/((x[1]-x[0])*(x[1]-x[2]))
q[2]=y[2]/((x[2]-x[0])*(x[2]-x[1]))

p[0]=q[0]*x[1]*x[2] + q[1]*x[0]*x[2] + q[2]*x[0]*x[1]
p[1]=-q[0]*(x[1]+x[2])-q[1]*(x[0]+x[2]) -q[2]*(x[0]+x[1])
p[2]=total(q)

return,p
end
