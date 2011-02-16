function hermite, n, x
  
; July 99 - Written by A. Refregier
;
; PURPOSE: compute the hermite polynomial Hn(x) of order n.
; INPUT: n,x: output Hn(x) where x can be a scalar, a vector or an array
; OUTPUT: Hn(x)
  
  nx = n_elements(x)
  
  case n of
     0: begin
          sx = size(x)
          case sx(0) of
             0: h = 1.
             1: h = replicate(1., sx(1))
             2: h = replicate(1., sx(1), sx(2))
          else: begin
                  print, 'hermite: x must be a scalar, a vector or an array'
                end
          endcase
        end
     1: h = 2.*x
     2: h = 4.*x^2-2.
     3: h = 8*x^3-12.*x
     4: h = 16*x^4-48.*x^2+12.
     5: h = 32*x^5-160.*x^3+120.*x
     6: h = 64.*x^6-480.*x^4+720.*x^2-120.
  else: begin
          c = hermitecof(n)     
          sx = size(x)
          if sx(0) eq 0 then h = total(c*x^findgen(n+1)) else begin
             h = c(0)+c(1)*x
             for i=2, n do h = h+c(i)*x^i
          endelse
        end
  endcase
     
  
  return, h
  end  
