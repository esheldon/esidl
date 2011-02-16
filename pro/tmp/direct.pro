pro direct, u, kappa, lx=lx, ly=ly

; Performs mass inversion using the direct method (Marco Lombardi, 1999).
;
; u:	  [NX, NY, 2] array with the u vector field (input)
; kappa:  [NX, NY] array containing the reconstruced mass map (output)
; lx, ly: sizes of the rectangle in some units (e.g., arcmin)

usize = size(u)
if usize[0] ne 3 then $
  message, "u parameter must have 3 indices"
if usize[3] ne 2 then $
  message, "u parameter does not represent a vector field"
nx = usize[1]
ny = usize[2]

nx2 = 2*nx - 1
ny2 = 2*ny - 1

if n_elements(lx) ne 1 then lx = nx
if n_elements(ly) ne 1 then ly = ny

kappa = dblarr(nx, ny)
u1 = dblarr(nx2, ny2)
u2 = u1
c = u1


u1[0:nx-1, 0:ny-1] = u[*, *, 0]
u1[0:nx-1, ny:*] = reverse(u[*, 1:*, 0], 2)
u1[nx:*, *] = -reverse(u1[1:nx-1, *], 1)

u2[0:nx-1, 0:ny-1] = u[*, *, 1]
u2[nx:*, 0:ny-1] = reverse(u[1:*, *, 1], 1)
u2[*, ny:*] = -reverse(u2[*, 1:ny-1], 2)

out1 = -imaginary(fft(u1, 1, /double))
out2 = -imaginary(fft(u2, 1, /double))

x = dindgen(nx) # replicate(1.0, ny)
y = replicate(1.0, nx) # dindgen(ny)
x[0, 0] = 1.0
y[0, 0] = 1.0
coeff1 = x / (((x/lx)^2 + (y/ly)^2) * $
  4.0*!DPI*nx2*ny2*lx)
coeff2 = y / (((x/lx)^2 + (y/ly)^2) * $
  4.0*!DPI*nx2*ny2*ly)

c[1:nx-1, 0:ny-1] =  (out1[1:nx-1, 0:ny-1] - $
  reverse(out1[nx:nx2-1, 0:ny-1], 1)) * $
  coeff1[1:*, *] + (out2[1:nx-1, 0:ny-1] + $
  reverse(out2[nx:nx2-1, 0:ny-1], 1)) * $
  coeff2[1:*, *]
c[0, 0:ny-1] = 2.0 * out2[0, 0:ny-1] * coeff2[0, *]
c[0, 0] = 0.0

c[nx:nx2-1, 0:ny-1] = reverse(c[1:nx-1, 0:ny-1], 1)
c[*, ny:ny2-1] = reverse(c[*, 1:ny-1], 2)

out1 = double(fft(c, 1, /double))

kappa[1:*, *] = out1[1:nx-1, 0:ny-1] + $
  reverse(out1[nx:*, 0:ny-1], 1)
kappa[0, *] = 2.0 * out1[0, 0:nx-1]

end



  
