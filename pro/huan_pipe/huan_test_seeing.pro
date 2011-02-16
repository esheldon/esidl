PRO huan_test_seeing, im, cat, wstar

  arcperpix = 0.3

  ;; Get some bright ones
  w=where(cat[wstar].Iiso lt 21, nw)

  w=wstar[w]


  ;; coord in sub-image
  cen = [10,10]

  nx = 21
  ny = 21

  index=lindgen(nx*ny)

  fx = findgen(nx) - 10.
  fy = findgen(ny) - 10.

  x=index MOD nx
  y=index/nx

  xp = x - cen[0]
  yp = y - cen[1]

  rr = xp^2 + yp^2
  r = sqrt(rr)

  s = sort(r)

  !p.multi=[0,0,2]
  FOR i=0L, nw-1 DO BEGIN 

      x = cat[w[i]].x_image - 1.0
      y = cat[w[i]].y_image - 1.0


      ;; extract a sub-image
      tim = im[ x-10:x+10, y-10:y+10 ]

      tvasinh, tim



      res = gauss2dfit(tim, AA, fx, fy,/tilt)

      maxs = max([AA[2],AA[3]], min=mins)
      seeing_admom = sqrt( (cat[w[i]].ixx+cat[w[i]].iyy)/2.)*2.35*0.3
      print
      print,'I_AB: ',cat[w[i]].Iiso
      print,'seeing x:   ',AA[2]*2.35*arcperpix
      print,'seeing y:   ',AA[3]*2.35*arcperpix
      print,'seeing admom',seeing_admom

      ;model = tim
      model = AA[0] + AA[1]*exp(-0.5*( (xp/AA[2])^2 + (yp/AA[3])^2 ) )

      plot, r[s], tim[s], psym=8, $
        xrange=[0, 10]
      ;; oplot, r[s], tim[s]
      oplot, r[s], model[s], psym=8, color=!red
      oplot, r[s], model[s], color=!red


      key=prompt_kbrd()
      IF key EQ 'q'THEN return


  ENDFOR 


END 
