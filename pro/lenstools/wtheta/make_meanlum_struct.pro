PRO make_meanlum_struct, infile

  ;; infile should be one of those matchgN5*lensum*, or matchrN5, etc 

  outfile = repstr(infile, 'lensum', 'meanlum')

  print,'input file: ',infile
  print,'output file: ',outfile

  t=mrdfits(infile,1)

  meanlum=fltarr(5)
  meanlumerr=meanlum

  s=create_struct('meanlum', meanlum, $
                  'meanlumerr', meanlumerr)

  units = 1.e10

  FOR i=0L, 4 DO BEGIN 

      w=where( (t.lum[i] GT 0.0) AND (t.lum[i]/1.e10 LT 20.0), nw)

      lensave, t[w], 'lum', mean, err, element=i
      meanlum[i] = mean/units
      meanlumerr[i] = err/units

  ENDFOR 

  forprint,meanlum,meanlumerr

  s=create_struct('meanlum', meanlum, $
                  'meanlumerr', meanlumerr)

  print,'outputting file: ',outfile
  mwrfits, s, outfile, /create

END 
