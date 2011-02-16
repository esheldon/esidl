PRO testphotomap, str

  order=2
  photomap, str, order, str.ra, str.dec, outrow, outcol

  diffrow=str.rowc[2]-outrow
  diffcol=str.colc[2]-outcol
  diff = sqrt(diffrow^2 + diffcol^2)
  w=where(diff LE 1.)
  plothist, diff[w], bin=.05

  key=get_kbrd(1)

  order=3
  photomap, str, order, str.ra, str.dec, outrow, outcol

  diffrow=str.rowc[2]-outrow
  diffcol=str.colc[2]-outcol
  diff = sqrt(diffrow^2 + diffcol^2)
  w=where(diff LE 1.)
  plothist, diff[w], bin=.05

  key=get_kbrd(1)

  order=4
  photomap, str, order, str.ra, str.dec, outrow, outcol

  diffrow=str.rowc[2]-outrow
  diffcol=str.colc[2]-outcol
  diff = sqrt(diffrow^2 + diffcol^2)
  w=where(diff LE 1.)
  plothist, diff[w], bin=.05

  key=get_kbrd(1)

  order=5
  photomap, str, order, str.ra, str.dec, outrow, outcol

  diffrow=str.rowc[2]-outrow
  diffcol=str.colc[2]-outcol
  diff = sqrt(diffrow^2 + diffcol^2)
  w=where(diff LE 1.)
  plothist, diff[w], bin=.05


return
END 
