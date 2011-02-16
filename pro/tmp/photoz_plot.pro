PRO photoz_plot, str, photoz



  w = where( str.petrocounts[2] LT 22. AND str.petrocounts[2] GT 18)
  help,w
  ww = where(str[w].type[2] EQ 3 AND str[w].type[1] EQ 3)
  
  w = w[ww]
  help,w
  ww = where(str[w].petrocountserr[0] LT 1. AND $
             str[w].petrocountserr[1] LT 1. AND $
             str[w].petrocountserr[2] LT 1. AND $
             str[w].petrocountserr[3] LT 1. )
  w = w[ww]
  help,w
  photoz = photozpoly(str[w].petrocounts[0], $
                      str[w].petrocounts[1], $
                      str[w].petrocounts[2], $
                      str[w].petrocounts[3] ) > 0. < 10

  xrange = [.01,1.8]
  plothist, photoz, bin=.1,xrange=xrange

  return
END 
