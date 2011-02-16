FUNCTION shapelet_decomp_png_name, str, clr, type, nmax, idString=idString

  c = ['u','g','r','i','z']
  idString = shapelet_idstring(str)

  nmaxstr = strn(nmax, len=2, padchar='0')
  name = $
    'shapeletDecomp-'+c[clr]+'-'+type+'-'+'nmax'+nmaxstr+'-' + $
    idString + '.png'

  return,name
END 
