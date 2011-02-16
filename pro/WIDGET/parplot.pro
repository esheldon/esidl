PRO parplot_event, ev

  COMMON par, command, npar, v
  widget_control, ev.id, get_uvalue=type

  IF type EQ npar THEN BEGIN
      widget_control,/destroy, ev.top 
  ENDIF ELSE BEGIN
    
    commstr = command
    FOR i=0, npar-1 DO BEGIN
      v[i] = (i EQ type)*(ev.value - v[i]) + v[i]
      commstr = commstr + ','+ntostr(v[i])
    ENDFOR
    result = execute(commstr)

  ENDELSE 

END

PRO parplot

  COMMON par, command, npar, v

  command = ' '
  print,format='($, "Command Name")'
  read, command
  print,format='($, "Number of Parameters")'
  read, npar
  npar=fix(npar)

  v=fltarr(npar)

  device, get_screen_size=s

  xsize=.5*s[0]
  ysize=.5*s[1]

  base = widget_base(/column, title='Parplot', xoffset=0,yoffset=0)
  
  draw = widget_draw(base, xsize=xsize, ysize=ysize)
  b = widget_button(base, value='DONE', uvalue=npar)

  widget_control, base, /realize, update=0
  widget_control, draw, get_value=win_id
  wset, win_id

  FOR i=0, npar-1 DO slide = cw_fslider(base, /drag, value=0, $
           minimum=0, maximum=1,title = 'Parameter '+ntostr(i+1), uvalue=i)
  
  widget_control, base, /map, /update

  xmanager, 'PARPLOT', base

END
