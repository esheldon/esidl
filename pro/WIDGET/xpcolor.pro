PRO xpcolor_event, event

  COMMON colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr

  widget_control, get_uvalue=type, event.id

  IF type EQ 0 THEN BEGIN
    IF event.value EQ 0 THEN widget_control, /destroy, event.top $
    ELSE xdisplayfile, title='XPcolor Help', $
      group=event.top, height=5, width=60, text='blah blah'
  ENDIF ELSE BEGIN
    r_curr(!p.color) = event.r
    g_curr(!p.color) = event.g
    b_curr(!p.color) = event.b
    tvlct, event.r, event.g, event.b, !p.color
  ENDELSE 

END 

PRO xpcolor

  COMMON colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr

  IF n_elements(r_orig) LE 0 THEN BEGIN
    tvlct, r_orig, g_orig, b_orig, /get
    r_curr=r_orig & b_curr=b_orig & g_curr=g_orig
  ENDIF 

  base = widget_base(/column, title='Set Plot Color')
  
  b = cw_bgroup(base, /row, ['Done', 'Help'], uvalue=0)

  draw = widget_draw(base, xsize=100, ysize=50)
  
  rgb = cw_rgbslider(base, /drag, uvalue=1)

  widget_control, rgb, set_value=[r_curr(!p.color), g_curr(!p.color), $
                                  b_curr(!p.color)]

  widget_control, /realize, base

  widget_control, draw, get_value=win_id

  save=!d.window

  wset,win_id
  erase, color=!p.color
  
  IF save NE -1 THEN wset, save

  xmanager, 'XPCOLOR', base

END 
