PRO  widget3_event, ev

  COMMON wid3, seed
  
  widget_control, ev.top, get_uvalue=drawid
  wset, drawid

  IF ( tag_names(ev, /structure_name) EQ 'WIDGET_TIMER')  THEN BEGIN
    table = fix(randomu(seed)*41)
    loadct, table

    widget_control, ev.id, timer=3.0
  ENDIF 

  IF ( tag_names(ev, /structure_name) EQ 'WIDGET_DROPLIST') THEN BEGIN
    
    CASE ev.index OF 
      0: plot, dist(500)
      1: surface, dist(500)
      2: shade_surf, dist(500)
      3: widget_control, ev.top, /destroy
    ENDCASE 

  ENDIF 

END

PRO widget3

  select = ['Plot', 'Surface', 'Shaded Surface', 'Done']
  base = widget_base(/column)
  draw=widget_draw(base, xsize=500, ysize=500)
  dlist = widget_droplist(base, value=select)

  widget_control, base, /realize
  widget_control, draw, get_value=drawid
  widget_control, base, set_uvalue=drawid
  widget_control, draw, timer=0.0

  widget_control, dlist, set_droplist_select=2
  wset, drawid

  shade_surf, dist(500)
  xmanager, 'widget3', base

END
