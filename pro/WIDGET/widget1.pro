PRO widget1_event, ev
  IF ev.select THEN widget_control, ev.top, /destroy
END 

PRO widget1
  base = widget_base()
  button = widget_button(base, value='Push this button to kill me')
  widget_control, base, /realize
  xmanager, 'Widget1', base
END
