PRO widget2_event, ev
  
  widget_control, ev.top, get_uvalue=textwid
  widget_control, ev.id, get_uvalue=uval

  CASE uval OF 
    'ONE' : widget_control, textwid, set_value='Button One Pressed'
    'TWO' : widget_control, textwid, set_value='Button Two Pressed'
    'DONE' : widget_control, ev.top, /destroy
  ENDCASE 

END 

PRO widget2
  
  base = widget_base(/column)
  button1 = widget_button(base, value='One', uvalue='ONE')
  button2 = widget_button(base, value='Two', uvalue='TWO')

  text = widget_text(base, xsize=20)

  button3 = widget_button(base, value='Done', uvalue='DONE')

  widget_control, base, set_uvalue = text
  widget_control, base, /realize

  xmanager, 'Widget2', base

END
