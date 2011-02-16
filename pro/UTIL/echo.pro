
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    ECHO
;       
; PURPOSE:
;    call system function echo, allow user to format the text the same
;    way as with echo (colors, underline, bold, etc)
;
; CALLING SEQUENCE:
;    echo, string, color=color, bcolor=bcolor, bold=bold, $
;          underscore=underscore, $
;          reverse=reverse, nonewline=nonewline
; COMMENTS:
;    This procedure uses the spawn function. If your on linux it will
;    be pretty quick, because I know echo takes the -e option, so I 
;    can use /noshell in spawn. Othersize I need to create a shell which
;    is slow. 
;
; INPUTS: 
;    string:  a string (or array of string) to ouput with echo
;
; OPTIONAL INPUTS:
;    color: One of the standard ansii colors (or an array)
;        'black','red','green','yellow','blue','magenta','cyan','white'
;        or 'none' for no color formatting
;    bclolor: same as color but the background for the text
;    bold: 1 or 0 for bold or not bold
;    underscore: " "
;    reverse: " "
;    nonewline: 1 or 0 for putting a newline character after string or not
;
; KEYWORD PARAMETERS:
;    bold,underscore,reverse,nonnewline above are like keywords, but they can
;    be arrays
;       
; OUTPUTS: 
;    echo to stdout
;
; OPTIONAL OUTPUTS:
;    NONE
;
;
; CALLED ROUTINES:
;    ECHO_CHECKPARAM (subroutine)
; 
;
; REVISION HISTORY:
;    Created 1-Mar-2002 Erin Scott Sheldon (UofM)
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


FUNCTION echo_checkparam, paramarr, np, nst, stind

  IF np NE nst THEN BEGIN 
      IF np EQ 1 THEN return,paramarr[0]
      message,'# elements in each keyword must be 1 or same as # in string array'
  ENDIF ELSE BEGIN 
      return,paramarr[stind]
  ENDELSE 

END 

PRO echo, stringarr, colorarr=colorarr, bcolorarr=bcolorarr, boldarr=boldarr, $
          underscorearr=underscorearr, $
          reversearr=reversearr, nonewlinearr=nonewlinearr

  IF n_params() LT 1 THEN BEGIN
      print,'-Syntax: echo, string, color=color, bcolor=bcolor, /bold, $'
      print,'         /underscore, /reverse, /nonewline'
      return
  ENDIF 

  ;; # of elements in string array
  nst=n_elements(stringarr)
  nc=n_elements(colorarr)
  IF nc EQ 0 THEN colorarr = 'none' 
  nc=nc>1

  nbc=n_elements(bcolorarr)
  IF nbc EQ 0 THEN bcolorarr = 'none'
  nbc=nbc>1

  nb=n_elements(boldarr)
  IF nb EQ 0 THEN boldarr=0
  nb=nb>1

  nu=n_elements(underscorearr)
  IF nu EQ 0 THEN underscorearr=0
  nu=nu>1

  nr=n_elements(reversearr)
  IF nr EQ 0 THEN reversearr=0
  nr=nr>1

  nnl=n_elements(nonewlinearr)
  IF nnl EQ 0 THEN nonewlinearr=0
  nnl=nnl>1

  ;; loop over strings
  FOR i=0L, nst-1 DO BEGIN 

      string = stringarr[i]

      color=echo_checkparam(colorarr, nc, nst, i)
      bcolor=echo_checkparam(bcolorarr, nbc, nst, i)
      bold=echo_checkparam(boldarr, nb, nst, i)
      underscore=echo_checkparam(underscorearr, nu, nst, i)
      reverse=echo_checkparam(reversearr, nr, nst, i)
      nonewline=echo_checkparam(nonewlinearr, nnl, nst, i)

      CASE strlowcase(color) OF
          'black': cstring='30'
          'red':   cstring='31'
          'green': cstring='32'
          'yellow': cstring='33'
          'blue': cstring='34'
          'magenta': cstring='35'
          'cyan': cstring='36'
          'white': cstring='37'
          'none': cstring=''
          ELSE: message,'Unknown color string: "'+color+'"'
      ENDCASE 

      CASE strlowcase(bcolor) OF
          'black': cstring=cstring+';40'
          'red':   cstring=cstring+';41'
          'green': cstring=cstring+';42'
          'yellow': cstring=cstring+';43'
          'blue': cstring=cstring+';44'
          'magenta': cstring=cstring+';45'
          'cyan': cstring=cstring+';46'
          'white': cstring=cstring+';47'
          'none': 
          ELSE: message,'Unknown background color string: "'+bcolor+'"'
      ENDCASE 

      IF (!version.os EQ 'Linux') OR (!version.os EQ 'OSF1') THEN BEGIN 
          ;; we know how to deal with escapes for /noshell in spawn (much
          ;; faster)
          command='echo'
          IF !version.os EQ 'Linux' THEN command = ['echo','-e']
          IF keyword_set(nonewline) THEN command = [command, '-n']
          IF keyword_set(bold) THEN BEGIN 
              cstr = '\033['+cstring+';01m'+string+'\033[0;00m'
          ENDIF ELSE IF keyword_set(underscore)THEN BEGIN 
              cstr='\033['+cstring+';04m'+string+'\033[0;00m'
          ENDIF ELSE IF keyword_set(reverse) THEN BEGIN 
              cstr='\033['+cstring+';07m'+string+'\033[0;00m'
          ENDIF ELSE BEGIN 
              cstr='\033['+cstring+'m'+string+'\033[0;00m'
          ENDELSE 
          command = [command,cstr]
          spawn,command,/noshell
      ENDIF ELSE BEGIN  
    
          ;; can't necessarily use -e, so use slower spawn with shell
          command = 'echo '
          IF keyword_set(nonewline) THEN command = command + ' -n '
          IF keyword_set(bold) THEN BEGIN 
              spawn,command+'"\033['+cstring+';01m'+string+'\033[0;00m"'
          ENDIF ELSE IF keyword_set(underscore)THEN BEGIN 
              spawn,command+'"\033['+cstring+';04m'+string+'\033[0;00m"'
          ENDIF ELSE IF keyword_set(reverse) THEN BEGIN 
              spawn,command+'"\033['+cstring+';07m'+string+'\033[0;00m"'
          ENDIF ELSE BEGIN 
              spawn,command+'"\033['+cstring+'m'+string+'\033[0;00m"'
          ENDELSE 
      ENDELSE 

  ENDFOR 

END 
