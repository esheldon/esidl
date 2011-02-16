
;+
; NAME:
;    PSLAYOUT
;
;
; PURPOSE:
;    Set up the system variable !PSLAYOUT which contains information for
;    setting up the X, Z, or PS devices.  The user, if not satisfied with these
;    settings, should copy pslayout into their code directory and make the
;    required modifications.
;
; CATEGORY:
;    Plotting routine

; CALLING SEQUENCE:
;    pslayout
;
; SIDE EFFECTS:
;    The !PSLAYOUT variable is created if not already in existence.  If in
;    existence, it will NOT be modified.
;
;
; MODIFICATION HISTORY:
;    Creation:  ??-??-2001 Erin Sheldon UofMichigan
;-
;
;
;
;  Copyright (C) 2005  Erin Sheldon, NYU.  erin dot sheldon at gmail dot com
;
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program; if not, write to the Free Software
;    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
;
;

PRO pslayout

  defsysv,'!pslayout',exist=exist

  IF NOT exist THEN BEGIN 
      
      pslayout = {runsetup:1,$  ;should begplot/endplot run setupplot?
                  name:'ps', $
                  bits_per_pixel:4, $
                  bold:0, $
                  book:0, $
                  close_file:0, $
                  color:0, $
                  demi:0, $
                  Encapsulated:0, $
                  filename:'idl.ps', $
                  font:'times', $ ;the default true-type or postscript font
                  isolatin1:1, $
                  true: 0, $    ;use true-type for all fonts? Obsolete, see below
                  font_index:6, $ 
                  font_size:12, $
                  italic:0, $
                  inches:1, $
                  landscape:0, $
                  light:0, $
                  medium:0, $
                  narrow:0, $
                  oblique:0, $
                  output:'', $
                  portrait:1, $
                  preview:1, $
                  scale_factor:1.0, $
                  xoffset:0.75, $
                  xsize:7.0, $
                  yoffset:1.0, $
                  ysize:9.0,$
                  invbw:1,$
                  $             ; default plotting parameters for 'ps' device
                  ps_true:0, $
                  ps_thick:5, $
                  ps_xthick:5, $
                  ps_ythick:5, $
                  ps_xticklen:0.02,$
                  ps_yticklen:0.02,$
                  ps_charsize:1.5,$
                  ps_charthick:4,$ ;charthick only has meaning for true=0 above
                  $             ; default plotting parameters for 'x' device
                  x_true:0,$
                  x_thick:1, $
                  x_xthick:1, $
                  x_ythick:1, $
                  x_xticklen:0.02,$
                  x_yticklen:0.02,$
                  x_charsize:1.5,$
                  x_charthick:1, $
                  $
                  z_true:0,$
                  z_thick:1, $
                  z_xthick:1, $
                  z_ythick:1, $
                  z_xticklen:0.02,$
                  z_yticklen:0.02,$
                  z_charsize:0.8,$ ;smaller charsize
                  z_charthick:1,$
                  z_resolution:[640,512], $;resolution of Z-buffer
                  load_simpctable: 0} ; run simpctable?


      defsysv,'!pslayout',pslayout

  ENDIF
  

END 
