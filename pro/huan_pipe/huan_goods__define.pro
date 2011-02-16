;; Dealing with the seeing values, files, and plotting for the goods
;; noise-added images. 
;;
;; Currently only has the convolution stuff

FUNCTION huan_goods::init, type

  IF n_elements(type) EQ 0 THEN BEGIN 
      message,'You must initialize the type',/inf
      message,'-Syntax: hg=obj_new("huan_goods",type)',/inf
      return,0
  ENDIF 

  self->set_parameters, type=type, status=status
  IF status NE 0 THEN return,0

  ;; Read the seeing
  self->read_seeing, status=status
  IF status NE 0 THEN return,0

  return,1

END 

PRO huan_goods::set_parameters, type=type, status=status

  status = 1
  IF n_elements(type) NE 0 THEN BEGIN 

      IF size(type, /tname) NE 'STRING' THEN BEGIN 
          message,'You must enter a string for the type',/inf
          return
      ENDIF 
      self.type = strlowcase( type[0] )
  ENDIF 

  status = 0

END 

FUNCTION huan_goods::type
  return,self.type
END 


FUNCTION huan_goods::dir

  dir = esheldon_config('des_sensitivity_dir')
  return,dir

END 

FUNCTION huan_goods::sensitivity_file, type=type
  
  IF n_elements(type) EQ 0 THEN type=self.type
  file = 'sensitivity/'+type+'_convolved.st'
  file = concat_dir(self->dir(), file)
  return,file
END 
FUNCTION huan_goods::read_sensitivity, type=type, status=status

  status=1
  file = self->sensitivity_file(type=type)
  IF NOT fexist(file) THEN BEGIN 
      message,'Sensitivity file does not exist: '+file,/inf
      return,-1
  ENDIF 

  st = read_idlstruct(file, error=error)
  IF error EQ 0 THEN status=0
  return,st

END 


FUNCTION huan_goods::seeing_file
  seeing_file = 'seeing_data/'+self.type+'_seeing.dat'
  seeing_file = concat_dir(self->dir(), seeing_file)
  return,seeing_file
END 


PRO huan_goods::read_seeing, status=status

  status = 1

  seeing_file = self->seeing_file()
  IF NOT fexist(seeing_file) THEN BEGIN 
      message,'The seeing file for type: '+self.type+' was not found: ',/inf
      message,seeing_file,/inf
      return
  ENDIF 

  ;; read the seeing data
  print,'Reading seeing file: ',seeing_file
  readcol, seeing_file, seeing, format='F'

  self.seeing = ptr_new(seeing, /no_copy)
  status = 0
  return

END 

FUNCTION huan_goods::seeing

  IF n_elements( *self.seeing ) NE 0 THEN BEGIN 
      return,*self.seeing
  ENDIF ELSE BEGIN 
      return,-1
  ENDELSE 

END 

;; Wrapper for the stand-alone procedure huan_goods_admom

PRO huan_goods::run_admom

  dir = self->dir()
  type = self.type
  seeing = self->seeing()
  nseeing=n_elements(seeing)
  
  secstr = ['22','23','33','34']
  FOR i=0L, nseeing-1 DO BEGIN 

      seestr = ntostr(seeing[i], 4, /round)

      catlist = concat_dir(dir, $
        'cat/'+type+'_'+seestr+'_h_goods_si_sect'+secstr+'_r1.0z_cat_mod.fit')
      imlist  = concat_dir(dir,  $
        'images/'+type+'_'+seestr+'_h_si_sect'+secstr+'_v1.0_drz_img.fits')
      
      huan_goods_admom, imlist, catlist

  ENDFOR 

END 


;; Wrapper for stand-alone procedure  huan_goods_sensitivity_convolved
;; Should incorporate that into this class
PRO huan_goods::calc_sensitivity

  seeing = self->seeing()
  type = self.type
  outDir = concat_dir(self->dir(),'sensitivity/')
  outFile = concat_dir(outDir, type+'_convolved.st')

  nseeing = n_elements(seeing)
  FOR i=0L, nseeing-1 DO BEGIN 

      tsens_struct = $
        huan_goods_sensitivity_convolved(seeing[i], $
                                         type=type, $
                                         noprompt=noprompt)

      IF i EQ 0 THEN BEGIN 
          sens_struct = replicate(tsens_struct, nseeing)
      ENDIF ELSE BEGIN 
          sens_struct[i] = tsens_struct
      ENDELSE 

  ENDFOR 

  print
  print,'Writing to file: ',outfile
  write_idlstruct, sens_struct, outfile, /ascii


END 


FUNCTION huan_goods::color_list

  IF !d.name EQ 'X' THEN BEGIN 
      colors = [!p.color, !green, !magenta, !yellow, !DodgerBlue, !DarkGreen]
  ENDIF ELSE BEGIN 
      colors = [!p.color, !blue, !red, !DarkGreen, !magenta, !dodgerBlue]
  ENDELSE 
  return,colors

END 

PRO huan_goods::plot_sensitivity, type=type, set=set, dops=dops

    if n_elements(seeing_range) eq 0 then seeing_range=[0.0,1.5]
    
    ;; This set number determines the range of seeing output    
    if n_elements(set) eq 0 then set=1 
    
    IF n_elements(type) EQ 0 THEN type = self.type

    if keyword_set(dops) then begin 
        dir='~/plots/des_sensitivity'
        file = type+'_convolved_set'+ntostr(set)+'.eps'
        file=concat_dir(dir,file)
        begplot,file,/encapsulated, /color
    endif
   
    case set of
        1: seeing_range = [0.0, 1.5]
        2: seeing_range = [0.5, 1.25]
        else: message,'Unsupported set number'+ntostr(set)
    endcase
    
    st = self->read_sensitivity(status=status, type=type)
    IF status NE 0 THEN return

    colors = self->color_list()

    IF type EQ 'des5yr' THEN BEGIN 
        names = ['des5yr 3 band (goods)', 'des5yr 1 band (goods)']

        med_seeing1 = st.med_seeing
        effdens1 = st.weffdens3

        
        med_seeing2 = st.med_seeing
        effdens2 = st.effdens

        w=where(med_seeing1 gt seeing_range[0] and $
                med_seeing1 lt seeing_range[1], nw)
        if nw eq 0 then message,'seing range too restrictive'
        med_seeing1 = med_seeing1[w]
        effdens1 = effdens1[w]
        med_seeing2 = med_seeing2[w]
        effdens2 = effdens2[w]

        psym = [8, 4]

    ENDIF ELSE IF type EQ 'cfh' THEN BEGIN 
        names = ['CFHT 900s (goods)', 'CFHT 900s']

        med_seeing1 = st.med_seeing
        effdens1 = st.effdens

        w=where(med_seeing1 gt seeing_range[0] and $
                med_seeing1 lt seeing_range[1], nw)
        if nw eq 0 then message,'seing range too restrictive'
        med_seeing1 = med_seeing1[w]
        effdens1 = effdens1[w]

        med_seeing2 =  [0.626836]
        effdens2 = [7.035764]
        
        symsize2 = 2

        psym = [8,2]
    ENDIF ELSE message,'Do not know the type'

    nnames = n_elements(names)
    colors = colors[0:nnames-1]

    xtitle = 'seeing [arcsec]'
    ytitle = 'Effective Density [#/arcmin!U2!N]'


    xrange = [0.9*seeing_range[0],1.1*seeing_range[1]] 
  
    pplot,med_seeing1,effdens1, psym=psym[0], $
        xrange=xrange, xstyle=1, $
        xtitle=xtitle, ytitle=ytitle, aspect=1
    IF nnames GE 2 THEN BEGIN 
        oplot,med_seeing2, effdens2, psym=psym[1],color=colors[1],$
            symsize=symsize2
    ENDIF 

    legend,names,color=colors,psym=psym, /right, box=0

    IF n_elements(dops) NE 0 THEN BEGIN 
        endplot, /trim_bbox
    ENDIF 

END 

PRO huan_goods::plot_dens

  

END 



PRO huan_goods::free

  ptr_free, self.seeing
  self.type = ''
  
END 


;; Create scripts to run the external procedure huan_goods_addnoise
PRO huan_goods::addnoise_scripts

  oseeing = self->seeing()
  ns = n_elements(oseeing)

  ;; split into 4 sections, each for one cheops
  neach = ns/4

  scriptDir = '/net/cheops2/home/esheldon/idlscripts/goods_addnoise/'

  FOR scriptnum=1,4 DO BEGIN 

      scriptName = scriptDir + self.type+'_goods_addnoise'+ntostr(scriptnum)+'.sh'

      CASE scriptnum OF
          1: seeing = oseeing[0:neach-1]
          2: seeing = oseeing[neach:2*neach-1]
          3: seeing = oseeing[2*neach:3*neach-1]
          4: seeing = oseeing[3*neach:ns-1]
          ELSE: message,'Unknown script number: '+ntostr(scriptnum)
      ENDCASE 

      nseeing = n_elements(seeing)

      openw, lun, scriptName, /get_lun

      printf,lun,'#!/bin/bash'
      printf,lun
      
      FOR i=0L, nseeing-1 DO BEGIN 
          
          sstr = ntostr(seeing[i])
          
          printf,lun,'idl<<EOF'
          printf,lun,'  huan_goods_addnoise, "'+self.type+'", seeing='+sstr
          printf,lun,'EOF'
          printf,lun
      ENDFOR 
      
      free_lun, lun
      spawn,['chmod','755',scriptName],/noshell

  ENDFOR 

END 




FUNCTION huan_goods::cleanup

  ptr_free, self.seeing
  return,1

END 

PRO huan_goods__define

  struct = { $
             huan_goods, $
             type: '', $
             seeing: ptr_new() $
           }

END 
