FUNCTION lrg_morad::init
  
  ;; No data initialization needed other than running setup to ensure
  ;; that the config file is read.
  sdssidl_setup

  defsysv, '!SDSSIDL_CONFIG_LOADED', exist=sloaded_exists
  IF NOT sloaded_exists THEN BEGIN 
      message,'!SDSSIDL_CONFIG_LOADED is not defined.  This probably means '+$
        'the config file was not loaded, and many of the methods will not '+$
        'work',/inf
  ENDIF ELSE BEGIN 
      IF NOT !SDSSIDL_CONFIG_LOADED THEN BEGIN 
          message,'The config file has not been loaded.  Many of the methods '+$
            'will not work',/inf
      ENDIF 
  ENDELSE 

  return,1

END 

FUNCTION lrg_morad::dir
  dir = esheldon_config('lensinput_dir')+'lrg_morad/'
  return,dir
END 
FUNCTION lrg_morad::file, random=random
  dir = self->dir()
  IF keyword_set(random) THEN BEGIN 
      file = dir + 'dr4plus-random-lrg-spectroscopic-for-Erin.fits'
  ENDIF ELSE BEGIN 
      file = dir + 'dr4plus-lrg-spectroscopic-for-Erin.fits'
  ENDELSE 
  return,file
END 
FUNCTION lrg_morad::get, status=status, random=random
  file = self->file(random=random)
  print,'Reading file: ',file
  lrg = mrdfits(file, 1, status=status)
  return,lrg
END  

PRO lrg_morad__define

  struct = { $
             lrg_morad, $
             dummy: 0 $
           }

END 
