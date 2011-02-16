;+
; NAME:
;  ESHELDON_CONFIG
;
; PURPOSE:
;  Return the value of a configuration variable
;
; CALLING SEQUENCE:
;  value = esheldon_config(varName, /struct, /reload, exists=)
;
; INPUTS:
;  varName: The name of the config. variable to get.
;
; KEYWORD PARAMETERS:
;  /struct: Return the structure containing all the configuration variables. In
;           this case the user need not send any arguments.
;  /reload: Reload the config file.
;
; OUTPUTS:
;  By default, either the value of the configuration varaible or -1 if it does
;    not exist. 
;  If /struct, then the configuration structure, or -1 if it cannot be read. 
;
; OPTIONAL OUTPUTS:
;  exists=exists:  1 if the variable exists or 0 if not
;
; COMMON BLOCKS:
;  esheldon_config_block, configStruct, tags
;
; SIDE EFFECTS:
;  If the config file has not been loaded, then esheldon_load_config is run and
;  the common block above is modified.
;
; EXAMPLE:
;  ;; Read the sdss data_dir config variable
;  data_dir = esheldon_config('data_dir')
;  
;  ;; List the available configuration variables
;  help,esheldon_config(/struct),/str
;
;  ;; Check for a variable and if it is defined, continue
;  lss_dir = esheldon_config('lss_dir', exists=lss_dir_exists)
;  if lss_dir_exists then begin 
;      file = lss_dir + '......'
;
;
; MODIFICATION HISTORY:
;  Created: 14-Jan-2005, Erin Sheldon, UChicago
;
;-

function _esheldon_config_constants

	rhounits = double(1.e18)
	rhocrit = double(2.77545e-7) ; Msolar/pc^3
	rhocrit = rhocrit*rhounits ; Msolar/Mpc^3


	sunmag = [6.38d,5.06d,4.64d,4.53d,4.52d] ; sun absmag
	mstar= [-18.34d, -20.04d, -20.83d, -21.26d, -21.55d]
	mstar_err= [0.08d, 0.04d, 0.03d, 0.04d,0.04d]
	constants = { $
		bands:['u','g','r','i','z'], $
		$
		sunmag: sunmag, $
		mstar: mstar, $
		mstar_err: mstar_err, $
		lstar: 10d^( (mstar - sunmag)/(-2.5d) ), $
		d2r: !dpi/180.0d0, $
		r2d:180.0d0/!dpi, $
		plusminus:string(177b), $ ; ASCII plus or minus
		shapenoise:0.32, $
		rhocrit: rhocrit, $
		$
		hardrcut: 0.8, $
		hardprobcut: 0.8, $
		gratio: (1d + sqrt(5d) )/2d, $
		flags_masked_basic:'200'x, $
		flags_masked_simple:'400'x, $
		flags_masked_bound:'1000'x, $
		flags_masked_combined:'2000'x $
	}

	notused={ $
		$
		$ ; sdsspix boundary flags
		minlamflag: 2b^0, $
		maxlamflag: 2b^1, $
		minetaflag: 2b^2, $
		maxetaflag: 2b^3, $
		$
		flags_masked: '1'x, $
		flags_quad1_masked: '2'x, $
		flags_quad2_masked: '4'x, $
		flags_quad3_masked: '8'x, $
		flags_quad4_masked: '10'x, $
		$
		flags_quad1_masked_monte: '20'x, $
		flags_quad2_masked_monte: '40'x, $
		flags_quad3_masked_monte: '80'x, $
		flags_quad4_masked_monte: '100'x, $
		$
		flags_masked_basic:'200'x, $
		flags_masked_simple:'400'x, $
		flags_masked_bound:'1000'x, $
		flags_masked_combined:'2000'x, $
		$
		prob_psf_default: -9999., $
		probflag_noiter_requested: 2b^0, $
		probflag_noiter: 2b^1, $
		probflag_noprob: 2b^2 }

	return, constants

end

function esheldon_config, varname, struct=struct, reload=reload, exists=exists

	common esheldon_config_block, configstruct, tags

	if (n_elements(varname) eq 0) and (not keyword_set(struct)) then begin 
        on_error,2
        print,'usage: esheldon_config(varname, /struct, /reload, /exists)'
        message,'Halting'
	endif 

	; Make sure the configuration is loaded

	if n_elements(configStruct) eq 0 or keyword_set(reload) then begin 
		esheldon_load_config
		if n_elements(configStruct) eq 0 then begin 
            on_error, 2
			message,'Could not load config file'
		endif 

		constants = _esheldon_config_constants()
		configStruct = create_struct(configStruct, constants)
		tags = tag_names(configStruct)
	endif 

	if keyword_set(struct) then begin 
		exists = 1
		return,configStruct
	endif 

	; Check for the input config var
	checkName = strupcase(string(varName[0]))

	w=where(tags eq checkName, nw)
	if nw eq 0 then begin 
        message,'config variable not found: "'+checkName+'"'
	endif else begin 
		exists = 1
		return,configStruct.(w[0])
	endelse 

end 
