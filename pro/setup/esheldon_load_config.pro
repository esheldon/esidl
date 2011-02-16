PRO esheldon_load_config, status=status

    status = 1

    COMMON esheldon_config_block, configStruct, tags

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; Read the config file
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    config_file = getenv("ESHELDON_CONFIG")

    if config_file[0] eq '' then begin 
        config_file = '~esheldon/.idl_config/esheldon.conf'
        message, "$ESHELDON_CONFIG undefined, trying "+config_file, /inf
    endif 

	config_file = expand_tilde(config_file)

    print,'Loading config file: ',config_file
    parse_config, config_file, keywords, values, struct=configstruct, $
        status=status

    if status ne 0 then begin 
        print,'Error loading esheldon config file'
        return
    endif

    tags = tag_names(configStruct)

end 
