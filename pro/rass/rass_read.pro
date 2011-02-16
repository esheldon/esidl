function rass_read_stdef, type, rassnames=rassnames

	stdef={ $
		rxs:'1RXS',$
		srcname:'J111639.4+871202',$
		ra:0d,$
		dec:0d,$
		pos_err:0,$
		flags:'F....',$
		flags2:'....',$
		cps:0.,$
		cps_err:0.,$
		bgr_cpsa:0., $
		exp:0,$
		hr1:0.,$
		hr1_err:0.,$
		hr2:0.,$
		hr2_err:0., $
		ext:0, $
		extl:0, $
		srcl:0, $
		extr:0, $
		priflge:'111111b',$
		vigf:0., $
		orgdat:'960531',$
		moddat:'000000',$
		id:0, $
		field_id_src_num:' 33001001_0013'}


	if strlowcase(type) eq 'fsc' then begin
		stdef_add = $
			{rct:0, $
			itb:'900916.08',$
			ite:'900921.14', $
			rl:0}

		stdef = create_struct(stdef, stdef_add)
	endif


	if keyword_set(rassnames) then begin
		; add rass_ to the front of all names.  This is useful when the
		; structure will be incorporate to another existing structure
		stdefold = stdef
		tags=tag_names(stdefold)
		stdef = create_struct('rass_'+tags[0], stdefold.(0))
		for i=1,n_elements(tags)-1 do begin
			stdef = create_struct(stdef, 'rass_'+tags[i], stdefold.(i))
		endfor
	endif


	return, stdef

end
function rass_read, type, rows=rows, columns=columns, rassnames=rassnames

	rass_dir = getenv('RASS_DIR')


	case strlowcase(type) of
		'fsc': begin
			name = getenv('RASS_FSC_NAME')
			skiplines=175
		end
		'bsc': begin
			name = getenv('RASS_BSC_NAME')
			skiplines=4
		end
		'all': begin

			bsc = rass_read('bsc', rassnames=rassnames)
			columns = tag_names(bsc)

			fsc = rass_read('fsc', columns=columns, rassnames=rassnames)

			concat_structs, fsc, bsc, struct

			newtag='cat'
			if keyword_set(rassnames) then begin
				newtag='rass_cat'
			endif
			struct = struct_addtags(struct, newtag, "''")

			tagnum=n_tags(struct)-1
			struct[0:n_elements(fsc)-1].(tagnum) = 'fsc'
			struct[n_elements(fsc):n_elements(struct)-1].(tagnum) = 'bsc'
			
			return, struct
		end
		else: message,'Unkown catalog type: '+string(type)
	endcase

	name = filepath(root=rass_dir, name)
	nl=file_lines(name)

	stdef = rass_read_stdef(type, rassnames=rassnames)
	struct = ascii_read(name, stdef, nl, $
		skiplines=skiplines, delim=' ', rows=rows, columns=columns)

	return, struct

end
