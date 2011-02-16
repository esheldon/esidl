pro pbs_bnl_idl, pbs_file, idl_commands, job_name=name, queue=queue, walltime=walltime, nodes=nodes, ppn=ppn, setup=setup, _extra=_extra, gdl=gdl

	default_setup="source /home/users/esheldon/.dotfiles/bash/riemann/idl_setup.sh &&  source /home/users/esheldon/.dotfiles/bash/riemann/gdl_setup.sh"

	if n_elements(pbs_file) eq 0 or n_elements(idl_commands) eq 0 then begin
		on_error,2
		print,'Usage: pbs_bnl_idl, pbs_file, idl_commands, job_name=, queue="batch", walltime="24:00:00", nodes=1, ppn=1, setup=, _extra='
		print,'The values for all options must be expressed as strings'
		print,'Send as many -l options as you want through _extra, e.g'
		print,"  cput='35:00'"
		print
		print,"if setup= is not sent then the following command is issued: "
		print,'  '+default_setup
		message,'Returning to caller'
	endif

	if n_elements(nodes) eq 0 then nodes_str='1' else nodes_str=strn(nodes)
	if n_elements(ppn) eq 0 then ppn_str='1' else ppn_str=strn(ppn)
	;if n_elements(queue) eq 0 then queue='batch'
	if n_elements(walltime) eq 0 then walltime='24:00:00'

	if n_elements(setup) eq 0 then begin
		setup=default_setup
	endif
 
	openw, lun, pbs_file, /get

	printf,lun,"#PBS -l nodes="+nodes_str+":ppn="+ppn_str
	printf,lun,"#PBS -l walltime="+walltime

	nextra = n_tags(_extra)
	if nextra ne 0 then begin
		tn=tag_names(_extra)
		for i=0L, nextra-1 do begin
			printf,lun,"#PBS -l ",tn[i],"=",_extra.(i),f='(a,a,a,a)'
		endfor
	endif

	;printf,lun,"#PBS -q "+queue
	if n_elements(name) ne 0 then begin
		printf,lun,"#PBS -N "+string(name)
	endif

	printf,lun,"#PBS -j oe"
	; not useful since it only writes the file after the job finished
	printf,lun,"#PBS -o "+pbs_file+'.pbslog'
	printf,lun,"#PBS -m a"
	printf,lun,"#PBS -V"
	printf,lun,"#PBS -r n"
	printf,lun,"#PBS -W umask=0022"

	printf,lun,"echo Running on `hostname`"
	;printf,lun,"umask 0022"
	printf,lun
	printf,lun, "# my startup has dependencies I don't want for this"
	printf,lun, 'export IDL_STARTUP=""'
	printf,lun
	printf,lun,"# These are between && so if one fails all fail"
    for i=0L, n_elements(setup)-1 do begin
        printf,lun, setup[i]
    endfor
	printf,lun,"if [[ $? != 0 ]]; then"
	printf,lun,"    echo 'setup failed'"
	printf,lun,"    exit 45"
	printf,lun,"fi"
	printf,lun
	printf,lun,"logf="+pbs_file+'.log'
	;printf,lun,'idl_commands="'
	if keyword_set(gdl) then begin
		exec="gdl"
	endif else begin
		exec="idl"
	endelse
	printf,lun,exec+' &> "$logf" <<EOF'

	for i=0L, n_elements(idl_commands)-1 do begin
		printf,lun,"    "+strtrim(idl_commands[i],2)
	endfor
	;printf,lun,'"'
	printf,lun,'EOF'

	printf,lun

	;printf,lun,'echo "$idl_commands" | '+exec+' &> '+pbs_file+'.log'
	free_lun, lun

end
