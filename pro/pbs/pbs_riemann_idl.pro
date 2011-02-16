pro pbs_riemann_idl, $
        pbs_file, $
        idl_commands, $
        job_name=name, $
        queue=queue, $
        walltime=walltime, $
        setup=setup, $
        nodes=nodes, $
        ppn=ppn, $
        exec=exec, $
        _extra=_extra

	default_setup="source /home/esheldon/.dotfiles/bash/riemann/idl_setup.sh"
	if n_elements(pbs_file) eq 0 or n_elements(idl_commands) eq 0 then begin
		on_error,2
		print,'Usage: pbs_riemann_idl, pbs_file, idl_commands, job_name=, queue="batch", walltime=, setup=, _extra='
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
	if n_elements(queue) eq 0 then queue='batch'

	if n_elements(setup) eq 0 then begin
		setuplist=default_setup
	endif else begin
        setuplist=['source /clusterfs/riemann/software/itt/idl/bin/idl_setup.bash', $
                   'source /home/products/eups/bin/setups.sh',$
                   setup]
    endelse
 
	openw, lun, pbs_file, /get

    if n_elements(walltime) ne 0 then begin
        printf,lun,"#PBS -l walltime="+walltime
    endif

	nextra = n_tags(_extra)
	if nextra ne 0 then begin
		tn=tag_names(_extra)
		for i=0L, nextra-1 do begin
			printf,lun,"#PBS -l ",tn[i],"=",_extra.(i),f='(a,a,a,a)'
		endfor
	endif

	printf,lun,"#PBS -q "+queue
	if n_elements(name) ne 0 then begin
		printf,lun,"#PBS -N "+string(name)
	endif

	printf,lun,"#PBS -l nodes="+nodes_str+":ppn="+ppn_str
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
    for i=0L, n_elements(setuplist)-1 do begin
        printf,lun, setuplist[i]
    endfor
	printf,lun

	printf,lun,"logf="+pbs_file+'.log'
    if n_elements(exec) eq 0 then begin
        if keyword_set(gdl) then begin
            exec="gdl"
        endif else begin
            exec='/clusterfs/riemann/software/itt/idl70/bin/idl'
        endelse
    endif
	printf,lun,exec+' &> "$logf" <<EOF'

	for i=0L, n_elements(idl_commands)-1 do begin
		printf,lun,"    "+strtrim(idl_commands[i],2)
	endfor
	printf,lun,'EOF'
	printf,lun
	free_lun, lun


    return

end
