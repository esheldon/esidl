;; Add photoid and make sure rerun is an int

FUNCTION _convert_datasweep_struct, instruct

  ;; Change rerun to an int16 from ..... string?
  ;; camcol and field only need be int16 as well
  alter_struct = $
    {rerun:0,$
     camcol:0,$
     field:0}
  newst = alter_tags(instruct[0], alter_struct)

;  newst = remove_tags(newst, 'aperflux')

  ;; add photoid
  newst = create_struct('photoid', 0LL, newst)

  num = n_elements(instruct)
  newst = replicate(newst, num)
  
  struct_assign, instruct, newst

  newst.photoid = long64( photoid(newst) )
  newst.rerun = fix( instruct.rerun )

  return, newst

END 


pro stuff_datasweeps, files, tmpdir=tmpdir, table=table, maxn=maxn

    dir = sdssidl_config('data_dir')
    dir = path_join(dir, subdir=['datasweep','dr6uber','137'])
    ;dir = '/global/data/sdss/redux/datasweep/dr4uber/137'

    files=file_search(concat_dir(dir,'calibObj-??????-?-gal.fits.gz'),count=nf)
  
    if n_elements(tmpdir) eq 0 then tmpdir = '~/tmp'
    pg = obj_new('postgres')

    ;conn = 'user=postgres'
    tm=systime(1)

    if n_elements(table) eq 0 then table = 'datasweep'

    if n_elements(maxn) ne 0 then begin
        if maxn lt nf then nf=maxn 
    endif
    for i=0l, nf-1 do begin 

        print,'Reading file: ',files[i]
        st = mrdfits(files[i], 1)

        outstruct = _convert_datasweep_struct(st)

        pg->struct2table, outstruct, table, $
            primary_key='photoid', $
            tmpdir=tmpdir

    endfor 
    ptime,systime(1)-tm

    pg->create_index, table, 'run,rerun,camcol,field,id'

;    pg->query, 'GRANT SELECT ON '+table+' TO sdss', conn=conn, $
;        status=status
;    if status ne pg->status_val('no_result') then message,'Could not grant select to sdss'
    obj_destroy, pg

end 
