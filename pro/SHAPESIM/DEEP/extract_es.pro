pro extract_es,im,cat,parfile=parfile
;just a wrapper for running deep_pipe 
;on one image

if n_params() eq 0 then begin
	print,'-syntax extract_es,im,cat,parfile=parfile'
	return
endif

file='test.fits'
write_dir='/usr/users/esheldon/'

writefits,write_dir+file,im
read_dir=write_dir

if n_elements(parfile) eq 0 then begin
	parfile='/usr/users/esheldon/idl.lib/deep_es.par'
endif

deep_pipe_es, file,write_dir,read_dir,parfile,/nosmear
cat=mrdfits(read_dir+'test_sobj.fits',1,hdr)

return
end



