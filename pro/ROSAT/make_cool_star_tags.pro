pro make_cool_star_tags,taglist

if n_params() eq 0 then begin
	print,'-syntax make_cool_star_tags, taglist'
	return
endif

taglist = ['PARENT','NCHILD','ID','OBJC_TYPE','OBJC_FLAGS',$
	'OBJC_ROWC','OBJC_ROWCERR','OBJC_COLC','OBJC_COLCERR',$
	'ROWC','COLC','FIBERCOUNTS','FIBERCOUNTSERR',$
	'PSFCOUNTS','PSFCOUNTSERR','PETRORAD',$
	'PETRORADERR','PETROR50','PETROR50ERR','PETROR90','PETROR90ERR',$
	'Q','QERR','U','UERR','FLAGS','RA','DEC','ISO_A','ISO_B','ISO_PHI']

return
end
