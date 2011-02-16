function path_sep
	case !version.os_family of
		'unix': return, '/'
		'Windows': return, '\'
		'MacOS': return, ':'
		else: message,'Unsupported os_family: '+!version.os_family
	endcase
end
