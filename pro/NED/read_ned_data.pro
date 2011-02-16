pro read_ned_data,file,struct
;
;+
; NAME:
;    READ_NED_DATA
; PURPOSE:
;    This program writes reads in a cleaned NED output file.
;-
 if N_params() eq 0 then begin
	print,'read_ned_data,file,struct'
	return
 endif

 get_lun,flun
 openr, flun, file

 string=''

 t=create_struct('type','','name','','ra',0.0,'dec',0.0,'z',0.0,'mag',0.0)
 struct=t
 while not eof(flun) do begin
	readf,flun,string,format='(A80)'
	ntype=strtrim(strmid(string,0,7),2)
	nname=strtrim(strmid(string,7,17),2)
	nra=strmid(string,24,11)
	ndec=strmid(string,36,11)
	nz=float(strtrim(strmid(string,49,6),2))

	nra_h=float(strmid(nra,0,2))
	nra_m=float(strmid(nra,3,2))
	nra_s=float(strmid(nra,6,4))
	;Don't forget conversion to degrees!
	ra_new=ten(nra_h,nra_m,nra_s)*15.0

	ndec_h=float(strmid(ndec,0,3))
	ndec_m=float(strmid(ndec,4,2))
	ndec_s=float(strmid(ndec,7,2))
	dec_new=ten(ndec_h,ndec_m,ndec_s)
	;Handle the sign bug....
	if (strmid(ndec,0,1) eq '-' and ndec_h eq 0.0) then begin
		dec_new=dec_new*(-1.0)
	endif

	nt=t
	nt.type=ntype
	nt.name=nname
	nt.ra=ra_new
	nt.dec=dec_new
	nt.z=nz

	readf,flun,string, format='(A80)'	
	nmag=strtrim(strmid(string,60,5),2)
	if (strmid(nmag,4,1) ne '/') then begin
	 if(strmid(nmag,0,1) ne '>') then begin
	  if (nmag ne '/') then begin
	    nmag=float(nmag)
	  endif else begin
	    nmag=0.0
	  endelse
	 endif else begin
	   nmag=float(strmid(nmag,4,1))
	 endelse
	endif else begin
	   nmag=float(strmid(nmag,0,4))
	endelse

	nt.mag=nmag
	struct=[struct,nt]

 endwhile

 close,flun
 free_lun,flun


 return
 end
 



;1st Line:
;
;HEADING    : CONTENTS
;-------      --------
;Type       : NED "Preferred" object type
;Object Name: One name of the object, continued on 2nd line if necessary
;Position   : Coordinates in your chosen equinox
;vel or z   : Radial velocity, or redshift if in square brackets.
;Ref        : The number of references in NED for this object
;Field 1    : NED Basic Data; varies according to "Type" on this line.
;             See below.
;
;
;2nd Line:
;
;HEADING    : CONTENTS
;-------      --------
;Object Name: Continued from 1st line if necessary
;Pos Unc    : Positional uncertainty, expressed as the major and minor axes
;             of the 95% confidence ellipse, in arc seconds
;PA         : Position angle of the uncertainty ellipse, in degrees, for B1950.0.
;unc        : Velocity or redshift uncertainty, expressed as the standard
;             deviation, in the units of the field on the 1st line
;Pht        : Number of photometric data points
;Fld2       : NED Basic Data; varies according to "Type" on 1st line.  See below.
;Fld3       : NED Basic Data; varies according to "Type" on 1st line.  See below.
;Fld4       : NED Basic Data; varies according to "Type" on 1st line.  See below.

