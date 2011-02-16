;+
; NAME:
;   datafile
;
; PURPOSE:
;   The datadir.pro and datafile.pro functions implement a generic data and file naming structure for data projects.
;   The basic structure is that of a base type directory within a project, such 
;   as 'output', and a file type.  These are augmented with sample 
;   directories, subsets, file numbers (such as for bins) and extra path
;   elements. File names aggregate elements from the path to specify a fairly
;   descriptive name.  Although such a complex name is not necessary it does
;   provide some information if needed.
;
; CALLING SEQUENCE:
;   file = datafile(basetype, filetype, root=, project=, sample=, subtype=,
;                   postfix=, number=, extension=, 
;                   /nodir, /createdir, dir=)
;
; OPTIONAL INPUTS:  All inputs are optional.
;   paths: A filetype such as 'raw', 'corrected', or 'log'
;   root: The root directory. Default is '' in which case all directories
;       are relative to the current working directory.
;   project: A project name.  This directory sits right above the root.
;       Default is ''
;   sample: A sample name, such as 'sample12', or 'bright-galaxies'. This
;       sits above the root and project.
;   subtype: A subtype name.  Specification of a subtype adds a separate
;       heirarchy to the directory structure, so that the filetype dir
;       will come after 'sub/subtype'
;   postfix; A postfix string for the file.
;   number: File numbers, such as 3 or [1,2,3].  These are added as a
;       part of the postfix. If an array is sent then there may be 
;       multiple files returned.
;   extension: The file extension.  Default is '.fits' (I am an astronomer
;       after all).
;
; KEYWORD PARAMETERS:
;   /nodir: Don't add the directory to the file name.
;   /createdir: Create the directory.
;
; OUTPUTS:
;   The file name or names.
;
; OPTIONAL OUTPUTS:
;   dir: The directory where the file lives.
;
; EXAMPLE:
;   print, datafile(['output', 'corrected'], root='~/data', project='correlate', 
;                   sample='sample12', number=[0,1,2])
;   ~/data/correlate/sample12/output/corrected/correlate-sample12-output-corrected-00.fits
;   ~/data/correlate/sample12/output/corrected/correlate-sample12-output-corrected-01.fits
;   ~/data/correlate/sample12/output/corrected/correlate-sample12-output-corrected-02.fits
;
; MODIFICATION HISTORY:
;   Created: 2007-04-23, Erin Sheldon, NYU
;
;-
function datafile, ftypes, root=root, project=project, sample=sample, number=number, postfix=postfix, extension=extension, nodir=nodir, createdir=createdir, dir=dir

    ; Work through the optional elements one at a time
    if n_elements(ftypes) eq 0 then ftypes=''

    ; The root directory. 
    if n_elements(root) eq 0 then root=''
    ; A project name; dir element to come right after root
    if n_elements(project) eq 0 then project=''
    ; extra postfix on each file
    if n_elements(postfix) eq 0 then pfix='' else pfix=postfix
    ; The sample number and string
    if n_elements(sample) eq 0 then sample=''
    ; file extension, defaults to fits. Why not?
    if n_elements(extension) eq 0 then extension='.fits'
 

    ; now create inputs for datafile
    add_arrval, project, dfinputs
    add_arrval, sample,  dfinputs
    add_arrval, ftypes, dfinputs

    ; Finalize the postfix
    nnum = n_elements(number)
    if nnum ne 0 then begin
        ; numbers can also go into the postfix. Note if postfix above
        ; will be truncated to a scalar
        npostfix = string(number, format='(i02)')

        if pfix[0] eq '' then begin
            pfix = npostfix
        endif else begin
            pfix = pfix[0] + '-' + npostfix[*]
        endelse
    endif


    ; loop over all possible postfixes
    npost = n_elements(pfix)
    for i=0L, npost-1 do begin
        pf = pfix[i]

        tfile = dfile(root=root, dfinputs, postfix=pf, ext=extension, nodir=nodir, dir=dir)
        add_arrval, tfile, files
    endfor

    if keyword_set(createdir) then file_mkdir, dir

    return, files
end


