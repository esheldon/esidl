;+
; NAME:
;   datadir
;
; PURPOSE:
;   The datadir.pro and datafile.pro functions implement a generic data and file naming structure for data projects.
;   The basic structure is that of a base type directory within a project, such 
;   as 'output', and a file type.  These are augmented with sample 
;   directories, subsets sample names and extra path elements.
;
; CALLING SEQUENCE:
;   dir = datadir(basetype, filetype, root=, project=, sample=, subtype=,
;                   pathextra=, /createdir)
;
; INPUTS:
;   basetype: A base type directory, such as 'output' or 
;       'input' or 'simulations'
;   filetype: A filetype such as 'raw', 'corrected', or 'log'
;
; OPTIONAL INPUTS:
;   root: The root directory. Default is '' in which case all directories
;       are relative to the current working directory.
;   project: A project name.  This directory sits right above the root.
;       Default is ''
;   sample: A sample name, such as 'sample12', or 'bright-galaxies'. This
;       sits above the root and project.
;   subtype: A subtype name.  Specification of a subtype adds a separate
;       heirarchy to the directory structure, so that the filetype dir
;       will come after 'sub/subtype'
;   pathextra: Extra path elements, possibly an array, to be added after
;       all previous paths elements.
;
; KEYWORD PARAMETERS:
;   /createdir: Create the directory.
;
; OUTPUTS:
;   The directory.
;
; EXAMPLE:
;   print, datadir('output', 'corrected', root='~/data', project='correlate', 
;                   sample='sample12')
;   ~/data/correlate/output/sample12/corrected
;
;   Include subsamples.
;   print, datadir('output', 'corrected', root='~/data', project='correlate', 
;                   sample='sample12', sub='ngals200_12')
;   ~/data/correlate/output/sample12/sub/ngals200_12/corrected/
;
; MODIFICATION HISTORY:
;   Created: 2007-04-23, Erin Sheldon, NYU
;
;-

function datadir, btype, ftype, root=root, project=project, sample=sample, subtype=subtype, pathextra=pathextra, createdir=createdir

    tf = datafile(btype, ftype, root=root, project=project, sample=sample, subtype=subtype, pathextra=pathextra, createdir=createdir, dir=dir)
    return, dir
end
