function datafiles::init, root=root, project=project, sample=sample, extension=extension
    if n_elements(root) eq 0 then root=''
    if n_elements(project) eq 0 then project=''
    if n_elements(sample) eq 0 then sample=''
    if n_elements(extension) eq 0 then extension='.fits'

    self.root = root
    self.project=project
    self.sample=sample
    self.extension=extension
    return, 1
end


; These can be overridden by an inheriting class or input above
function datafiles::root
    return, self.root
end
function datafiles::project
    return, self.project
end
function datafiles::sample
    return, self.sample
end
function datafiles::extension, type
    return, self.extension
end




function datafiles::dir, btype, ftype, extra1, extra2, extra3, extra4, root=root, project=project, sample=sample, subtype=subtype, createdir=createdir
    tf=self->file(btype, ftype, extra1, extra2, extra3, extra4, root=root, project=project, sample=sample, subtype=subtype, createdir=createdir, dir=dir)
    return, dir
end

function datafiles::file, btype, ftype, extra1, extra2, extra3, extra4, root=root, project=project, sample=sample, subtype=subtype, number=number, postfix=postfix, extension=extension, nodir=nodir, createdir=createdir, dir=dir

    ; Set the defaults

    ; The root directory. 
    if n_elements(root) eq 0 then root = self->root()

    ; A project name; dir element to come right after root
    if n_elements(project) eq 0 then project = self->project()

    ; The sample number and string
    if n_elements(sample) eq 0 then sample = self->sample()

    if n_elements(extension) eq 0 then extension=self->extension()

    files = datafile(btype, ftype, extra1, extra2, extra3, extra4, root=root, project=project, sample=sample, subtype=subtype, number=number, postfix=postfix, extension=extension, nodir=nodir, createdir=createdir, dir=dir)

    return, files
end




pro datafiles__define
    struct = {                 $
            datafiles,         $
            root: '',          $
            project: '',       $
            sample: '',        $
            extension: ''      $
        }
end
