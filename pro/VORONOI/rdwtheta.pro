PRO rdwtheta, file, str

IF n_params() LT 1 THEN BEGIN
    print,'-Syntax: rdwtheta, file, str'
    return
ENDIF 

s = create_struct('meanr', 0., $
                  'area', 0., $
                  'shear', 0., $
                  'shearerr', 0., $
                  'realngal', 0., $
                  'realerr', 0., $
                  'randngal', 0., $
                  'randerr', 0., $
                  'relative', 0., $
                  'relativerr', 0.)

readcol,file,meanr,area,shear,shearerr,$
        realngal,realerr,randngal,randerr,relative,relativerr

n=n_elements(meanr)

str = replicate(s, n)

str.meanr = meanr
str.area = area
str.shear = shear
str.shearerr = shearerr
str.realngal = realngal
str.realerr = realerr
str.randngal = randngal
str.randerr = randerr
str.relative = relative
str.relativerr = relativerr

return
END 
