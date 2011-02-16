function m2lfunc, r, p
    index = p[2]
    ratio = (r/p[1])^index
    func = p[0]*( ratio/(1+ratio) )
    return, func
end
