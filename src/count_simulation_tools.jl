"""Binning x and y with binsz. 
    x and y must be of length 2^N
    binsz must be a power of 2
"""
function bin_xy(x, y, binsz)
    # raise error if the length of x and binsz are not a power of 2
    if !ispow2(binsz)
        error("binsz must be a power of 2")
    end
    if !ispow2(length(x))
        error("length of x must be a power of 2")
    end
    binindx = 1:binsz:length(x)
    xbin = x[binindx]
    indx = searchsortedlast.(Ref(xbin), x)
    ybinned = [sum(y[indx .== ii]) for ii in unique(indx)]
    xbinned = [mean(x[indx .== ii]) for ii in unique(indx)] # bin x as well

    return xbinned, ybinned
end