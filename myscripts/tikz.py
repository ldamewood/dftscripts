import numpy

#plot(x[idx],dos[0,idx],x[idx],-dos[1,idx])

def tikz_dos(x, dos, efermi, elim = [-0.2204959523877739,0.2204959523877739]):
    out = []
    data = numpy.vstack([x-efermi,dos]).T
    idx = numpy.logical_and(x-efermi>=elim[0],x-efermi<=elim[1])
    if data.shape[1] == 2:
        out.append("e\ttot+\ttot-\n")
    else:
        out.append("e\ttot\n") 
    for i in data[idx]:
        out.append("\t".join(map(str,i)))
    return out