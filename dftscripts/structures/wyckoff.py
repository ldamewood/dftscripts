import numpy

def wyckoff(spgrp, site, x = 0., y = 0., z = 0., coordinate = 0):
    coords = 230*[{}]
    coords[129] = {
    	'2a' : [[.0,.0,.0],[.5,.5,.0]],
        '2b' : [[.0,.0,.5],[.5,.5,.5]],
        '2c' : [[.0,.5,numpy.mod(z,1.0)],[.5,.0,numpy.mod(-z,1.0)]],
    }
    coords[216] = {
        '4a' : [[.00,.00,.00]],
        '4b' : [[.50,.50,.50]],
        '4c' : [[.25,.25,.25]],
        '4d' : [[.75,.75,.75]],
    }
    return coords[spgrp][site][coordinate]