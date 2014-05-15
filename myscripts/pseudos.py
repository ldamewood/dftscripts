# http://www.abinit.org/downloads/PAW2/JTH-TABLE/ATOMICDATA/H.GGA_PBE-JTH-paw.xml

import os
import urllib2

from pymatgen.core.periodic_table import Element

__all__ = [
    'get_paw'
]

def get_paw(symbol):
    
    if isinstance(symbol, Element):
        symbol = symbol.symbol

    dirpath = os.path.join(os.getenv("HOME"), ".abinit", "paw")
    urlpath = 'http://www.abinit.org/downloads/PAW2/JTH-TABLE/ATOMICDATA/'
    pawfilename = "%s.GGA_PBE-JTH-paw.xml" % symbol
    
    url = urllib2.urlparse.urljoin(urlpath,pawfilename)
    path = os.path.join(dirpath, pawfilename)
    
    if not os.path.exists(path):
        response = urllib2.urlopen(url)
        with open(path, 'w') as fh:
            fh.write(response.read())
    return path