# http://www.abinit.org/downloads/PAW2/JTH-TABLE/ATOMICDATA/H.GGA_PBE-JTH-paw.xml

import os
import urllib2

from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import IStructure
from collections import OrderedDict

__all__ = [ 'get_psp' ]

paw_urlpath = 'http://www.abinit.org/downloads/PAW2/JTH-TABLE/ATOMICDATA/'
paw_filenamepattern = '%(symbol)s.GGA_PBE-JTH-paw.xml'
paw_localpath = os.path.join(os.getenv("HOME"), ".abinit", "paw")

nc_urlpath = 'ftp://ftp.abinit.org/pub/abinitio/Psps/GGA_FHI/'
nc_filenamepattern = '%(Z)02d-%(symbol)s.GGA.fhi'
nc_localpath = os.path.join(os.getenv("HOME"), ".abinit", "nc")

def get_psp(element, method = 'paw'):
    
    def unique_order(x):
        return list(OrderedDict.fromkeys(x))
    
    if isinstance(element, IStructure):
        return list([get_psp(element) for element in unique_order(element.species)])
    
    if isinstance(element, list):
        return list([get_psp(Element(atom), method) for atom in set(element)])
    
    if method == 'paw':
        dirpath = paw_localpath
        urlpath = paw_urlpath
        filename = paw_filenamepattern % dict(Z = element.Z, symbol = element.symbol)
    elif method == 'nc':
        dirpath = nc_localpath
        urlpath = nc_urlpath
        filename = nc_filenamepattern  % dict(Z = element.Z, symbol = element.symbol)
    else:
        raise NotImplemented()
    
    path = os.path.join(dirpath, filename)
    
    if not os.path.exists(path):
        url = urllib2.urlparse.urljoin(urlpath,filename)
        response = urllib2.urlopen(url)
        with open(path, 'w') as fh:
            fh.write(response.read())
    return path