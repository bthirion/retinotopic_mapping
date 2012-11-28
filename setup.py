#! /usr/bin/env python
#
# Copyright (C) 2012 Bertrand Thirion <bertrand.thirion@inria.fr>

descr = """retinotopic mapping: Tools for analyzing retinotopy datasets."""

import os

DISTNAME = 'retinotopic_mapping'
DESCRIPTION = descr
MAINTAINER = 'Bertrand Thirion'
MAINTAINER_EMAIL = 'bertrand.thirion@inria.fr'
LICENSE = 'BSD (3-clause)'
DOWNLOAD_URL = 'git@github.com:bthirion/retinotopic_mapping.git'
VERSION = '0.1.dev'

from numpy.distutils.core import setup


if __name__ == "__main__":
    if os.path.exists('MANIFEST'):
        os.remove('MANIFEST')

    setup(name=DISTNAME,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        license=LICENSE,
        version=VERSION,
        download_url=DOWNLOAD_URL,
        packages=['retino']    
)
