# fuzzy-onions - FlyBase scRNAseq scripts
# Copyright © 2021 Damien Goutte-Gattat
#
# This file is part of the Fuzzy-onions project and distributed under
# the terms of the MIT license. See the LICENSE.md file in that project
# for the detailed conditions.

from setuptools import setup
from fbcam.fuzzyonions import __version__

setup(
    name='fuzzy-onions',
    version=__version__,
    description='FlyBase scRNAseq scripts',
    author='Damien Goutte-Gattat',
    author_email='dpg44@cam.ac.uk',
    classifiers=[
        'Development Status :: 1 - Planning',
        'Environment :: Console',
        'Programming Language :: Python :: 3.9',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    install_requires=['requests', 'click_shell', 'pandas', 'IPython'],
    packages=['fbcam', 'fbcam.fuzzyonions'],
    entry_points={'console_scripts': ['fzo = fbcam.fuzzyonions.main:main']},
    command_options={
        'build_sphinx': {
            'version': ('setup.py', __version__),
            'release': ('setup.py', __version__),
        }
    },
)
