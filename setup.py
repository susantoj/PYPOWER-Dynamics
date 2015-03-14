# Copyright (C) 2014-2015 Julius Susanto
#
# PYPOWER-Dynamics is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# PYPOWER-Dynamics is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PYPOWER-Dynamics. If not, see <http://www.gnu.org/licenses/>.

from setuptools import setup, find_packages

setup(
    name='pypower-dynamics',
    version='1.0.1',
    author='Julius Susanto',
    author_email='susanto@ieee.org',
    description='Time-domain simulation (transient stability) module for PYPOWER',
    long_description='PYPOWER-Dynamics',
    url='https://github.com/susantoj/PYPOWER-Dynamics',
    license='GPLv3',
    install_requires=[
        # Deactivated to avoid problems with system packages.
        # Manual installation of PYPOWER, NumPy and SciPy required.
        # 'numpy>=1.6',
        # 'scipy>=0.9',
    ],
    #packages=find_packages(exclude=['screenshots']),
    packages=['pydyn'],
    include_package_data=True,
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Topic :: Scientific/Engineering',
    ],
)
