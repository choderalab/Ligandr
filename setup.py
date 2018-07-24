"""
Ligandr
Set up systems to superpose ligands on different conformations
"""
from setuptools import setup
import versioneer

DOCLINES = __doc__.split("\n")

setup(
    # Self-descriptive entries which should always be present
    name='Ligandr',
    author='choderalab/Hersh Gupta',
    description=DOCLINES[0],
    long_description="\n".join(DOCLINES[2:]),
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='LGPLv3',

    # Which Python importable modules should be included when your package is installed
    packages=['Ligandr', "Ligandr.tests"],

    # Optional include package data to ship with your package
    # Comment out this line to prevent the files from being packaged with your software
    # Extend/modify the list to include/exclude other items as need be
    package_data={'Ligandr': ["data/*.dat"]
                  }, install_requires=['mdtraj', 'numpy', 'openeye', 'requests', 'openmmtools', 'openmm'],

    # Additional entries you may want simply uncomment the lines you want and fill in the data
    # author_email='me@place.org',      # Author email
    # url='http://www.my_package.com',  # Website
    # install_requires=[],              # Required packages, pulls from pip if needed; do not use for Conda deployment
    # platforms=['Linux',
    #            'Mac OS-X',
    #            'Unix',
    #            'Windows'],            # Valid platforms your code works on, adjust to your flavor
    python_requires=">=3.5"  # Python version restrictions

    # Manual control if final package is compressible or not, set False to prevent the .egg from being made
    # zip_safe=False,

)
