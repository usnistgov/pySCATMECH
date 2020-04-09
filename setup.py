from setuptools import setup, Extension
#from distutils.core import setup, Extension
import platform

with open("README.md", "r") as fh:
    long_description = fh.read()

print("platform.system() = ",platform.system())

if platform.system() == 'Windows':
    extra_compile_args = []
    extra_link_args =[]
if platform.system() == 'Linux':
    extra_compile_args = ['-Wno-deprecated-declarations']
    extra_link_args =[]
elif 'CYGWIN' in platform.system():
    extra_compile_args = ['-Wno-deprecated-declarations']
    extra_link_args =[]

scatmech_source = ['scatmech/allrough.cpp',
                   'scatmech/crossgrating2.cpp',
                   'scatmech/fresnel.cpp',
                   'scatmech/models.cpp',
                   'scatmech/rcw.cpp',
                   'scatmech/scatmatrix.cpp',
                   'scatmech/two_source.cpp',
                   'scatmech/askuser.cpp',
                   'scatmech/crossrcw.cpp',
                   'scatmech/gcross.cpp',
                   'scatmech/mueller.cpp',
                   'scatmech/reflectance.cpp',
                   'scatmech/scattabl.cpp',
                   'scatmech/twoface.cpp',
                   'scatmech/axifree.cpp',
                   'scatmech/crough.cpp',
                   'scatmech/grating.cpp',
                   'scatmech/nsphere.cpp',
                   'scatmech/reg_brdf.cpp',
                   'scatmech/sizedistribution.cpp',
                   'scatmech/urough.cpp',
                   'scatmech/axipart1.cpp',
                   'scatmech/dielfunc.cpp',
                   'scatmech/inherit.cpp',
                   'scatmech/oasphere.cpp',
                   'scatmech/reg_facet.cpp',
                   'scatmech/sphdfct.cpp',
                   'scatmech/vector3d.cpp',
                   'scatmech/axipart2.cpp',
                   'scatmech/diffuse.cpp',
                   'scatmech/instrument.cpp',
                   'scatmech/onelayer.cpp',
                   'scatmech/reg_instrument.cpp',
                   'scatmech/sphprt.cpp',
                   'scatmech/zernike.cpp',
                   'scatmech/axisym.cpp',
                   'scatmech/facet.cpp',
                   'scatmech/jmatrix.cpp',
                   'scatmech/phasefunction.cpp',
                   'scatmech/reg_lambert.cpp',
                   'scatmech/sphrscat.cpp',
                   'scatmech/zernikeexpansion.cpp',
                   'scatmech/bobvlieg1.cpp',
                   'scatmech/fft.cpp',
                   'scatmech/jvector.cpp',
                   'scatmech/polydisperse.cpp',
                   'scatmech/reg_local.cpp',
                   'scatmech/stokes.cpp',
                   'scatmech/bobvlieg2.cpp',
                   'scatmech/filmtran.cpp',
                   'scatmech/lambert.cpp',
                   'scatmech/random.cpp',
                   'scatmech/reg_rough.cpp',
                   'scatmech/subbobvlieg.cpp',
                   'scatmech/bobvlieg3.cpp',
                   'scatmech/finiteaperture.cpp',
                   'scatmech/local.cpp',
                   'scatmech/raygscat.cpp',
                   'scatmech/reg_sphrscat.cpp',
                   'scatmech/subsphere.cpp',
                   'scatmech/brdf.cpp',
                   'scatmech/firstdiffuse.cpp',
                   'scatmech/matrixmath.cpp',
                   'scatmech/rayinst.cpp',
                   'scatmech/rough.cpp',
                   'scatmech/tmatrix.cpp',
                   'scatmech/coatedmie.cpp',
                   'scatmech/flake.cpp',
                   'scatmech/matrixmath2.cpp',
                   'scatmech/rayscat.cpp',
                   'scatmech/roughnes.cpp',
                   'scatmech/torrspar.cpp',
                   'scatmech/crossgrating.cpp',
                   'scatmech/focussedbeam.cpp',
                   'scatmech/miescat.cpp',
                   'scatmech/raystack.cpp',
                   'scatmech/scateval.cpp',
                   'scatmech/transmit.cpp']
    
module1 = Extension('SCATPY',
                    include_dirs = ['scatmech'],
                    library_dirs = [],
                    libraries = [],
                    sources = ['SCATPYmodule.cpp'] + scatmech_source,
                    extra_compile_args = extra_compile_args,
                    extra_link_args =extra_link_args)

#setuptools.setup(
setup(
    name="pySCATMECH",
    version="0.0.1a",
    author="Thomas A. Germer",
    author_email="thomas.germer@nist.gov",
    description="A Python interface to the SCATMECH library",
    long_description=long_description,
    #long_description_content_type="text/markdown",
    url="https://github.com/thomas-germer/pyscatmech",
    packages=["pySCATMECH"],
    ext_modules=[module1],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: Public Domain",
        "Topic :: Scientific/Engineering :: Physics",
        "Operating System :: Microsoft :: Windows"
    ]
    #, python_requires='>=3.6'
)

