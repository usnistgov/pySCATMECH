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

scatmech_headers  = ['scatmech/allrough.h',
                     'scatmech/crossrcw.h',
                     'scatmech/flake.h',
                     'scatmech/matrix3d.h',
                     'scatmech/polydisperse.h',
                     'scatmech/roughnes.h',
                     'scatmech/subbobvlieg.h',
                     'scatmech/zernike.h',
                     'scatmech/askuser.h',
                     'scatmech/crough.h',
                     'scatmech/focussedbeam.h',
                     'scatmech/matrixmath.h',
                     'scatmech/random.h',
                     'scatmech/scateval.h',
                     'scatmech/subsphere.h',
                     'scatmech/zernikeexpansion.h',
                     'scatmech/axifree.h',
                     'scatmech/dielfunc.h',
                     'scatmech/fresnel.h',
                     'scatmech/miescat.h',
                     'scatmech/raygscat.h',
                     'scatmech/scatmatrix.h',
                     'scatmech/tmatrix.h',
                     'scatmech/axipart.h',
                     'scatmech/diffuse.h',
                     'scatmech/gcross.h',
                     'scatmech/mueller.h',
                     'scatmech/rayinst.h',
                     'scatmech/scatmech.h',
                     'scatmech/torrspar.h',
                     'scatmech/bobvlieg.h',
                     'scatmech/facet.h',
                     'scatmech/grating.h',
                     'scatmech/nsphere.h',
                     'scatmech/rayscat.h',
                     'scatmech/scattabl.h',
                     'scatmech/transmit.h',
                     'scatmech/brdf.h',
                     'scatmech/fft.h',
                     'scatmech/inherit.h',
                     'scatmech/oasphere.h',
                     'scatmech/raystack.h',
                     'scatmech/sizedistribution.h',
                     'scatmech/two_source.h',
                     'scatmech/coatedmie.h',
                     'scatmech/filmtran.h',
                     'scatmech/instrument.h',
                     'scatmech/onelayer.h',
                     'scatmech/rcw.h',
                     'scatmech/sphdfct.h',
                     'scatmech/twoface.h',
                     'scatmech/crossgrating.h',
                     'scatmech/finiteaperture.h',
                     'scatmech/lambert.h',
                     'scatmech/optconst.h',
                     'scatmech/reflectance.h',
                     'scatmech/sphprt.h',
                     'scatmech/urough.h',
                     'scatmech/crossgrating2.h',
                     'scatmech/firstdiffuse.h',
                     'scatmech/local.h',
                     'scatmech/phasefunction.h',
                     'scatmech/rough.h',
                     'scatmech/sphrscat.h',
                     'scatmech/vector3d.h']
                     

module1 = Extension('SCATPY',
                    include_dirs = ['scatmech'],
                    library_dirs = [],
                    libraries = [],
                    sources = ['SCATPYmodule.cpp'] + scatmech_source,
                    depends = scatmech_headers,
                    extra_compile_args = extra_compile_args,
                    extra_link_args =extra_link_args)

#setuptools.setup(
setup(
    name="pySCATMECH",
    version="0.1.5",  # NOTE: If you change this, change it in index.rst AND conf.py too!
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
    ],
    python_requires='>=3.6'
)
