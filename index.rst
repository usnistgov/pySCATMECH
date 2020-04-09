.. pySCATMECH documentation master file, created by
   sphinx-quickstart on Tue Mar 17 07:29:04 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pySCATMECH : A Python Interface to the SCATMECH Library
=======================================================

pySCATMECH is a Python interface to the SCATMECH library of scattering codes.

pySCATMECH is a Python interface to many of the features in the
SCATMECH library. If you are unfamiliar with the SCATMECH library, be
sure to check out its documentation at `SCATMECH Polarized Light
Scattering C++ Class Library <https://pages.nist.gov/SCATMECH/docs/index.htm>`_.

To install, use ``pip install pySCATMECH`` at a command prompt.

:Author: Thomas A. Germer
:E-mail: thomas.germer@nist.gov
:Telephone: 301-975-2876
:Address:
   | Sensor Science Division
   | National Institute of Standards Technology
   | 100 Bureau Drive Stop 8443
   | Gaitherburg, MD 20899-8443
:Current Version: 0.0.1a0 (8 April 2020)

.. toctree::
   :maxdepth: 2
   :caption: Modules
	     
   mueller
   fresnel
   model
   brdf
   local
   integrate
   scatterer
   rcw
   crossrcw

.. toctree::
   :maxdepth: 2
   :caption: Examples

   UsingMueller
   UsingFresnel
   UsingBRDF
   UsingLocal
   UsingIntegrate
   UsingScatterer
   UsingRCW
   UsingCrossRCW

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
