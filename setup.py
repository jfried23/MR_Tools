from distutils.core import setup
from Cython.Build import cythonize
import numpy
import scipy


setup(
    ext_modules = cythonize(["./src/Sim/Bloch/BlochSim.pyx",
			     "./src/Sim/Bloch/PulseSeq/PulseElem.pyx",

			     "./src/Sim_cython/Spin.pyx",
			     "./src/Sim_cython/Bloch/BlochSim.pyx",
			     "./src/Sim_cython/Bloch/PulseSeq/PulseElem.pyx" 
			     ]),

    include_dirs=[numpy.get_include()]
)

