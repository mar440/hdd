import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

mpi_compile_args = os.popen("mpic++ --showme:compile").read().strip().split(' ')
mpi_link_args    = os.popen("mpic++ --showme:link").read().strip().split(' ')

print "++++++++++++++++++++++++++++++++++++++++++++++++++++"
print mpi_compile_args
print "++++++++++++++++++++++++++++++++++++++++++++++++++++"
print mpi_link_args
print "++++++++++++++++++++++++++++++++++++++++++++++++++++"

#######################################
extra_compile_args = []
for ii in mpi_compile_args:
    extra_compile_args.append(ii)
extra_compile_args.append("-I/home/mar440/intel/mkl/include")
extra_compile_args.append("-I/home/mar440/WorkSpace/hfls/mpfeti_separate/")
extra_compile_args.append("-I/home/mar440/Eigen3")
extra_compile_args.append("-I/home/mar440/usr/src/boost_1_72_0")
extra_compile_args.append("-fopenmp")
extra_compile_args.append("-O3")

#######################################
extra_link_args = []
for ii in mpi_link_args:
    extra_link_args.append(ii)
extra_link_args.append("-L/home/mar440/intel/mkl/lib/intel64")
extra_link_args.append("-L/home/mar440/WorkSpace/hfls/mpfeti_separate/build/")


ext_modules_mpi =[
        Extension("pythonHdd",
                  sources=["HddApi.pyx"],
                  libraries=["hdd"],
                  language="c++",
                  extra_compile_args=extra_compile_args,
                  extra_link_args=extra_link_args
                  )
        ]


setup(
  name = "pympi",
  cmdclass = {"build_ext": build_ext},
  ext_modules = ext_modules_mpi
)
