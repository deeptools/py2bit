from setuptools import setup, Extension
import glob
import sysconfig

srcs = [x for x in 
    glob.glob("lib2bit/*.c")]
srcs.append("py2bit.c")

additional_libs = [sysconfig.get_config_var("LIBDIR"), sysconfig.get_config_var("LIBPL")]

module1 = Extension('py2bit',
                    sources = srcs,
                    library_dirs = additional_libs, 
                    include_dirs = ['lib2bit', sysconfig.get_config_var("INCLUDEPY")])

setup_args = dict(ext_modules=[module1])
setup(**setup_args)
