import sys
from setuptools import setup, Extension
from Cython.Build import cythonize

if sys.platform == "win32":
    compile_args = [
        "/std:c++17", 
        "/utf-8", 
        "/O2",         
        "/Ot",         
        "/Ob2",        
        "/fp:fast",    
        "/arch:AVX2",  
        "/GL",         
    ]
    link_args = ["/LTCG"]
else:
    #para Linux/Mac (GCC/Clang)
    compile_args = ["-std=c++17", "-O3", "-Ofast", "-march=native", "-flto"]
    link_args = ["-flto"]

ext_modules = [
    Extension(
        "gsim",          # Nome do módulo final
        sources=["gsim_corrigido (2).pyx"], 
        language="c++",
        extra_compile_args=compile_args,
        extra_link_args=link_args,
    )
]

setup(
    name="_gsim_backend",
    ext_modules=cythonize(
        ext_modules, 
        annotate=True,       #pra gerar o html
        compiler_directives={
            'language_level': "3",
            'boundscheck': False,
            'wraparound': False,
            'cdivision': True
        }
    ),
)