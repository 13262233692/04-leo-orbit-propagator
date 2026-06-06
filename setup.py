import os
import sys
import subprocess
import shutil
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        super().__init__(name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        extdir = os.path.join(extdir, "leo_propagator")

        cfg = "Debug" if self.debug else "Release"

        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",
        ]

        build_args = ["--config", cfg]

        if sys.platform.startswith("win"):
            cmake_args += [f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}"]
            build_args += ["--", "/m"]
        else:
            build_args += ["--", "-j4"]

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        subprocess.check_call(
            ["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp
        )
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args, cwd=self.build_temp
        )

        src_file = os.path.join(extdir, "_leo_propagator" + self._get_so_suffix())
        if not os.path.exists(src_file):
            for f in os.listdir(extdir):
                if f.startswith("_leo_propagator") and f.endswith(self._get_so_suffix()):
                    src_file = os.path.join(extdir, f)
                    break

    def _get_so_suffix(self):
        if sys.platform.startswith("win"):
            return ".pyd"
        elif sys.platform == "darwin":
            return ".so"
        else:
            return ".so"


with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

with open("requirements.txt", "r", encoding="utf-8") as f:
    requirements = [line.strip() for line in f if line.strip() and not line.startswith("#")]

setup(
    name="leo-orbit-propagator",
    version="1.0.0",
    author="LEO Satellite Team",
    description="High precision LEO satellite orbit propagation engine",
    long_description=long_description,
    long_description_content_type="text/markdown",
    package_dir={"": "python"},
    packages=find_packages(where="python"),
    ext_modules=[CMakeExtension("_leo_propagator")],
    cmdclass=dict(build_ext=CMakeBuild),
    install_requires=requirements,
    python_requires=">=3.8",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Programming Language :: Python :: 3",
        "Programming Language :: C++",
    ],
)
