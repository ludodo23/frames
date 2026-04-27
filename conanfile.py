import os
import re
from pathlib import Path
from conan import ConanFile
from conan.tools.files import copy
from conan.tools.cmake import cmake_layout, CMake, CMakeToolchain, CMakeDeps
from conan.tools.build import check_min_cppstd


class FramesConan(ConanFile):
    name = "frames"
    settings = "os", "arch", "compiler", "build_type"
    exports_sources = "include/*", "tests/*"
    no_copy_source = True
    package_type = "header-library"
    options = {
        "with_eigen": [True, False],
    }
    default_options = {
        "with_eigen": True,
    }

    def set_version(self):
        header = Path(self.recipe_folder) / "include" / "frames.hpp"
        content = header.read_text()
        version = re.search(r'FRAMES_VERSION\s+"(.+)"', content)
        if version:
            self.version = version.group(1)

    def requirements(self):
        self.requires("interpolation/0.2.0")
        if self.options.with_eigen:
            self.requires("eigen/[^3.4.0]")

    def build_requirements(self):
        self.test_requires("catch2/3.4.0")

    def validate(self):
        check_min_cppstd(self, 14)

    def layout(self):
        cmake_layout(self)

    def generate(self):
        tc = CMakeToolchain(self)
        # Transmet l'option with_eigen au CMakeLists.txt
        tc.variables["FRAMES_WITH_EIGEN"] = self.options.with_eigen
        tc.generate()
        deps = CMakeDeps(self)
        deps.generate()

    def build(self):
        if self.conf.get("tools.build:skip_test", default=False):
            return
        cmake = CMake(self)
        cmake.configure()   # prend la racine, plus build_script_folder
        cmake.build()
        cmake.ctest(cli_args=["--output-on-failure"])

    def package(self):
        copy(self, "*.hpp", self.source_folder, self.package_folder)

    def package_info(self):
        self.cpp_info.bindirs = []
        self.cpp_info.libdirs = []
        if self.options.with_eigen:
            self.cpp_info.defines = ["FRAMES_WITH_EIGEN"]

    def package_id(self):
        self.info.clear()