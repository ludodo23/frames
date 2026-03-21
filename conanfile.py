import os
import re
from pathlib import Path
from conan import ConanFile
from conan.tools.files import copy
from conan.tools.cmake import cmake_layout, CMake
from conan.tools.build import check_min_cppstd


class FramesConan(ConanFile):
    name = "frames"
    settings = "os", "arch", "compiler", "build_type"
    exports_sources = "include/*", "test/*"
    no_copy_source = True
    generators = "CMakeToolchain", "CMakeDeps"
    package_type = "header-library"

    def set_version(self):
        header = Path(self.recipe_folder) / "include" / "frames.hpp"
        content = header.read_text()

        version = re.search(r'FRAMES_VERSION\s+"(.+)"', content)
        if version:
            self.version = version.group(1)

    def requirements(self):
        self.requires("Eigen/[^3.4.0]")

    def build_requirements(self):
        self.test_requires("catch2/3.4.0")

    def validate(self):
        check_min_cppstd(self, 14)

    def layout(self):
        cmake_layout(self)

    def build(self):
        if not self.conf.get("tools.build:skip_test", default=False):
            cmake = CMake(self)
            cmake.configure(build_script_folder=os.path.join(self.source_folder, "tests"))
            cmake.build()
            self.run(os.path.join(self.cpp.build.bindir, "test_interpolation"))

    def package(self):
        # This will also copy the "include" folder
        copy(self, "*.hpp", self.source_folder, self.package_folder)

    def package_info(self):
        # For header-only packages, libdirs and bindirs are not used
        # so it's necessary to set those as empty.
        self.cpp_info.bindirs = []
        self.cpp_info.libdirs = []

    def package_id(self):
        self.info.clear()