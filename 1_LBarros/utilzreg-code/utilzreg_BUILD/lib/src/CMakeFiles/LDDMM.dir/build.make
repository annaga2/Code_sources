# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/utilzreg_BUILD

# Include any dependencies generated for this target.
include lib/src/CMakeFiles/LDDMM.dir/depend.make

# Include the progress variables for this target.
include lib/src/CMakeFiles/LDDMM.dir/progress.make

# Include the compile flags for this target's objects.
include lib/src/CMakeFiles/LDDMM.dir/flags.make

lib/src/CMakeFiles/LDDMM.dir/LDDMM.cc.o: lib/src/CMakeFiles/LDDMM.dir/flags.make
lib/src/CMakeFiles/LDDMM.dir/LDDMM.cc.o: ../lib/src/LDDMM.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/utilzreg_BUILD/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object lib/src/CMakeFiles/LDDMM.dir/LDDMM.cc.o"
	cd /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/utilzreg_BUILD/lib/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/LDDMM.dir/LDDMM.cc.o -c /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/lib/src/LDDMM.cc

lib/src/CMakeFiles/LDDMM.dir/LDDMM.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LDDMM.dir/LDDMM.cc.i"
	cd /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/utilzreg_BUILD/lib/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/lib/src/LDDMM.cc > CMakeFiles/LDDMM.dir/LDDMM.cc.i

lib/src/CMakeFiles/LDDMM.dir/LDDMM.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LDDMM.dir/LDDMM.cc.s"
	cd /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/utilzreg_BUILD/lib/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/lib/src/LDDMM.cc -o CMakeFiles/LDDMM.dir/LDDMM.cc.s

lib/src/CMakeFiles/LDDMM.dir/LDDMM.cc.o.requires:

.PHONY : lib/src/CMakeFiles/LDDMM.dir/LDDMM.cc.o.requires

lib/src/CMakeFiles/LDDMM.dir/LDDMM.cc.o.provides: lib/src/CMakeFiles/LDDMM.dir/LDDMM.cc.o.requires
	$(MAKE) -f lib/src/CMakeFiles/LDDMM.dir/build.make lib/src/CMakeFiles/LDDMM.dir/LDDMM.cc.o.provides.build
.PHONY : lib/src/CMakeFiles/LDDMM.dir/LDDMM.cc.o.provides

lib/src/CMakeFiles/LDDMM.dir/LDDMM.cc.o.provides.build: lib/src/CMakeFiles/LDDMM.dir/LDDMM.cc.o


# Object files for target LDDMM
LDDMM_OBJECTS = \
"CMakeFiles/LDDMM.dir/LDDMM.cc.o"

# External object files for target LDDMM
LDDMM_EXTERNAL_OBJECTS =

lib/src/libLDDMM.a: lib/src/CMakeFiles/LDDMM.dir/LDDMM.cc.o
lib/src/libLDDMM.a: lib/src/CMakeFiles/LDDMM.dir/build.make
lib/src/libLDDMM.a: lib/src/CMakeFiles/LDDMM.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/utilzreg_BUILD/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libLDDMM.a"
	cd /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/utilzreg_BUILD/lib/src && $(CMAKE_COMMAND) -P CMakeFiles/LDDMM.dir/cmake_clean_target.cmake
	cd /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/utilzreg_BUILD/lib/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/LDDMM.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/src/CMakeFiles/LDDMM.dir/build: lib/src/libLDDMM.a

.PHONY : lib/src/CMakeFiles/LDDMM.dir/build

lib/src/CMakeFiles/LDDMM.dir/requires: lib/src/CMakeFiles/LDDMM.dir/LDDMM.cc.o.requires

.PHONY : lib/src/CMakeFiles/LDDMM.dir/requires

lib/src/CMakeFiles/LDDMM.dir/clean:
	cd /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/utilzreg_BUILD/lib/src && $(CMAKE_COMMAND) -P CMakeFiles/LDDMM.dir/cmake_clean.cmake
.PHONY : lib/src/CMakeFiles/LDDMM.dir/clean

lib/src/CMakeFiles/LDDMM.dir/depend:
	cd /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/utilzreg_BUILD && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/lib/src /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/utilzreg_BUILD /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/utilzreg_BUILD/lib/src /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/utilzreg_BUILD/lib/src/CMakeFiles/LDDMM.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/src/CMakeFiles/LDDMM.dir/depend
