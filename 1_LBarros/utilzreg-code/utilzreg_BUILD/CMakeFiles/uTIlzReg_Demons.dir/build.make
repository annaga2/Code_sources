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
include CMakeFiles/uTIlzReg_Demons.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/uTIlzReg_Demons.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/uTIlzReg_Demons.dir/flags.make

CMakeFiles/uTIlzReg_Demons.dir/applications/uTIlzReg_Demons.cc.o: CMakeFiles/uTIlzReg_Demons.dir/flags.make
CMakeFiles/uTIlzReg_Demons.dir/applications/uTIlzReg_Demons.cc.o: ../applications/uTIlzReg_Demons.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/utilzreg_BUILD/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/uTIlzReg_Demons.dir/applications/uTIlzReg_Demons.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/uTIlzReg_Demons.dir/applications/uTIlzReg_Demons.cc.o -c /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/applications/uTIlzReg_Demons.cc

CMakeFiles/uTIlzReg_Demons.dir/applications/uTIlzReg_Demons.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/uTIlzReg_Demons.dir/applications/uTIlzReg_Demons.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/applications/uTIlzReg_Demons.cc > CMakeFiles/uTIlzReg_Demons.dir/applications/uTIlzReg_Demons.cc.i

CMakeFiles/uTIlzReg_Demons.dir/applications/uTIlzReg_Demons.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/uTIlzReg_Demons.dir/applications/uTIlzReg_Demons.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/applications/uTIlzReg_Demons.cc -o CMakeFiles/uTIlzReg_Demons.dir/applications/uTIlzReg_Demons.cc.s

CMakeFiles/uTIlzReg_Demons.dir/applications/uTIlzReg_Demons.cc.o.requires:

.PHONY : CMakeFiles/uTIlzReg_Demons.dir/applications/uTIlzReg_Demons.cc.o.requires

CMakeFiles/uTIlzReg_Demons.dir/applications/uTIlzReg_Demons.cc.o.provides: CMakeFiles/uTIlzReg_Demons.dir/applications/uTIlzReg_Demons.cc.o.requires
	$(MAKE) -f CMakeFiles/uTIlzReg_Demons.dir/build.make CMakeFiles/uTIlzReg_Demons.dir/applications/uTIlzReg_Demons.cc.o.provides.build
.PHONY : CMakeFiles/uTIlzReg_Demons.dir/applications/uTIlzReg_Demons.cc.o.provides

CMakeFiles/uTIlzReg_Demons.dir/applications/uTIlzReg_Demons.cc.o.provides.build: CMakeFiles/uTIlzReg_Demons.dir/applications/uTIlzReg_Demons.cc.o


# Object files for target uTIlzReg_Demons
uTIlzReg_Demons_OBJECTS = \
"CMakeFiles/uTIlzReg_Demons.dir/applications/uTIlzReg_Demons.cc.o"

# External object files for target uTIlzReg_Demons
uTIlzReg_Demons_EXTERNAL_OBJECTS =

uTIlzReg_Demons: CMakeFiles/uTIlzReg_Demons.dir/applications/uTIlzReg_Demons.cc.o
uTIlzReg_Demons: CMakeFiles/uTIlzReg_Demons.dir/build.make
uTIlzReg_Demons: lib/src/libDemons.a
uTIlzReg_Demons: lib/src/libDemons_plain.a
uTIlzReg_Demons: lib/src/libLDDMM.a
uTIlzReg_Demons: lib/src/libLDDMM_plain.a
uTIlzReg_Demons: lib/src/libGeoShoot.a
uTIlzReg_Demons: lib/src/libLIDM.a
uTIlzReg_Demons: lib/src/libSciCalcPack.a
uTIlzReg_Demons: CMakeFiles/uTIlzReg_Demons.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/utilzreg_BUILD/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable uTIlzReg_Demons"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/uTIlzReg_Demons.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/uTIlzReg_Demons.dir/build: uTIlzReg_Demons

.PHONY : CMakeFiles/uTIlzReg_Demons.dir/build

CMakeFiles/uTIlzReg_Demons.dir/requires: CMakeFiles/uTIlzReg_Demons.dir/applications/uTIlzReg_Demons.cc.o.requires

.PHONY : CMakeFiles/uTIlzReg_Demons.dir/requires

CMakeFiles/uTIlzReg_Demons.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/uTIlzReg_Demons.dir/cmake_clean.cmake
.PHONY : CMakeFiles/uTIlzReg_Demons.dir/clean

CMakeFiles/uTIlzReg_Demons.dir/depend:
	cd /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/utilzreg_BUILD && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/utilzreg_BUILD /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/utilzreg_BUILD /mnt/c/Users/Admin/Desktop/stage/Code_sources/1_LBarros/utilzreg-code/utilzreg_BUILD/CMakeFiles/uTIlzReg_Demons.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/uTIlzReg_Demons.dir/depend

