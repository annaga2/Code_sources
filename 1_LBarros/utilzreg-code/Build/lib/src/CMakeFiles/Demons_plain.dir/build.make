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
CMAKE_SOURCE_DIR = /home/risser/Softwares/utilzreg-code

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/risser/Softwares/utilzreg-code/Build

# Include any dependencies generated for this target.
include lib/src/CMakeFiles/Demons_plain.dir/depend.make

# Include the progress variables for this target.
include lib/src/CMakeFiles/Demons_plain.dir/progress.make

# Include the compile flags for this target's objects.
include lib/src/CMakeFiles/Demons_plain.dir/flags.make

lib/src/CMakeFiles/Demons_plain.dir/Demons_plain.cc.o: lib/src/CMakeFiles/Demons_plain.dir/flags.make
lib/src/CMakeFiles/Demons_plain.dir/Demons_plain.cc.o: ../lib/src/Demons_plain.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/risser/Softwares/utilzreg-code/Build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object lib/src/CMakeFiles/Demons_plain.dir/Demons_plain.cc.o"
	cd /home/risser/Softwares/utilzreg-code/Build/lib/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Demons_plain.dir/Demons_plain.cc.o -c /home/risser/Softwares/utilzreg-code/lib/src/Demons_plain.cc

lib/src/CMakeFiles/Demons_plain.dir/Demons_plain.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Demons_plain.dir/Demons_plain.cc.i"
	cd /home/risser/Softwares/utilzreg-code/Build/lib/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/risser/Softwares/utilzreg-code/lib/src/Demons_plain.cc > CMakeFiles/Demons_plain.dir/Demons_plain.cc.i

lib/src/CMakeFiles/Demons_plain.dir/Demons_plain.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Demons_plain.dir/Demons_plain.cc.s"
	cd /home/risser/Softwares/utilzreg-code/Build/lib/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/risser/Softwares/utilzreg-code/lib/src/Demons_plain.cc -o CMakeFiles/Demons_plain.dir/Demons_plain.cc.s

lib/src/CMakeFiles/Demons_plain.dir/Demons_plain.cc.o.requires:

.PHONY : lib/src/CMakeFiles/Demons_plain.dir/Demons_plain.cc.o.requires

lib/src/CMakeFiles/Demons_plain.dir/Demons_plain.cc.o.provides: lib/src/CMakeFiles/Demons_plain.dir/Demons_plain.cc.o.requires
	$(MAKE) -f lib/src/CMakeFiles/Demons_plain.dir/build.make lib/src/CMakeFiles/Demons_plain.dir/Demons_plain.cc.o.provides.build
.PHONY : lib/src/CMakeFiles/Demons_plain.dir/Demons_plain.cc.o.provides

lib/src/CMakeFiles/Demons_plain.dir/Demons_plain.cc.o.provides.build: lib/src/CMakeFiles/Demons_plain.dir/Demons_plain.cc.o


# Object files for target Demons_plain
Demons_plain_OBJECTS = \
"CMakeFiles/Demons_plain.dir/Demons_plain.cc.o"

# External object files for target Demons_plain
Demons_plain_EXTERNAL_OBJECTS =

lib/src/libDemons_plain.a: lib/src/CMakeFiles/Demons_plain.dir/Demons_plain.cc.o
lib/src/libDemons_plain.a: lib/src/CMakeFiles/Demons_plain.dir/build.make
lib/src/libDemons_plain.a: lib/src/CMakeFiles/Demons_plain.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/risser/Softwares/utilzreg-code/Build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libDemons_plain.a"
	cd /home/risser/Softwares/utilzreg-code/Build/lib/src && $(CMAKE_COMMAND) -P CMakeFiles/Demons_plain.dir/cmake_clean_target.cmake
	cd /home/risser/Softwares/utilzreg-code/Build/lib/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Demons_plain.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/src/CMakeFiles/Demons_plain.dir/build: lib/src/libDemons_plain.a

.PHONY : lib/src/CMakeFiles/Demons_plain.dir/build

lib/src/CMakeFiles/Demons_plain.dir/requires: lib/src/CMakeFiles/Demons_plain.dir/Demons_plain.cc.o.requires

.PHONY : lib/src/CMakeFiles/Demons_plain.dir/requires

lib/src/CMakeFiles/Demons_plain.dir/clean:
	cd /home/risser/Softwares/utilzreg-code/Build/lib/src && $(CMAKE_COMMAND) -P CMakeFiles/Demons_plain.dir/cmake_clean.cmake
.PHONY : lib/src/CMakeFiles/Demons_plain.dir/clean

lib/src/CMakeFiles/Demons_plain.dir/depend:
	cd /home/risser/Softwares/utilzreg-code/Build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/risser/Softwares/utilzreg-code /home/risser/Softwares/utilzreg-code/lib/src /home/risser/Softwares/utilzreg-code/Build /home/risser/Softwares/utilzreg-code/Build/lib/src /home/risser/Softwares/utilzreg-code/Build/lib/src/CMakeFiles/Demons_plain.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/src/CMakeFiles/Demons_plain.dir/depend

