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
include lib/src/CMakeFiles/Demons.dir/depend.make

# Include the progress variables for this target.
include lib/src/CMakeFiles/Demons.dir/progress.make

# Include the compile flags for this target's objects.
include lib/src/CMakeFiles/Demons.dir/flags.make

lib/src/CMakeFiles/Demons.dir/Demons.cc.o: lib/src/CMakeFiles/Demons.dir/flags.make
lib/src/CMakeFiles/Demons.dir/Demons.cc.o: ../lib/src/Demons.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/risser/Softwares/utilzreg-code/Build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object lib/src/CMakeFiles/Demons.dir/Demons.cc.o"
	cd /home/risser/Softwares/utilzreg-code/Build/lib/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Demons.dir/Demons.cc.o -c /home/risser/Softwares/utilzreg-code/lib/src/Demons.cc

lib/src/CMakeFiles/Demons.dir/Demons.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Demons.dir/Demons.cc.i"
	cd /home/risser/Softwares/utilzreg-code/Build/lib/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/risser/Softwares/utilzreg-code/lib/src/Demons.cc > CMakeFiles/Demons.dir/Demons.cc.i

lib/src/CMakeFiles/Demons.dir/Demons.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Demons.dir/Demons.cc.s"
	cd /home/risser/Softwares/utilzreg-code/Build/lib/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/risser/Softwares/utilzreg-code/lib/src/Demons.cc -o CMakeFiles/Demons.dir/Demons.cc.s

lib/src/CMakeFiles/Demons.dir/Demons.cc.o.requires:

.PHONY : lib/src/CMakeFiles/Demons.dir/Demons.cc.o.requires

lib/src/CMakeFiles/Demons.dir/Demons.cc.o.provides: lib/src/CMakeFiles/Demons.dir/Demons.cc.o.requires
	$(MAKE) -f lib/src/CMakeFiles/Demons.dir/build.make lib/src/CMakeFiles/Demons.dir/Demons.cc.o.provides.build
.PHONY : lib/src/CMakeFiles/Demons.dir/Demons.cc.o.provides

lib/src/CMakeFiles/Demons.dir/Demons.cc.o.provides.build: lib/src/CMakeFiles/Demons.dir/Demons.cc.o


# Object files for target Demons
Demons_OBJECTS = \
"CMakeFiles/Demons.dir/Demons.cc.o"

# External object files for target Demons
Demons_EXTERNAL_OBJECTS =

lib/src/libDemons.a: lib/src/CMakeFiles/Demons.dir/Demons.cc.o
lib/src/libDemons.a: lib/src/CMakeFiles/Demons.dir/build.make
lib/src/libDemons.a: lib/src/CMakeFiles/Demons.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/risser/Softwares/utilzreg-code/Build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libDemons.a"
	cd /home/risser/Softwares/utilzreg-code/Build/lib/src && $(CMAKE_COMMAND) -P CMakeFiles/Demons.dir/cmake_clean_target.cmake
	cd /home/risser/Softwares/utilzreg-code/Build/lib/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Demons.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/src/CMakeFiles/Demons.dir/build: lib/src/libDemons.a

.PHONY : lib/src/CMakeFiles/Demons.dir/build

lib/src/CMakeFiles/Demons.dir/requires: lib/src/CMakeFiles/Demons.dir/Demons.cc.o.requires

.PHONY : lib/src/CMakeFiles/Demons.dir/requires

lib/src/CMakeFiles/Demons.dir/clean:
	cd /home/risser/Softwares/utilzreg-code/Build/lib/src && $(CMAKE_COMMAND) -P CMakeFiles/Demons.dir/cmake_clean.cmake
.PHONY : lib/src/CMakeFiles/Demons.dir/clean

lib/src/CMakeFiles/Demons.dir/depend:
	cd /home/risser/Softwares/utilzreg-code/Build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/risser/Softwares/utilzreg-code /home/risser/Softwares/utilzreg-code/lib/src /home/risser/Softwares/utilzreg-code/Build /home/risser/Softwares/utilzreg-code/Build/lib/src /home/risser/Softwares/utilzreg-code/Build/lib/src/CMakeFiles/Demons.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/src/CMakeFiles/Demons.dir/depend
