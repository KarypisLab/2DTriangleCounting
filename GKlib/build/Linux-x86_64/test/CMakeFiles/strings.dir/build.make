# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Produce verbose output by default.
VERBOSE = 1

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
CMAKE_COMMAND = /panfs/roc/msisoft/cmake/3.10.2/bin/cmake

# The command to remove a file.
RM = /panfs/roc/msisoft/cmake/3.10.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/build/Linux-x86_64

# Include any dependencies generated for this target.
include test/CMakeFiles/strings.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/strings.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/strings.dir/flags.make

test/CMakeFiles/strings.dir/strings.c.o: test/CMakeFiles/strings.dir/flags.make
test/CMakeFiles/strings.dir/strings.c.o: ../../test/strings.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/build/Linux-x86_64/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object test/CMakeFiles/strings.dir/strings.c.o"
	cd /home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/build/Linux-x86_64/test && /panfs/roc/msisoft/gcc/8.2.0/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/strings.dir/strings.c.o   -c /home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/test/strings.c

test/CMakeFiles/strings.dir/strings.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/strings.dir/strings.c.i"
	cd /home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/build/Linux-x86_64/test && /panfs/roc/msisoft/gcc/8.2.0/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/test/strings.c > CMakeFiles/strings.dir/strings.c.i

test/CMakeFiles/strings.dir/strings.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/strings.dir/strings.c.s"
	cd /home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/build/Linux-x86_64/test && /panfs/roc/msisoft/gcc/8.2.0/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/test/strings.c -o CMakeFiles/strings.dir/strings.c.s

test/CMakeFiles/strings.dir/strings.c.o.requires:

.PHONY : test/CMakeFiles/strings.dir/strings.c.o.requires

test/CMakeFiles/strings.dir/strings.c.o.provides: test/CMakeFiles/strings.dir/strings.c.o.requires
	$(MAKE) -f test/CMakeFiles/strings.dir/build.make test/CMakeFiles/strings.dir/strings.c.o.provides.build
.PHONY : test/CMakeFiles/strings.dir/strings.c.o.provides

test/CMakeFiles/strings.dir/strings.c.o.provides.build: test/CMakeFiles/strings.dir/strings.c.o


# Object files for target strings
strings_OBJECTS = \
"CMakeFiles/strings.dir/strings.c.o"

# External object files for target strings
strings_EXTERNAL_OBJECTS =

test/strings: test/CMakeFiles/strings.dir/strings.c.o
test/strings: test/CMakeFiles/strings.dir/build.make
test/strings: libGKlib.a
test/strings: test/CMakeFiles/strings.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/build/Linux-x86_64/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable strings"
	cd /home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/build/Linux-x86_64/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/strings.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/strings.dir/build: test/strings

.PHONY : test/CMakeFiles/strings.dir/build

test/CMakeFiles/strings.dir/requires: test/CMakeFiles/strings.dir/strings.c.o.requires

.PHONY : test/CMakeFiles/strings.dir/requires

test/CMakeFiles/strings.dir/clean:
	cd /home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/build/Linux-x86_64/test && $(CMAKE_COMMAND) -P CMakeFiles/strings.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/strings.dir/clean

test/CMakeFiles/strings.dir/depend:
	cd /home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/build/Linux-x86_64 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib /home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/test /home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/build/Linux-x86_64 /home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/build/Linux-x86_64/test /home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/build/Linux-x86_64/test/CMakeFiles/strings.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/strings.dir/depend
