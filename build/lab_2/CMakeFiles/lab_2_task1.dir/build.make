# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /workspaces/trs_labs

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /workspaces/trs_labs/build

# Include any dependencies generated for this target.
include lab_2/CMakeFiles/lab_2_task1.dir/depend.make

# Include the progress variables for this target.
include lab_2/CMakeFiles/lab_2_task1.dir/progress.make

# Include the compile flags for this target's objects.
include lab_2/CMakeFiles/lab_2_task1.dir/flags.make

lab_2/CMakeFiles/lab_2_task1.dir/task_1.cpp.o: lab_2/CMakeFiles/lab_2_task1.dir/flags.make
lab_2/CMakeFiles/lab_2_task1.dir/task_1.cpp.o: ../lab_2/task_1.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/workspaces/trs_labs/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object lab_2/CMakeFiles/lab_2_task1.dir/task_1.cpp.o"
	cd /workspaces/trs_labs/build/lab_2 && /usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lab_2_task1.dir/task_1.cpp.o -c /workspaces/trs_labs/lab_2/task_1.cpp

lab_2/CMakeFiles/lab_2_task1.dir/task_1.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lab_2_task1.dir/task_1.cpp.i"
	cd /workspaces/trs_labs/build/lab_2 && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /workspaces/trs_labs/lab_2/task_1.cpp > CMakeFiles/lab_2_task1.dir/task_1.cpp.i

lab_2/CMakeFiles/lab_2_task1.dir/task_1.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lab_2_task1.dir/task_1.cpp.s"
	cd /workspaces/trs_labs/build/lab_2 && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /workspaces/trs_labs/lab_2/task_1.cpp -o CMakeFiles/lab_2_task1.dir/task_1.cpp.s

# Object files for target lab_2_task1
lab_2_task1_OBJECTS = \
"CMakeFiles/lab_2_task1.dir/task_1.cpp.o"

# External object files for target lab_2_task1
lab_2_task1_EXTERNAL_OBJECTS =

lab_2/lab_2_task1: lab_2/CMakeFiles/lab_2_task1.dir/task_1.cpp.o
lab_2/lab_2_task1: lab_2/CMakeFiles/lab_2_task1.dir/build.make
lab_2/lab_2_task1: lab_2/CMakeFiles/lab_2_task1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/workspaces/trs_labs/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable lab_2_task1"
	cd /workspaces/trs_labs/build/lab_2 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lab_2_task1.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lab_2/CMakeFiles/lab_2_task1.dir/build: lab_2/lab_2_task1

.PHONY : lab_2/CMakeFiles/lab_2_task1.dir/build

lab_2/CMakeFiles/lab_2_task1.dir/clean:
	cd /workspaces/trs_labs/build/lab_2 && $(CMAKE_COMMAND) -P CMakeFiles/lab_2_task1.dir/cmake_clean.cmake
.PHONY : lab_2/CMakeFiles/lab_2_task1.dir/clean

lab_2/CMakeFiles/lab_2_task1.dir/depend:
	cd /workspaces/trs_labs/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /workspaces/trs_labs /workspaces/trs_labs/lab_2 /workspaces/trs_labs/build /workspaces/trs_labs/build/lab_2 /workspaces/trs_labs/build/lab_2/CMakeFiles/lab_2_task1.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lab_2/CMakeFiles/lab_2_task1.dir/depend
