# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = /home/mihail/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/192.6817.32/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/mihail/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/192.6817.32/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/mihail/Documents/workspace/02_TSISA/TSISA_lab_06

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/mihail/Documents/workspace/02_TSISA/TSISA_lab_06/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/TSISA_lab_06.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/TSISA_lab_06.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/TSISA_lab_06.dir/flags.make

CMakeFiles/TSISA_lab_06.dir/main.cpp.o: CMakeFiles/TSISA_lab_06.dir/flags.make
CMakeFiles/TSISA_lab_06.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mihail/Documents/workspace/02_TSISA/TSISA_lab_06/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/TSISA_lab_06.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/TSISA_lab_06.dir/main.cpp.o -c /home/mihail/Documents/workspace/02_TSISA/TSISA_lab_06/main.cpp

CMakeFiles/TSISA_lab_06.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TSISA_lab_06.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mihail/Documents/workspace/02_TSISA/TSISA_lab_06/main.cpp > CMakeFiles/TSISA_lab_06.dir/main.cpp.i

CMakeFiles/TSISA_lab_06.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TSISA_lab_06.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mihail/Documents/workspace/02_TSISA/TSISA_lab_06/main.cpp -o CMakeFiles/TSISA_lab_06.dir/main.cpp.s

# Object files for target TSISA_lab_06
TSISA_lab_06_OBJECTS = \
"CMakeFiles/TSISA_lab_06.dir/main.cpp.o"

# External object files for target TSISA_lab_06
TSISA_lab_06_EXTERNAL_OBJECTS =

TSISA_lab_06: CMakeFiles/TSISA_lab_06.dir/main.cpp.o
TSISA_lab_06: CMakeFiles/TSISA_lab_06.dir/build.make
TSISA_lab_06: CMakeFiles/TSISA_lab_06.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/mihail/Documents/workspace/02_TSISA/TSISA_lab_06/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable TSISA_lab_06"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/TSISA_lab_06.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/TSISA_lab_06.dir/build: TSISA_lab_06

.PHONY : CMakeFiles/TSISA_lab_06.dir/build

CMakeFiles/TSISA_lab_06.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/TSISA_lab_06.dir/cmake_clean.cmake
.PHONY : CMakeFiles/TSISA_lab_06.dir/clean

CMakeFiles/TSISA_lab_06.dir/depend:
	cd /home/mihail/Documents/workspace/02_TSISA/TSISA_lab_06/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mihail/Documents/workspace/02_TSISA/TSISA_lab_06 /home/mihail/Documents/workspace/02_TSISA/TSISA_lab_06 /home/mihail/Documents/workspace/02_TSISA/TSISA_lab_06/cmake-build-debug /home/mihail/Documents/workspace/02_TSISA/TSISA_lab_06/cmake-build-debug /home/mihail/Documents/workspace/02_TSISA/TSISA_lab_06/cmake-build-debug/CMakeFiles/TSISA_lab_06.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/TSISA_lab_06.dir/depend

