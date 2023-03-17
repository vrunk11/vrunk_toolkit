# Copyright 2023 vrunk-toolkit Project
#
# This file is part of vrunk-toolkit
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

include(GNUInstallDirs)

########################################################################
# Set up Windows DLL resource files
########################################################################
IF(MSVC)
    include(${CMAKE_SOURCE_DIR}/cmake/Modules/Version.cmake)

    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/vrunk11-toolkit.rc.in
        ${CMAKE_CURRENT_BINARY_DIR}/vrunk11-toolkit.rc
    @ONLY)
ENDIF(MSVC)

########################################################################
# Setup libraries used in executables
########################################################################
if(WIN32 AND NOT MINGW)
add_library(libgetopt_static STATIC
    getopt/getopt.c
)
endif()

########################################################################
# Build utility
########################################################################
add_executable(yuv2ycbcr yuv2ycbcr.c)
add_executable(yc2cvbs yc2cvbs.c)
set(INSTALL_TARGETS yuv2ycbcr)
set(INSTALL_TARGETS yc2cvbs)

target_link_libraries(yuv2ycbcr 
    ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(yc2cvbs 
    ${CMAKE_THREAD_LIBS_INIT})


if(UNIX)
target_link_libraries(yuv2ycbcr m)
target_link_libraries(yc2cvbs m)
endif()

if(WIN32 AND NOT MINGW)
target_link_libraries(yuv2ycbcr libgetopt_static)
target_link_libraries(yc2cvbs libgetopt_static)
endif()

########################################################################
# Install built library files & utilities
########################################################################
install(TARGETS ${INSTALL_TARGETS}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # .so/.dylib file
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} # .lib file
    RUNTIME DESTINATION bin              # .dll file
)
