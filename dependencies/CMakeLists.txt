cmake_minimum_required(VERSION 3.12)

project(MitsubaDependencies)

set(CMAKE_CXX_STANDARD 17)

# Following lines are configurations need to modify.
include(ProcessorCount)
ProcessorCount(COMPILE_THREAD)
message("Thread Count: ${COMPILE_THREAD}")

set(ENV{http_proxy}  "http://127.0.0.1:10809") # for http
set(ENV{https_proxy} "http://127.0.0.1:10809") # for https

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
		"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
		FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(DEPENDENCIES_DOWNLOAD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Downloads")
set(DEPENDENCIES_BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

set(BUILD_COMMAND_OPTS --target install --config ${CMAKE_BUILD_TYPE})
file(MAKE_DIRECTORY ${DEPENDENCIES_DOWNLOAD_DIR})
file(MAKE_DIRECTORY ${DEPENDENCIES_BINARY_DIR}/include)
file(MAKE_DIRECTORY ${DEPENDENCIES_BINARY_DIR}/lib)
file(MAKE_DIRECTORY ${DEPENDENCIES_BINARY_DIR}/bin)

function(DownloadExternal name url target md5)
	set(file_download_path "${DEPENDENCIES_DOWNLOAD_DIR}/${target}")
	if(NOT EXISTS ${file_download_path})
		message(STATUS "Downloading ${name}...")
		file(DOWNLOAD ${url} ${file_download_path} SHOW_PROGRESS EXPECTED_MD5 ${md5})
	endif()
	message(STATUS "${name} downloaded.")
endfunction(DownloadExternal)

function(CloneExternal name url)
	if(NOT EXISTS "${DEPENDENCIES_DOWNLOAD_DIR}/${name}")
		message(STATUS "Cloning ${name}...")
		execute_process(COMMAND "git" "clone" "--recurse-submodules" ${url} WORKING_DIRECTORY ${DEPENDENCIES_DOWNLOAD_DIR})
	endif()
	message(STATUS "${name} cloned")
endfunction()

# ########
# # zlib #
# ########
set(ZLIB_PROJECT_DIR "${DEPENDENCIES_DOWNLOAD_DIR}/zlib")
# file(MAKE_DIRECTORY ${ZLIB_PROJECT_DIR})
CloneExternal("zlib" "https://github.com/madler/zlib.git")
execute_process(COMMAND ${CMAKE_COMMAND} 
	-DCMAKE_INSTALL_PREFIX=${DEPENDENCIES_BINARY_DIR}
	-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
	-G ${CMAKE_GENERATOR}
	-DCMAKE_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}
	-DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
	${ZLIB_PROJECT_DIR}
	WORKING_DIRECTORY ${ZLIB_PROJECT_DIR}
)
execute_process(COMMAND ${CMAKE_COMMAND} --build ${ZLIB_PROJECT_DIR} ${BUILD_COMMAND_OPTS} -j${COMPILE_THREAD})

#############################################################################################

# ###########
# # OpenEXR #
# ###########
set(OPENEXR_PROJECT_DIR "${DEPENDENCIES_DOWNLOAD_DIR}/OpenEXR")
# file(MAKE_DIRECTORY ${OPENEXR_PROJECT_DIR})
CloneExternal("OpenEXR" "https://github.com/AcademySoftwareFoundation/openexr.git")
execute_process(COMMAND ${CMAKE_COMMAND} 
	-DCMAKE_INSTALL_PREFIX=${DEPENDENCIES_BINARY_DIR}
	-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
	-G ${CMAKE_GENERATOR}
	-DCMAKE_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}
	-DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
	${OPENEXR_PROJECT_DIR}
	WORKING_DIRECTORY ${OPENEXR_PROJECT_DIR}
)
execute_process(COMMAND ${CMAKE_COMMAND} --build ${OPENEXR_PROJECT_DIR} ${BUILD_COMMAND_OPTS} -j${COMPILE_THREAD})

# ##########################
# # libjpeg(libjpeg-turbo) #
# ##########################
set(LIBJPEG_PROJECT_DIR "${DEPENDENCIES_DOWNLOAD_DIR}/libjpeg-turbo")
# file(MAKE_DIRECTORY ${LIBJPEG_PROJECT_DIR})
CloneExternal("libjpeg" "https://github.com/libjpeg-turbo/libjpeg-turbo.git")
execute_process(COMMAND ${CMAKE_COMMAND} 
	-DCMAKE_INSTALL_PREFIX=${DEPENDENCIES_BINARY_DIR}
	-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
	-G ${CMAKE_GENERATOR}
	-DCMAKE_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}
	-DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
	${LIBJPEG_PROJECT_DIR}
	WORKING_DIRECTORY ${LIBJPEG_PROJECT_DIR}
)
execute_process(COMMAND ${CMAKE_COMMAND} --build ${LIBJPEG_PROJECT_DIR} ${BUILD_COMMAND_OPTS} -j${COMPILE_THREAD})

##########
# libpng #
##########
set(LIBPNG_PROJECT_DIR "${DEPENDENCIES_DOWNLOAD_DIR}/libpng")
# file(MAKE_DIRECTORY ${LIBPNG_PROJECT_DIR})
CloneExternal("libpng" "https://github.com/glennrp/libpng.git")
execute_process(COMMAND ${CMAKE_COMMAND} 
	-DCMAKE_INSTALL_PREFIX=${DEPENDENCIES_BINARY_DIR}
	-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
	-G ${CMAKE_GENERATOR}
	-DZLIB_ROOT=${DEPENDENCIES_BINARY_DIR}
	-DCMAKE_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}
	-DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
	${LIBPNG_PROJECT_DIR}
	WORKING_DIRECTORY ${LIBPNG_PROJECT_DIR}
)
execute_process(COMMAND ${CMAKE_COMMAND} --build ${LIBPNG_PROJECT_DIR} ${BUILD_COMMAND_OPTS} -j${COMPILE_THREAD})

# #########
# # boost #
# #########
set(BOOST_PROJECT_DIR "${DEPENDENCIES_DOWNLOAD_DIR}/boost")
if(NOT EXISTS ${BOOST_PROJECT_DIR})
	message(STATUS "Cloning boost...")
	execute_process(COMMAND "git" "clone" "--recurse-submodules" "--depth=1" "--branch=boost-1.77.0" "https://github.com/boostorg/boost.git" WORKING_DIRECTORY ${DEPENDENCIES_DOWNLOAD_DIR})
endif()
message(STATUS "boost cloned")
execute_process(COMMAND ${CMAKE_COMMAND} 
	-DCMAKE_INSTALL_PREFIX=${DEPENDENCIES_BINARY_DIR}
	-DCMAKE_BUILD_TYPE=Debug
	-G ${CMAKE_GENERATOR}
	-DZLIB_ROOT=${DEPENDENCIES_BINARY_DIR}
	-DCMAKE_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}
	-DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
	-DBUILD_SHARED_LIBS=on
	-DBOOST_ENABLE_PYTHON=on
	-DBOOST_INSTALL_INCLUDE_SUBDIR=/
	-DBOOST_LOCALE_ENABLE_ICU=off
	${BOOST_PROJECT_DIR}
	WORKING_DIRECTORY ${BOOST_PROJECT_DIR}
)
execute_process(COMMAND ${CMAKE_COMMAND} --build ${BOOST_PROJECT_DIR} --target install --config Debug -j${COMPILE_THREAD})
execute_process(COMMAND ${CMAKE_COMMAND} 
	-DCMAKE_INSTALL_PREFIX=${DEPENDENCIES_BINARY_DIR}
	-DCMAKE_BUILD_TYPE=Release
	-G ${CMAKE_GENERATOR}
	-DZLIB_ROOT=${DEPENDENCIES_BINARY_DIR}
	-DCMAKE_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}
	-DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
	-DBUILD_SHARED_LIBS=on
	-DBOOST_ENABLE_PYTHON=on
	-DBOOST_INSTALL_INCLUDE_SUBDIR=/
	-DBOOST_LOCALE_ENABLE_ICU=off
	${BOOST_PROJECT_DIR}
	WORKING_DIRECTORY ${BOOST_PROJECT_DIR}
)
execute_process(COMMAND ${CMAKE_COMMAND} --build ${BOOST_PROJECT_DIR} --target install --config Release -j${COMPILE_THREAD})

###########
# xerces-c #
###########
set(XERCES_PROJECT_DIR "${DEPENDENCIES_DOWNLOAD_DIR}/xerces-c")
CloneExternal("xerces" "https://github.com/apache/xerces-c.git")
execute_process(COMMAND ${CMAKE_COMMAND} 
	-DCMAKE_INSTALL_PREFIX=${DEPENDENCIES_BINARY_DIR}
	-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
	-G ${CMAKE_GENERATOR}
	-DCMAKE_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}
	-DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
	-Dtranscoder=windows
	${XERCES_PROJECT_DIR}
	WORKING_DIRECTORY ${XERCES_PROJECT_DIR}
)
execute_process(COMMAND ${CMAKE_COMMAND} --build ${XERCES_PROJECT_DIR} ${BUILD_COMMAND_OPTS} -j${COMPILE_THREAD})

########
# glew #
########
set(GLEW_PROJECT_DIR "${DEPENDENCIES_DOWNLOAD_DIR}/glew-2.2.0")
file(MAKE_DIRECTORY ${GLEW_PROJECT_DIR})
DownloadExternal("glew-2.2.0" "https://github.com/nigels-com/glew/releases/download/glew-2.2.0/glew-2.2.0-win32.zip" "glew-2.2.0-win32.zip" "1FEDDFE8696C192FA46A0DF8EAC7D4BF")
execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf glew-2.2.0-win32.zip WORKING_DIRECTORY ${DEPENDENCIES_DOWNLOAD_DIR})
file(GLOB GLEW_INC "${GLEW_PROJECT_DIR}/include/*")
file(GLOB GLEW_LIB "${GLEW_PROJECT_DIR}/lib/Release/x64/*")
file(GLOB GLEW_BIN "${GLEW_PROJECT_DIR}/bin/Release/x64/*")
file(COPY ${GLEW_INC} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/include)
file(COPY ${GLEW_LIB} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib)
file(COPY ${GLEW_BIN} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin)

#########
# Eigen #
#########
set(EIGEN_PROJECT_DIR "${DEPENDENCIES_DOWNLOAD_DIR}/eigen-3.3.9")
set(EIGEN_INCLUDE_DIR "${DEPENDENCIES_BINARY_DIR}/include/Eigen")
DownloadExternal("Eigen" "https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.tar.gz" "eigen-3.3.9.tar.gz" "609286804B0F79BE622CCF7F9FF2B660")
execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf eigen-3.3.9.tar.gz WORKING_DIRECTORY ${DEPENDENCIES_DOWNLOAD_DIR})
file(MAKE_DIRECTORY ${EIGEN_INCLUDE_DIR})
file(GLOB EIGEN_INC "${EIGEN_PROJECT_DIR}/Eigen/*")
list(REMOVE_ITEM EIGEN_INC "${EIGEN_PROJECT_DIR}/Eigen/CMakeLists.txt")
file(COPY ${EIGEN_INC} DESTINATION ${EIGEN_INCLUDE_DIR})

# ##########
# # OpenGL #
# ##########
DownloadExternal("glext" "https://www.khronos.org/registry/OpenGL/api/GL/glext.h" "glext.h" "3CF5021BC0371C754E5CD9157E431AAE")
DownloadExternal("khr" "https://www.khronos.org/registry/EGL/api/KHR/khrplatform.h" "khrplatform.h" "D03191518AC2CFC3C10D22DF034B154A")
file(MAKE_DIRECTORY ${DEPENDENCIES_BINARY_DIR}/include/KHR)
file(COPY ${DEPENDENCIES_DOWNLOAD_DIR}/khrplatform.h DESTINATION ${DEPENDENCIES_BINARY_DIR}/include/KHR/)
file(COPY ${DEPENDENCIES_DOWNLOAD_DIR}/glext.h DESTINATION ${DEPENDENCIES_BINARY_DIR}/include/GL/)

########
# FFTW #
########
set(FFTW_PROJECT_DIR "${DEPENDENCIES_DOWNLOAD_DIR}/fftw-3.3.5-dll64")
file(MAKE_DIRECTORY ${FFTW_PROJECT_DIR})
DownloadExternal("fftw" "https://fftw.org/pub/fftw/fftw-3.3.5-dll64.zip" "fftw-3.3.5-dll64.zip" "CB3C5AD19A89864F036E7A2DD5BE168C")
execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf ../fftw-3.3.5-dll64.zip WORKING_DIRECTORY ${FFTW_PROJECT_DIR})
execute_process(COMMAND lib /def:libfftw3f-3.def WORKING_DIRECTORY ${FFTW_PROJECT_DIR})
execute_process(COMMAND lib /def:libfftw3-3.def WORKING_DIRECTORY ${FFTW_PROJECT_DIR})
execute_process(COMMAND lib /def:libfftw3l-3.def WORKING_DIRECTORY ${FFTW_PROJECT_DIR})
file(GLOB FFTW_INC "${FFTW_PROJECT_DIR}/*.h")
file(GLOB FFTW_LIB "${FFTW_PROJECT_DIR}/*.lib")
file(GLOB FFTW_BIN "${FFTW_PROJECT_DIR}/*.dll")
file(COPY ${FFTW_INC} DESTINATION ${DEPENDENCIES_BINARY_DIR}/include)
file(COPY ${FFTW_LIB} DESTINATION ${DEPENDENCIES_BINARY_DIR}/lib)
file(COPY ${FFTW_BIN} DESTINATION ${DEPENDENCIES_BINARY_DIR}/bin)
