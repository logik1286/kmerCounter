cmake_minimum_required (VERSION 3.2)
project (kmerCounter CXX)

set(CMAKE_BUILD_TYPE Release)

if (MSVC)
        message( FATAL_ERROR "Windows is not supported." )
endif()


if (CMAKE_COMPILER_IS_GNUCXX)

	if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8.4")
		message(FATAL_ERROR	"g++ version is ${CMAKE_CXX_COMPILER_VERSION}. Package was tested with 4.8.4. Please upgrade to at least 4.8.4")
	endif()

endif()

set(SRC_FILES
	kmerError.cpp
	fastqReader.cpp
	kmerUtils.cpp
	merizer.cpp
	main.cpp
	)

set(INC_FILES
	kmerError.h
	fastqReader.h
	kmerUtils.h
	kmerCounter.h
	merizer.h
	keyCache.h
	)

add_executable(kmerCounter ${SRC_FILES} ${INC_FILES})
target_link_libraries(kmerCounter kyotocabinet)
target_compile_features(kmerCounter PUBLIC cxx_constexpr)
