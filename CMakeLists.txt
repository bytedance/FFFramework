cmake_minimum_required(VERSION 3.10)

set(PROJECT_NAME FFFramework)
if(NOT ffmpeg_include)
    set(ffmpeg_include "..\\Dependencies\\win32\\32\\include")
endif()
if(NOT ffmpeg_avutil_lib)
    set(ffmpeg_avutil_lib "..\\..\\Dependencies\\win32\\32\\lib\\avutil.lib")
endif()
#add the executable

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
file(GLOB SOURCE_SET
     "*.h"
     "*.cpp"
)

add_library(${PROJECT_NAME} SHARED ${SOURCE_SET})

target_include_directories(${PROJECT_NAME}
PUBLIC .
${ffmpeg_include}
)

if(MSVC)
    if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "ARM")
        target_compile_options(${PROJECT_NAME} PRIVATE
            "$<$<CONFIG:Debug>:/MDd>"
            "$<$<CONFIG:Release>:/MD>"
        )
    elseif("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "Win32")
        target_compile_options(${PROJECT_NAME} PRIVATE
            "$<$<CONFIG:Debug>:/MDd>"
            "$<$<CONFIG:Release>:/MD>"
        )
    elseif("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
        target_compile_options(${PROJECT_NAME} PRIVATE
            "$<$<CONFIG:Debug>:/MDd>"
            "$<$<CONFIG:Release>:/MD>"
        )
    endif()

    target_link_options(${PROJECT_NAME} PRIVATE
                $<$<CONFIG:Debug>:
                    /INCREMENTAL;
                >
                $<$<CONFIG:Release>:
                    /OPT:REF;
                    /OPT:ICF;
                    /INCREMENTAL:NO;
                >
            )
endif()
message("default ffmpeg_include ${ffmpeg_include}")
message("default ffmpeg_avutil_lib ${ffmpeg_avutil_lib}")


set(ADDITIONAL_LIBRARY_DEPENDENCIES
    ${ffmpeg_avutil_lib}
)


target_link_libraries(${PROJECT_NAME} PUBLIC "${ADDITIONAL_LIBRARY_DEPENDENCIES}")

if(MSVC)
  add_definitions(/MP)
endif()