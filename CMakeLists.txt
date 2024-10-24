cmake_minimum_required(VERSION 3.28)

project(cfb)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)

if (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif()


file(GLOB SRCS "${PROJECT_SOURCE_DIR}/src/*.cpp")
file(COPY ${CMAKE_SOURCE_DIR}/share/input.txt 
     DESTINATION ${CMAKE_BINARY_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/share/FFCM-2.yaml 
     DESTINATION ${CMAKE_BINARY_DIR})

if (APPLE)
    find_program(BREW_EXECUTABLE brew)

    if (BREW_EXECUTABLE)
        execute_process(
            COMMAND ${BREW_EXECUTABLE} --prefix
            OUTPUT_VARIABLE HOMEBREW_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        
        if (HOMEBREW_PREFIX)
            set(HOMEBREW_INCLUDE_DIR "${HOMEBREW_PREFIX}/include")
            set(HOMEBREW_LIB_DIR "${HOMEBREW_PREFIX}/lib")
            message(STATUS "Homebrew include directory: ${HOMEBREW_INCLUDE_DIR}")
            message(STATUS "Homebrew lib directory: ${HOMEBREW_LIB_DIR}")
        else()
            message(FATAL_ERROR "Homebrew prefix not found")
        endif()
    else()
        message(FATAL_ERROR "Homebrew not found")
    endif()
endif()

# set(CMAKE_BUILD_RPATH "${PROJECT_SOURCE_DIR}/lib/cantera/lib")

add_executable(cfb ${SRCS})

if (APPLE)
target_include_directories(
    cfb 
    PUBLIC ${PROJECT_SOURCE_DIR}/include
    PUBLIC ${PROJECT_SOURCE_DIR}/lib/cantera/include
    PUBLIC ${HOMEBREW_INCLUDE_DIR}
)

target_link_directories(
    cfb 
    PUBLIC ${PROJECT_SOURCE_DIR}/lib/cantera/lib
    PUBLIC ${HOMEBREW_LIB_DIR}
)


target_link_libraries(
    cfb 
    PUBLIC cantera_shared 
    PUBLIC "-framework Accelerate"
)
elseif(WIN32)
file(COPY ${CMAKE_SOURCE_DIR}/lib/cantera/bin/cantera_shared.dll 
     DESTINATION ${CMAKE_BINARY_DIR})

target_include_directories(
    cfb 
    PUBLIC ${PROJECT_SOURCE_DIR}/include
    PUBLIC ${PROJECT_SOURCE_DIR}/lib/cantera/include
)

target_link_directories(
    cfb 
    PUBLIC ${PROJECT_SOURCE_DIR}/lib/cantera/lib
    PUBLIC ${PROJECT_SOURCE_DIR}/lib/cantera/bin
)


target_link_libraries(
    cfb 
    PUBLIC cantera_shared 
)
    
endif()

