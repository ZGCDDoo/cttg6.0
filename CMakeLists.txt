#Inspired by https://github.com/dmonopoly/gtest-cmake-example
#           https://github.com/misaka-10032/boost-mpi-examples

cmake_minimum_required(VERSION 2.8.2)

option(TEST "Build all tests." ON) # Makes boolean 'test' available.
option(HOME "home build" ON) #builds for linux, default
option(MAC "mac build" OFF)
option(GRAHAM "graham build" OFF)
option(MP2 "mp2 build" OFF)
option(MPI_BUILD "build with mpi" OFF)

if (${MP2})
    SET(CMAKE_CXX_COMPILER g++)
elseif (${MPI_BUILD})
    SET(CMAKE_CXX_COMPILER mpic++)
else ()
    SET(CMAKE_CXX_COMPILER clang++)
endif()


set(PROJECT_NAME CT-INT)
project("${PROJECT_NAME}")

set(CMAKE_CXX_FLAGS " -Wall -Wextra -std=c++14 -pipe -Ofast -march=native -fomit-frame-pointer") 
if(${MPI_BUILD})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVEMPI")
endif()

set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/src/Includes)

find_package(LAPACK REQUIRED)

#-------------HOME--------------------------------------------------
if (${HOME})
    find_package(Boost REQUIRED mpi serialization filesystem system)
    set(LIBRARIES_EXEC ${LAPACK_LIBRARIES} ${Boost_LIBRARIES} armadillo)
endif ()
#-------------HOME--------------------------------------------------

if(${MAC})
    find_package(Boost REQUIRED mpi serialization filesystem system)
    set(LIBRARIES_EXEC lapack blas ${Boost_LIBRARIES} armadillo)
endif()
 


#-------------GRAHAM--------------------------------------------------
if (${GRAHAM})
    set(TEST OFF)
    set(CMAKE_EXE_LINKER_FLAGS "-L${MKLROOT}/lib/intel64")
    set(LIBRARIES_EXEC armadillo mkl_gf_lp64 mkl_sequential mkl_core boost_filesystem boost_mpi)
endif ()
#-------------GRAHAM--------------------------------------------------



#-------------MP2--------------------------------------------------
if (${MP2})
    set(TEST OFF)
    SET(CMAKE_EXE_LINKER_FLAGS "-L${MKLROOT}/lib/intel64 -L/opt/boost64/1.65.1/lib -L/opt/armadillo/8.300.0/usr/lib64")
    include_directories(/opt/armadillo/8.300.0/usr/include /opt/boost64/1.65.1/include/boost/  ${MKLROOT}/include)
    set(LIBRARIES_EXEC armadillo mkl_gf_lp64 mkl_sequential mkl_core boost_filesystem boost_mpi)
endif ()
#-------------MP2--------------------------------------------------


################################

# Executables

################################
set(EXECECUTABLES
    dmft
    dmft_aux_sub
    dmft_aux
    dmft_sub
    cdmft_triangle2x2
    cdmft_triangle2x2_aux
    cdmft_triangle2x2_sub
    cdmft_triangle2x2_aux_sub
    cdmft_square2x2
    cdmft_square2x2_aux
    cdmft_square2x2_sub
    cdmft_square2x2_aux_sub
    cdmft_square4x4
    cdmft_square4x4_aux
    cdmft_square4x4_aux_sub
    cdmft_square4x4_sub
)

foreach(executable ${EXECECUTABLES})
    list(APPEND ${executable}_src ${PROJECT_SOURCE_DIR}/src/${executable}.cpp)
    add_executable(${executable}
                ${${executable}_src}
                  )
endforeach()

foreach(executable ${EXECECUTABLES})
    target_link_libraries(${executable}  ${LIBRARIES_EXEC})
    install(TARGETS ${executable} DESTINATION $ENV{HOME}/bin)
endforeach()


#--------------------end Executables ----------------------------------------------------------


################################

# Testing

################################

if (${TEST})
    enable_testing()
    
    # Download and unpack googletest at configure time
    configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
                    RESULT_VARIABLE result
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
    if(result)
        message(FATAL_ERROR "CMake step for googletest failed: ${result}")
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
                    RESULT_VARIABLE result
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
    if(result)
        message(FATAL_ERROR "Build step for googletest failed: ${result}")
    endif()

    # Prevent overriding the parent project's compiler/linker
    # settings on Windows
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    # Add googletest directly to our build. This defines
    # the gtest and gtest_main targets.
    add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
                    ${CMAKE_BINARY_DIR}/googletest-build
                    EXCLUDE_FROM_ALL)

    # The gtest/gtest_main targets carry header search path
    # dependencies automatically when using CMake 2.8.11 or
    # later. Otherwise we have to add them here ourselves.
    if (CMAKE_VERSION VERSION_LESS 2.8.11)
        include_directories("${gtest_SOURCE_DIR}/include")
    endif()

    set(TESTS
    FourierTests
    FillingAndDoccTests
    GreenBinningTests
    GreenMatTests
    GreenTauTests
    H0TriangleTests
    IntegratorTests
    ModelTriangle2x2Tests
    UtilitiesTests
    SelfConsistencyTests
    MarkovChainTests
    MarkovChainSquare2x2Tests
    MarkovChainAuxTests
    MarkovChainSubMatrixSquare2x2Tests
    MarkovChainAuxSubMatrixTests
    MarkovChainSubMatrixTests
    MatrixTests
    ObservablesTests
    )

    set(TESTS_DIR ${PROJECT_SOURCE_DIR}/test)

    foreach(test ${TESTS})
        list(APPEND ${test}_src ${TESTS_DIR}/${test}.cpp)
        add_executable(${test}
                    ${${test}_src}
                    )
    endforeach()

    set(LIBRARIES_TEST ${LIBRARIES_EXEC} gtest gtest_main)
    foreach(test ${TESTS})
        target_link_libraries(${test} ${LIBRARIES_TEST})
    endforeach()

    foreach(test ${TESTS})
        add_test(NAME ${test} COMMAND ./${test})
    endforeach()

    #add_test(NAME ObservablesDMFTTests COMMAND "python3" ../test/obs_dmft_tests.py)
    #add_test(NAME ObservablesCDMFTTests COMMAND "python3" ../test/obs_cdmft_tests.py)
    #add_test(NAME ObservablesCDMFTBinningTests COMMAND "python3" ../test/obs_cdmft_binning_tests.py)
endif()


#--------------------end Testing ----------------------------------------------------------

