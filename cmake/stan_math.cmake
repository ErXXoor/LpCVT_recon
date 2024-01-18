include(FetchContent)
FetchContent_Declare(
        stanmath
        GIT_REPOSITORY https://github.com/stan-dev/math.git
        GIT_TAG v4.7.0
)
FetchContent_MakeAvailable(stanmath)

find_package(Boost 1.81.0 REQUIRED)
find_package(TBB 2021 REQUIRED)
find_package(Eigen3 3.4.0 REQUIRED)
option(BUILD_TESTING OFF)
option(EIGEN_BUILD_PKGCONFIG OFF)

add_library(stan-math INTERFACE)
target_compile_definitions(stan-math INTERFACE NO_FPRINTF_OUTPUT BOOST_DISABLE_ASSERTS TBB_INTERFACE_NEW _REENTRANT STAN_MATH_REV_CORE_INIT_CHAINABLESTACK_HPP STAN_THREADS)
target_link_libraries(stan-math INTERFACE Eigen3::Eigen TBB::tbb -lstdc++ -lm)
target_include_directories(stan-math INTERFACE ${stanmath_SOURCE_DIR} ${Boost_INCLUDE_DIRS})
add_library(stan::math ALIAS stan-math)