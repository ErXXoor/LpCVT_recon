include(FetchContent)
FetchContent_Declare(
        geogram
        GIT_REPOSITORY https://github.com/BrunoLevy/geogram.git
        GIT_TAG        v1.8.6
)
FetchContent_MakeAvailable(geogram)
target_compile_features(geogram PUBLIC cxx_std_11)
