project (murlib)

add_library(${PROJECT_NAME} "src/murlib.cpp")
add_library(murlib::interpolation ALIAS ${PROJECT_NAME})

target_include_directories( ${PROJECT_NAME}
    PUBLIC ${PROJECT_SOURCE_DIR}/include
)