project(lab2)

add_executable(${PROJECT_NAME} "main.cpp")

set_property(TARGET lab2 PROPERTY CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -Wall -O2")

target_link_libraries(${PROJECT_NAME} murlib sfml-system sfml-window sfml-graphics ImGui-SFML::ImGui-SFML implot effolkronium_random)
