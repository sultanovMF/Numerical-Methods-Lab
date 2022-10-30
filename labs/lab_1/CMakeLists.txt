project(lab1)

add_executable(${PROJECT_NAME} "main.cpp")

set(CMAKE_CXX_STANDARD 20)

target_link_libraries(${PROJECT_NAME} murlib sfml-system sfml-window sfml-graphics ImGui-SFML::ImGui-SFML implot)