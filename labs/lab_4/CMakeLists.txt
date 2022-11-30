project(lab4)

add_executable(${PROJECT_NAME} "main.cpp")

set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 20)

target_link_libraries(${PROJECT_NAME} murlib sfml-system sfml-window sfml-graphics ImGui-SFML::ImGui-SFML implot)