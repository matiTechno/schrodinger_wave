add_executable(brick_crusher main.cpp game_scene.cpp)
target_link_libraries(brick_crusher mits -lfftw3)
file(COPY data DESTINATION .)
