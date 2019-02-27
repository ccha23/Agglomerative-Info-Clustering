cmake_minimum_required(VERSION 2.8)

project( AIC )
add_compile_options(-std=c++11)

include_directories(.)
find_package( OpenCV REQUIRED )
add_executable( AIC cart_demo.cpp )
target_link_libraries( AIC ${OpenCV_LIBS} )
