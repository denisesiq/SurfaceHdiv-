# Project to implement Hdiv++ on surface 
file(GLOB header *.h)
file(GLOB source *.cpp)

set(SurfaceHdiv++SRC ${header} ${source})

add_executable(Surface ${SurfaceHdiv++SRC})

target_link_libraries(Surface pz)

configure_file(BasicMesh.msh BasicMesh.msh COPYONLY)
