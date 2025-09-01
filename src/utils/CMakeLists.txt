# set(SOURCE
#   ${SOURCE}
#   ${CMAKE_CURRENT_SOURCE_DIR}/Vec.cpp
#   PARENT_SCOPE
# )

# set(HEADERS
#    ${HEADERS}
#    ${CMAKE_CURRENT_SOURCE_DIR}/BBox.hpp
#    ${CMAKE_CURRENT_SOURCE_DIR}/Dims.hpp
#    ${CMAKE_CURRENT_SOURCE_DIR}/Vec.hpp
#    ${CMAKE_CURRENT_SOURCE_DIR}/ProgressbarReportHelper.hpp
#    PARENT_SCOPE
# )

# add_library(utils
# ${SOURCE}
# )

add_library(utils
BBox.hpp
Dims.hpp
Direction.hpp
Vec.hpp
Vec.cpp
ProgressbarReportHelper.hpp
)

target_include_directories(utils PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)