# set(SOURCE
#   ${SOURCE}
#   PARENT_SCOPE
# )

# set(HEADERS
#    ${HEADERS}
#    ${CMAKE_CURRENT_SOURCE_DIR}/RawField.hpp
#    PARENT_SCOPE
# )

add_library(field
RawField.hpp
RawField.cpp
)

target_include_directories(field PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)