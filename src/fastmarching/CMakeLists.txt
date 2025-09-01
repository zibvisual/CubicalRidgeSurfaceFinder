# set(SOURCE
#   ${SOURCE}
#   ${CMAKE_CURRENT_SOURCE_DIR}/file1.cpp
#   ${CMAKE_CURRENT_SOURCE_DIR}/file2.cpp
#   PARENT_SCOPE
# )

set(HEADERS
   ${HEADERS}
   ${CMAKE_CURRENT_SOURCE_DIR}/FastMarching.hpp
   ${CMAKE_CURRENT_SOURCE_DIR}/FastMarchingAggregateObserver.hpp
   ${CMAKE_CURRENT_SOURCE_DIR}/FastMarchingFlowObserver.hpp
   ${CMAKE_CURRENT_SOURCE_DIR}/FastMarchingGeodesicDistanceObserver.hpp
   ${CMAKE_CURRENT_SOURCE_DIR}/FastMarchingObserver.hpp
   ${CMAKE_CURRENT_SOURCE_DIR}/FastMarchingUtilities.hpp
   PARENT_SCOPE
)