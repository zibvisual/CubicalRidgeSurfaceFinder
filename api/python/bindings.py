import ctypes

# Load the shared library (C API)
lib = ctypes.CDLL("/srv/public/ploba/cubicalRidgeSurfaceFinder/build/librsfapi.so")

# Structure for points
class VecFloat(ctypes.Structure):
    _fields_ = [
        ("x", ctypes.c_float),
        ("y", ctypes.c_float),
        ("z", ctypes.c_float)
    ]

# Structure for triangles
class Triangle(ctypes.Structure):
    _fields_ = [
        ("v0", ctypes.c_size_t),
        ("v1", ctypes.c_size_t),
        ("v2", ctypes.c_size_t),
    ]

# Structure for surfaces
class Surface(ctypes.Structure):
    _fields_ = [
        ("points", ctypes.POINTER(VecFloat)), # Pointer to array of points
        ("num_points", ctypes.c_size_t),
        ("triangles", ctypes.POINTER(Triangle)), # Pointer to array of triangles
        ("num_triangles", ctypes.c_size_t),
    ]


# Incomplete class for `CubicalRidgeSurfaceFinder`
class CRSF(ctypes.Structure):
    pass

# Pointer to a CRSF object
CRSF_p = ctypes.POINTER(CRSF)

# Set the argument types and result types for all the API functions
lib.free_surface.argtypes = [ctypes.POINTER(Surface)]

lib.CRSF_new.restype = CRSF_p 

lib.CRSF_delete.argtypes = [CRSF_p]

lib.CRSF_setImage.argtypes = [
    CRSF_p,                  
    ctypes.POINTER(ctypes.c_float),   
    ctypes.c_int, ctypes.c_int, ctypes.c_int 
]

lib.CRSF_setThreshold.argtypes = [
    CRSF_p,
    ctypes.c_float,
    ctypes.c_float
]

lib.CRSF_getMinThreshold.argtypes = [CRSF_p]
lib.CRSF_getMinThreshold.restype = ctypes.c_float

lib.CRSF_getMaxThreshold.argtypes = [CRSF_p]
lib.CRSF_getMaxThreshold.restype = ctypes.c_float

lib.CRSF_addSeed.argtypes = [CRSF_p, ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_float]

lib.CRSF_addAutomaticSeeds.argtypes = [CRSF_p, ctypes.c_float, ctypes.c_float]
lib.CRSF_addAutomaticSeeds.restype = ctypes.c_int

lib.CRSF_clearSeedPoints.argtypes = [CRSF_p]

lib.CRSF_recalculate.argtypes = [CRSF_p]

lib.CRSF_compute_finalize.argtypes = [CRSF_p]
lib.CRSF_compute_finalize.restype = ctypes.POINTER(Surface)