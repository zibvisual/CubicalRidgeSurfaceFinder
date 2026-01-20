import ctypes
import numpy as np

lib = ctypes.CDLL("/srv/public/ploba/cubicalRidgeSurfaceFinder/build/librsfapi.so")

class VecFloat(ctypes.Structure):
    _fields_ = [
        ("x", ctypes.c_float),
        ("y", ctypes.c_float),
        ("z", ctypes.c_float),
    ]

class Triangle(ctypes.Structure):
    _fields_ = [
        ("v0", ctypes.c_size_t),
        ("v1", ctypes.c_size_t),
        ("v2", ctypes.c_size_t),
        ("patch", ctypes.c_uint64),
    ]

class Surface(ctypes.Structure):
    _fields_ = [
        ("points", ctypes.POINTER(VecFloat)),
        ("num_points", ctypes.c_size_t),
        ("points_capacity", ctypes.c_size_t),
        ("triangles", ctypes.POINTER(Triangle)),
        ("num_triangles", ctypes.c_size_t),
        ("triangles_capacity", ctypes.c_size_t),
    ]

# the output of the function is set as a void pointer
lib.CRSF_new.restype = ctypes.c_void_p 

# set the argument data types
lib.CRSF_setThreshold.argtypes = [
    ctypes.c_void_p,
    ctypes.c_float,
    ctypes.c_float
]

# set the inputs and output of setImage
lib.CRSF_setImage.argtypes = [
    ctypes.c_void_p,                  
    ctypes.POINTER(ctypes.c_float),   
    ctypes.c_int, ctypes.c_int, ctypes.c_int 
]
lib.CRSF_setImage.restype = None


lib.CRSF_getMinThreshold.argtypes = [ctypes.c_void_p]
lib.CRSF_getMinThreshold.restype = ctypes.c_float

lib.CRSF_delete.argtypes = [ctypes.c_void_p]

# addSeed
lib.CRSF_addSeed.argtypes = [ctypes.c_void_p, ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_float]
lib.CRSF_addSeed.restype = None

# addAutomaticSeed
lib.CRSF_addAutomaticSeeds.argtypes = [ctypes.c_void_p, ctypes.c_float, ctypes.c_float]
lib.CRSF_addAutomaticSeeds.restype = ctypes.c_int

# free surface memory
lib.free_surface.argtypes = [ctypes.POINTER(Surface)]
lib.free_surface.restype = None

# compute
lib.CRSF_compute.argtypes = [ctypes.c_void_p]
lib.CRSF_compute.restype = None

# recalculate
lib.CRSF_recalculate.argtypes = [ctypes.c_void_p]
lib.CRSF_recalculate.restype = None

# getPatchedSurface
lib.CRSF_getPatchedSurface.argtypes = [ctypes.c_void_p] 
lib.CRSF_getPatchedSurface.restype = ctypes.POINTER(Surface)

# finalize
lib.CRSF_finalize.argtypes = [ctypes.c_void_p]
lib.CRSF_finalize.restype = ctypes.POINTER(Surface)