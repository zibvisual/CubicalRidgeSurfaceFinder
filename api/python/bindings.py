import ctypes
import numpy as np

lib = ctypes.CDLL("/srv/public/ploba/cubicalRidgeSurfaceFinder/build/librsfapi.so")

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
    ctypes.c_void_p,                  # puntero al objeto C++
    ctypes.POINTER(ctypes.c_float),   # puntero a los datos de la imagen
    ctypes.c_int, ctypes.c_int, ctypes.c_int  # width, height, depth
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

# Test the code
rsf = lib.CRSF_new()
lib.CRSF_setThreshold(rsf, 0.2, 0.8)

#print(f"min threshold: {lib.CRSF_getMinThreshold(rsf):.2f}" )

# Create an 3D image
depth, height, width = 5, 4, 3
img = np.random.rand(depth, height, width).astype(np.float32)
img = np.ascontiguousarray(img) # Make sure that it is contiguous in memory

# Obtain the pointer
img_ptr = img.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

# Call the function
lib.CRSF_setImage(rsf, img_ptr, width, height, depth)
#print('Function called without errors')

# add a seed point
x, y, z = 1.5, 2.5, 3.5
distance = 2.0
lib.CRSF_addSeed(rsf, x, y, z, distance)
#print(f"Seed added at ({x}, {y}, {z}) with distance {distance}.")

#add automatic seed
index = lib.CRSF_addAutomaticSeeds(rsf, 0.5, 2)
#if index >= 0:
    #print(f"Automatic seed added at index {index}")
#else:
    #print("No automatic seed added")

# clear seed points
lib.CRSF_clearSeedPoints(rsf)

# compute and recalculate
lib.CRSF_compute(rsf)
lib.CRSF_recalculate(rsf)

# Clean
lib.CRSF_delete(rsf)