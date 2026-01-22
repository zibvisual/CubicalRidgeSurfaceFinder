import ctypes

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

class CubicalRidgeSurfaceFinderAPI:
    """
    This class provides an interface to the C++ API for Ridge Surface Finder.
    """
    def __init__(self, lib_path):
        """
        Initializes the API by loading the specified shared library.

        :param lib_path: Path to the compiled C++ library.
        :raises RuntimeError: If the C++ library cannot be loaded.
        """
        try:
            self.lib = ctypes.CDLL(lib_path)
            self._setup_functions()
        except OSError as e:
            raise RuntimeError(f"Failed to load the C++ library: {e}")

    def _setup_functions(self):
        """
        Sets up argument types and return types for all the C++ functions
        exposed in the API to make them compatible with ctypes.
        """
        self.lib.free_surface.argtypes = [ctypes.POINTER(Surface)]

        self.lib.CRSF_new.restype = CRSF_p 

        self.lib.CRSF_delete.argtypes = [CRSF_p]

        self.lib.CRSF_setImage.argtypes = [
            CRSF_p,                  
            ctypes.POINTER(ctypes.c_float),   
            ctypes.c_int, ctypes.c_int, ctypes.c_int 
        ]

        self.lib.CRSF_setThreshold.argtypes = [
            CRSF_p,
            ctypes.c_float,
            ctypes.c_float
        ]

        self.lib.CRSF_getMinThreshold.argtypes = [CRSF_p]
        self.lib.CRSF_getMinThreshold.restype = ctypes.c_float

        self.lib.CRSF_getMaxThreshold.argtypes = [CRSF_p]
        self.lib.CRSF_getMaxThreshold.restype = ctypes.c_float

        self.lib.CRSF_addSeed.argtypes = [CRSF_p, ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_float]

        self.lib.CRSF_addAutomaticSeeds.argtypes = [CRSF_p, ctypes.c_float, ctypes.c_float]
        self.lib.CRSF_addAutomaticSeeds.restype = ctypes.c_int

        self.lib.CRSF_clearSeedPoints.argtypes = [CRSF_p]

        self.lib.CRSF_recalculate.argtypes = [CRSF_p]

        self.lib.CRSF_compute_finalize.argtypes = [CRSF_p]
        self.lib.CRSF_compute_finalize.restype = ctypes.POINTER(Surface)

    def free_surface(self, surface_p):
        """
        Frees the memory allocated for a `Surface_C` object.
        
        :param surface_p: Pointer to the `Surface_C` object to be freed.
        :raises ValueError: If the surface pointer is invalid or None.
        :raises RuntimeError: If the memory freeing fails.
        """
        if surface_p is None:
            raise ValueError("Surface pointer is None. Cannot free memory.")

        try:
            self.lib.free_surface(surface_p)
        except Exception as e:
            raise RuntimeError(f"Failed to free the surface memory: {e}")

    def create_instance(self):
        """
        Creates a new instance of `CubicalRidgeSurfaceFinder`.

        :return: Pointer to the newly created `CubicalRidgeSurfaceFinder` instance.
        :raises RuntimeError: If the instance creation fails.
        """
        try:
            return self.lib.CRSF_new()
        except Exception as e:
            raise RuntimeError(f"Error creating CRSF instance: {e}")

    def delete_instance(self, crsf_p):
        """
        Deletes an instance of `CubicalRidgeSurfaceFinder`.

        :param crsf_p: Pointer to the `CubicalRidgeSurfaceFinder` instance to delete.
        :raises RuntimeError: If instance deletion fails.
        """
        try:
            self.lib.CRSF_delete(crsf_p)
        except Exception as e:
            raise RuntimeError(f"Error deleting CRSF instance: {e}")

    def set_image(self, crsf_p, image, width, height, depth):
        """
        Sets the input image for the `CubicalRidgeSurfaceFinder`.

        :param crsf_p: Pointer to the `CubicalRidgeSurfaceFinder` instance.
        :param image: Pointer to the image data (1D array of floats).
        :param width: Width of the image in voxels.
        :param height: Height of the image in voxels.
        :param depth: Depth of the image in voxels.
        :raises ValueError: If setting the image fails.
        """
        try:
            self.lib.CRSF_setImage(crsf_p, image, width, height, depth)
        except Exception as e:
            raise ValueError(f"Error setting the image: {e}")

    def set_threshold(self, crsf_p, min_threshold, max_threshold):
        """
        Sets the minimum and maximum thresholds for the `CubicalRidgeSurfaceFinder`.

        :param crsf_p: Pointer to the `CubicalRidgeSurfaceFinder` instance.
        :param min_threshold: Minimum threshold value.
        :param max_threshold: Maximum threshold value.
        :raises ValueError: If setting the thresholds fails.
        """
        try:
            self.lib.CRSF_setThreshold(crsf_p, min_threshold, max_threshold)
        except Exception as e:
            raise ValueError(f"Error setting the thresholds: {e}")

    def get_min_threshold(self, crsf_p):
        """
        Gets the current minimum threshold for the `CubicalRidgeSurfaceFinder`.

        :param crsf_p: Pointer to the `CubicalRidgeSurfaceFinder` instance.
        :return: The minimum threshold value.
        :raises RuntimeError: If getting the minimum threshold fails.
        """
        try:
            return self.lib.CRSF_getMinThreshold(crsf_p)
        except Exception as e:
            raise RuntimeError(f"Error getting the minimum threshold: {e}")

    def get_max_threshold(self, crsf_p):
        """
        Gets the current maximum threshold for the `CubicalRidgeSurfaceFinder`.

        :param crsf_p: Pointer to the `CubicalRidgeSurfaceFinder` instance.
        :return: The maximum threshold value.
        :raises RuntimeError: If getting the maximum threshold fails.
        """
        try:
            return self.lib.CRSF_getMaxThreshold(crsf_p)
        except Exception as e:
            raise RuntimeError(f"Error getting the maximum threshold: {e}")

    def add_seed(self, crsf_p, x, y, z, distance):
        """
        Adds a seed point to the `CubicalRidgeSurfaceFinder`.

        :param crsf_p: Pointer to the `CubicalRidgeSurfaceFinder` instance.
        :param x: X-coordinate of the seed point in 3D space.
        :param y: Y-coordinate of the seed point in 3D space.
        :param z: Z-coordinate of the seed point in 3D space.
        :param distance: The maximum distance for Fast Marching propagation.
        :raises ValueError: If adding the seed fails.
        """
        try:
            self.lib.CRSF_addSeed(crsf_p, x, y, z, distance)
        except Exception as e:
            raise ValueError(f"Error adding seed point: {e}")

    def add_automatic_seeds(self, crsf_p, minDistanceRel, maxDistanceNewSeed):
        """
        Docstring for add_automatic_seeds
        
        :param crsf_p: Pointer to the `CubicalRidgeSurfaceFinder` instance.
        :param minDistanceRel: The minimum relative distance required to consider a candidate as a valid seed point.
        :param maxDistanceNewSeed: The maximum distance to propagate during the Fast Marching algorithm by the new seed.
        :raises ValueError: If no valid seed point is found.
        :raises RuntimeError: If adding the automatic seeds fails.
        """
        try:
            index = self.lib.CRSF_addAutomaticSeeds(crsf_p, minDistanceRel, maxDistanceNewSeed)
            if index == -1:
                raise ValueError("No valid seed point found.")
            return index
        except Exception as e:
            raise RuntimeError(f"Error adding automatic seed points: {e}")

    def clear_seed_points(self, crsf_p):
        """
        Clears all seed points from the `CubicalRidgeSurfaceFinder`.

        :param crsf_p: Pointer to the `CubicalRidgeSurfaceFinder` instance.
        :raises RuntimeError: If clearing the seed points fails.
        """
        try:
            self.lib.CRSF_clearSeedPoints(crsf_p)
        except Exception as e:
            raise RuntimeError(f"Error clearing seed points: {e}")
    def recalculate(self, crsf_p):
        """
        Recalculates the ridge surface patches based on the current seed points.

        :param crsf_p: Pointer to the `CubicalRidgeSurfaceFinder` instance.
        :raises RuntimeError: If the recalculation fails.
        """
        try:
            self.lib.CRSF_recalculate(crsf_p)
        except Exception as e:
            raise RuntimeError(f"Error recalculating ridge surface: {e}")

    def compute_finalize(self, crsf_p):
        """
        Computes the surface using Fast Marching algorithm and merges all the patches.
        
        :param crsf_p: Pointer to the `CubicalRidgeSurfaceFinder` instance.
        :raises RuntimeError: If the computation fails.
        """
        try:
            return self.lib.CRSF_compute_finalize(crsf_p)
        except Exception as e:
            raise RuntimeError(f"Error computing ridge surface: {e}")
