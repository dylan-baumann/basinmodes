import ctypes
import os

# Load the installed Gmsh shared library (provided by libgmsh-dev)
libgmsh = ctypes.CDLL("libgmsh.so")

# Define return types and argument types for the C API functions
libgmsh.gmshInitialize.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_char_p), ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_int)]
libgmsh.gmshInitialize.restype = None

libgmsh.gmshOpen.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_int)]
libgmsh.gmshOpen.restype = None

libgmsh.gmshOptionSetNumber.argtypes = [ctypes.c_char_p, ctypes.c_double, ctypes.POINTER(ctypes.c_int)]
libgmsh.gmshOptionSetNumber.restype = None

libgmsh.gmshModelMeshGenerate.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_int)]
libgmsh.gmshModelMeshGenerate.restype = None

libgmsh.gmshWrite.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_int)]
libgmsh.gmshWrite.restype = None

libgmsh.gmshFinalize.argtypes = [ctypes.POINTER(ctypes.c_int)]
libgmsh.gmshFinalize.restype = None


def generate_msh_from_geo(input_filename: str) -> str:
    """
    Loads a .geo file via libgmsh ctypes binding and exports 
    a 2D triangular mesh in Gmsh 2.2 format.
    """
    if not input_filename.endswith(".geo"):
        raise ValueError("Input file must have a .geo extension.")

    # Derives 'foo.msh' from 'foo.geo'
    output_filename = os.path.splitext(input_filename)[0] + ".msh"
    
    ierr = ctypes.c_int(0)

    # 1. Initialize Gmsh API
    libgmsh.gmshInitialize(0, None, 1, 0, ctypes.byref(ierr))

    try:
        # Suppress INFO (5) and WARNING (2) messages; show ERRORS (1) only
        libgmsh.gmshOptionSetNumber(b"General.Verbosity", 1.0, ctypes.byref(ierr))

        # 2. Open input .geo script
        libgmsh.gmshOpen(input_filename.encode("utf-8"), ctypes.byref(ierr))

        # 3. Set algorithm to Frontal-Delaunay (2D triangles)
        libgmsh.gmshOptionSetNumber(b"Mesh.Algorithm", 6.0, ctypes.byref(ierr))

        # 4. Set format to Gmsh v2.2
        libgmsh.gmshOptionSetNumber(b"Mesh.MshFileVersion", 2.2, ctypes.byref(ierr))

        # 5. Generate 2D Mesh
        libgmsh.gmshModelMeshGenerate(2, ctypes.byref(ierr))

        # 6. Write to .msh file
        libgmsh.gmshWrite(output_filename.encode("utf-8"), ctypes.byref(ierr))

    finally:
        # 7. Ensure finalizing API call runs even if errors occur
        libgmsh.gmshFinalize(ctypes.byref(ierr))

    return output_filename