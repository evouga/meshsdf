Converts an .obj into a signed distance field.

First translated and scales a mesh into a cube of size [0,L]^3. Then computes the SDF on an m x m x m grid within that cube.

## Dependencies

Eigen (header-only): http://eigen.tuxfamily.org/index.php?title=Main_Page
Libigl (header-only): get it using

    git clone --recursive https://github.com/libigl/libigl.git

and set LIBIGL environment variable to point to the place where you clone libigl, so that CMake can find it.
