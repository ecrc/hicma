# Dataset

## Mesh Deformation Application

Dataset is available in KAUST repository: https://repository.kaust.edu.sa/handle/10754/664938. 
Add mesh file name to `--mesh_file=<mesh file name>` parameter.

## Acoustic Scattering Application

Dataset is available in KAUST repository: https://repository.kaust.edu.sa/handle/10754/664400.
Add mesh file name to `--mesh_file=<mesh file name>` parameter and the file containing interpolation points information to `--nipp=<number of interpolation points>` parameter.


For more information on the dataset please refer to the readme files in the data repositories.

## Testing Mesh Deformation and Acoustic Scattering Applications

This [cmake file](timing/CMakeLists.txt) contains sample commands to run the mesh deformation and acoustic scattering applications.
