#include <iostream>
#include <cstring>
#include <fstream>
#include <cstdlib>
#include <vector>


int main(int argc, char *argv[])
{


  if (argc != 3)
  {
    std::cout << "  \033[1;33mUse as " << argv[0]
              << " inputfile (plot 3D format ) outputfile (gmsh format) \033[1;m\n";
    exit(2);
  }

  //=============================================================
  //                          Variables
  //=============================================================

  std::string infilename = std::string(argv[1]);
  std::string outfilename = std::string(argv[2]);

  std::ifstream infile;
  std::ofstream outfile;

  std::vector<double> x_coord(0), y_coord(0);
  // Number of layers
  unsigned nl;
  // Number of points in the x-direction
  unsigned idim;
  // Number of points in the y-direction
  unsigned jdim;

  //=============================================================
  //                          File reading
  //=============================================================

  infile.open(infilename.c_str());

  // Read the number of layers in the third (z-) dimension
  infile >> nl;

  if (nl != 1)
  {
    std::cerr << "Number of layers has to be 1 (2D mesh), but here nl = " << nl << ". Aborting."
              << std::endl;
    exit(3);
  }

  infile >> idim;
  infile >> jdim;

  // Total number of nodes in the mesh
  const unsigned npoints = idim * jdim;

  std::cout << "Mesh size = [" << idim << " x " << jdim << " x " << nl << "]" << std::endl;

  x_coord.resize(npoints);
  y_coord.resize(npoints);

  for (unsigned i = 0; i < npoints; ++i)
  {
    infile >> x_coord[i];
  }

  for (unsigned j = 0; j < npoints; ++j)
  {
    infile >> y_coord[j];
  }

  infile.close();

  outfile.open(outfilename.c_str());

  outfile << "$MeshFormat" << std::endl;
  outfile << "2.2 0 8" << std::endl;
  outfile << "$EndMeshFormat" << std::endl;
  outfile << "$PhysicalNames" << std::endl;
  outfile << "6" << std::endl;
  outfile << "1 1 \"inlet\"" << std::endl;
  outfile << "1 2 \"bottom\"" << std::endl;
  outfile << "1 3 \"wall\"" << std::endl;
  outfile << "1 4 \"outlet\"" << std::endl;
  outfile << "1 5 \"farfield\"" << std::endl;
  outfile << "2 6 \"fluid\"" << std::endl;
  outfile << "$EndPhysicalNames" << std::endl;
  outfile << "$Nodes" << std::endl;
  outfile << npoints << std::endl;

  outfile.precision(14);
  outfile.setf(std::ios::fixed);

  for (unsigned i = 0; i < npoints; ++i)
  {
    outfile << i + 1 << " " << x_coord[i] << " " << y_coord[i] << " 0.0" << std::endl;
  }

  outfile << "$EndNodes" << std::endl;
  outfile << "$Elements" << std::endl;

  // Total number of quadrilateral elements
  const unsigned nelem = (idim - 1) * (jdim - 1);

  // Total number of boundary segments
  const unsigned nbdry_lines = 2 * ((idim - 1) + (jdim - 1));

  outfile << nelem + nbdry_lines << std::endl;

  unsigned ielem = 1;

  // Write the connectivity for the inlet boundary
  for (int j = (jdim - 2); j >= 0; --j)
  {
    const unsigned v0 = (j + 1) * idim;
    const unsigned v1 = j * idim;
    outfile << ielem << " 1 2 1 1 " << v0 + 1 << " " << v1 + 1 << std::endl;
    ielem++;
  }

  // The input file contains first all x-coordinates, then all y-coordinates
  // going from left to right in rows (top to bottom).
  // This means that the first idim points are points on the bottom if the
  // domain and some of them lie before the plate ( x <= 0 ) and some of them
  // lie on the plate ( x>= 0)

  // Find the index of the first point on the plate (x = 0)

  unsigned first_plate_pt = 0;

  while ((first_plate_pt < idim) && (x_coord[first_plate_pt] < -1.e-14))
  {
    first_plate_pt++;
  }

  // Write the connectivity for the bottom boundary
  for (unsigned i = 0; i < first_plate_pt; ++i)
  {
    const unsigned v0 = i;
    const unsigned v1 = i + 1;
    outfile << ielem << " 1 2 2 2 " << v0 + 1 << " " << v1 + 1 << std::endl;
    ielem++;
  }

  // Write the connectivity for the wall (plate) boundary
  for (unsigned i = first_plate_pt; i < (idim - 1); ++i)
  {
    const unsigned v0 = i;
    const unsigned v1 = i + 1;
    outfile << ielem << " 1 2 3 3 " << v0 + 1 << " " << v1 + 1 << std::endl;
    ielem++;
  }

  // Write the connectivity for the outlet boundary
  for (unsigned j = 0; j < (jdim - 1); ++j)
  {
    const unsigned v0 = j * idim + (idim - 1);
    const unsigned v1 = (j + 1) * idim + (idim - 1);
    ;
    outfile << ielem << " 1 2 4 4 " << v0 + 1 << " " << v1 + 1 << std::endl;
    ielem++;
  }

  // Write the connectivity for the farfield boundary
  for (int i = (idim - 1); i > 0; --i)
  {
    const unsigned v0 = (jdim - 1) * idim + i;
    const unsigned v1 = (jdim - 1) * idim + i - 1;
    outfile << ielem << " 1 2 5 5 " << v0 + 1 << " " << v1 + 1 << std::endl;
    ielem++;
  }

  for (unsigned j = 0; j < (jdim - 1); ++j)
  {
    for (unsigned i = 0; i < (idim - 1); ++i)
    {
      // Four indexes of each quad: the southwest, southeast, northeast and northwest corner
      const unsigned sw = (j * idim) + i;
      const unsigned se = sw + 1;
      const unsigned nw = ((j + 1) * idim) + i;
      const unsigned ne = nw + 1;
      // When outputting the vertex indices, add 1 because gmsh numbers all entities starting from 1
      // and not 0 like C/C++
      outfile << ielem << " 3 2 6 6 " << sw + 1 << " " << se + 1 << " " << ne + 1 << " " << nw + 1
              << std::endl;
      ielem++;
    }
  }

  outfile << "$EndElements" << std::endl;

  outfile.close();

  return 0;
}
