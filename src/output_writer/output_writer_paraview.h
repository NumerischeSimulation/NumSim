#pragma once

#include "output_writer/output_writer.h"
#include "discretization/1_discretization.h"

#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>

#include <memory>

/** Write *.vti files that can be viewed with ParaView.
 *  The mesh that can be visualized in ParaView corresponds to the mesh of the computational domain.
 *  All values are given for the nodes of the mesh, i.e., the corners of each cell.
 *  This means, values will be interpolated because the values are stored at positions given by the staggered grid.
 */
class OutputWriterParaview : 
  public OutputWriter
{
public:
  //! constructor
  OutputWriterParaview(std::shared_ptr<Discretization> discretization, std::shared_ptr<Partitioning> partitioning);

  //! write current velocities to file, filename is output_<count>.vti
  void writeFile(double currentTime);

private:

  vtkSmartPointer<vtkXMLImageDataWriter> vtkWriter_;   //< vtk writer to write ImageData
};
