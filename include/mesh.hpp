#pragma once

#include <vtkUnstructuredGrid.h>
#include <string>
#include <vector>

class vtkUnstructuredGrid;

class Mesh
{

  public:
    Mesh(){}
    ~Mesh();
    void print();
    int generateMesh(int ne, int ns, int nl);
    vtkUnstructuredGrid* squareMesh(double,double,int,int);
    void writeMesh(vtkUnstructuredGrid *vtkVolumeMesh,
        std::string filename, bool asciiOrBinaryVtu);
    vtkUnstructuredGrid* m_mesh;

    std::vector<int>& getDirDOFs(){return m_DirichletDofs;}

  private:
      int m_ne, m_ns, m_nl;
      int m_rank;
      std::vector<int> m_DirichletDofs;

};
