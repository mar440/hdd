#pragma once

#include <vtkUnstructuredGrid.h>
#include <string>
#include <vector>
#include <map>

#include "types.hpp"

class vtkUnstructuredGrid;

class Mesh
{

  public:
    Mesh();
    ~Mesh();
    void print();
    int generateMesh(int ne, int ns, int nl);
    vtkUnstructuredGrid* squareMesh(double,double,double=0, double=0);
    void writeMesh(vtkUnstructuredGrid *vtkVolumeMesh,
        std::string filename, bool asciiOrBinaryVtu);

    std::vector<int>& getDirDOFs(){return m_DirichletDofs;}
    void addSolution(Eigen::VectorXd& solution);
    void addSolution(vtkUnstructuredGrid* ug, Eigen::VectorXd& solution);
    std::vector<int> getNeighboursRanks();


    void extractSubdomainMesh();
    int createMappingVectors();
    vtkUnstructuredGrid* getGlobalMesh(){return m_mesh;}
    vtkUnstructuredGrid* getSubdomainMesh(){return m_subdomainMesh;}

  private:
      vtkUnstructuredGrid* m_mesh;
      vtkUnstructuredGrid* m_subdomainMesh;
      int m_nex, m_ney;
      int m_nsx, m_nsy;
      int m_nlx, m_nly;
      int m_rank;
      std::vector<int> m_DirichletDofs;
      std::vector<int> m_l2g;
      std::map<int,int>m_g2l;

};
