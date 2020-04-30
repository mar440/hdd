#pragma once

#include <vtkUnstructuredGrid.h>
#include <string>
#include <vector>
#include <map>
#include <boost/property_tree/ptree.hpp>

#include "types.hpp"

class vtkUnstructuredGrid;
class Data;

class Mesh
{
  public:
    Mesh(){}
    ~Mesh();
    void print();
    int GenerateMesh(int rank, boost::property_tree::ptree);
    vtkUnstructuredGrid* squareMesh(int);
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
    void SaveDecomposedMesh();

  private:
      vtkUnstructuredGrid* m_mesh;
      vtkUnstructuredGrid* m_subdomainMesh;
      double m_Lx2D, m_Ly2D;
      double m_x0, m_y0;
      int m_nex, m_ney;
      int m_nsx, m_nsy;
      int m_nsxOneSub, m_nsyOneSub;
      int m_nlx, m_nly;
      int m_rank;
      std::vector<int> m_DirichletDofs;
      std::vector<int> m_l2g;
      std::map<int,int>m_g2l;
      void ddm_metis();

};
