#include "../hddApi.hpp"



void HddApi::ParseJsonFile(std::string path2file)
{
  m_data.PathToSolverOptionFile(path2file);
}

void HddApi::SymbolicAssembling(int vals[], int size)
{
  m_data.SymbolicAssembling(vals, size);
}

void HddApi::SetDirichletDOFs(int vals[], int size)
{
  m_data.SetDirichletDOFs(vals, size);
}

int HddApi::FinalizeSymbolicAssembling()
{
  int neqSub = m_data.FinalizeSymbolicAssembling();
  return neqSub;
}


void HddApi::NumericAssembling(int glbIds[],
        double valK[], double valF[], int nDOF)
{
  m_data.NumericAssembling(glbIds, valK, valF, nDOF);
}

void HddApi::FinalizeNumericAssembling()
{
  m_data.FinalizeNumericAssembling();
}


void HddApi::Solve(double solution[], int nrows, int ncols)
{
//  Eigen::Map<Eigen::MatrixXd> solutionXd(solution, nrows, ncols);
//  m_data.Solve(solutionXd);
    m_data.Solve(solution, nrows, ncols);
}

void HddApi::Finalize()
{
  m_data.Finalize();
}
