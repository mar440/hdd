#include "../include/hddTime.hpp"
#include <iostream>



HddTime::HddTime(std::string name)
{
  m_open = true;
  m_name = name;
  std::cout << "Start  " << m_name << "...\n";
  m_start =  std::chrono::steady_clock::now();
  m_restart = m_start;
}

HddTime::~HddTime()
{
  CloseTime();
}

void HddTime::CloseTime()
{
  if (m_open!=0)
  {
    m_end = std::chrono::steady_clock::now();
    m_elapsedTime  = m_end - m_start;
    std::cout << "Finish " << m_name << ": ";
    std::cout << std::fixed << std::setprecision(2) <<
      m_elapsedTime.count() << " s\n";
    m_open = 0;
  }
}


void HddTime::Capture(std::string msg)
{
  if (m_open == 0)
  {
    m_start =  std::chrono::steady_clock::now();
    m_restart = m_start;
    m_open = 1;
  }
  else
  {

    // previous restart (or start)
    m_end = std::chrono::steady_clock::now();
    m_elapsedTime  = m_end - m_restart;

    m_restart = std::chrono::steady_clock::now(); 

    // from beginning
    m_elapsedTime1  = m_end - m_start;

    std::cout << "  " << m_open << " " << msg << " - " << m_name << " -- ";
    std::cout << std::fixed << std::setprecision(2) << 
       m_elapsedTime.count() << " s  // " <<
       m_elapsedTime1.count() << " s\n";

    m_open++;
  }
}
