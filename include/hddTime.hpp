#pragma once

#include <chrono>
#include <iomanip>
#include <string>



class HddTime
{

  public:
    HddTime(std::string);
    ~HddTime();
    void CloseTime();
    void Capture();
  private:
    std::chrono::steady_clock::time_point m_start;
    std::chrono::steady_clock::time_point m_restart;
    std::chrono::steady_clock::time_point m_end;
    std::chrono::duration<double> m_elapsedTime;
    std::chrono::duration<double> m_elapsedTime1;
    int m_open;
    std::string m_name;

};
