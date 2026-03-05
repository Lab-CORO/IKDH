#pragma once

#include <hupf/Input.h>

namespace LibHUPF
{

class Calculate
{
public:
  static std::vector<std::vector<double> > InverseKin(Input inputParms, std::vector<double> &time);
};

}
