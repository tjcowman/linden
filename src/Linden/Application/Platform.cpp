
#include "Platform.hpp"

#include <iostream>

namespace Linden::Application
{
    void Platform::Read(const std::string& data)
    {
        std::cout<<data<<std::endl;
    }
}
