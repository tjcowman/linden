
#include "Platform.hpp"

#include <iostream>
#include <nlohmann/json.hpp>

namespace Linden::Application
{
    void Platform::Initialize()
    {
        m_ping.SetUpdateHandler(m_updateHandler);
    }

    void Platform::Read(const std::string& data, const std::function<void(const std::string&)>& responseHandler)
    {
        auto json = nlohmann::json::parse(data, nullptr, false, false);

        if (!json.is_object())
        {
            std::cerr<<"request not an object"<<std::endl;
        }
        else
        {
            auto it = json.find("request");
            if (it == json.end())
            {
                std::cerr<<"no request specified"<<std::endl;
            }
            else
            {
                if (!it->is_string())
                {
                    std::cerr<<"request not a string"<<std::endl;
                }
                else
                {
                    auto request = it->get<std::string>();

                    if (request == "init")
                    {
                        responseHandler(RequestInit());
                    }
                }
            }
        }
    }

    std::string Platform::RequestInit() const
    {
        return "testHello";
    }
}
