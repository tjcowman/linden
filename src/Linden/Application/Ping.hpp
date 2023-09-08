
#pragma once

#include <boost/asio.hpp>
#include <chrono>
#include <functional>
#include <iostream>
#include <string>

namespace Linden::Application
{
    class Ping
    {
    public:
        inline explicit Ping(boost::asio::io_context& ioc) :
            m_ioc(ioc),
            m_timer(m_ioc),
            m_tick(0)
        { 
            Tick();

        }

        void Tick()
        {
            ++m_tick;
            m_timer.expires_after(std::chrono::milliseconds(1000));
            m_timer.async_wait([this](const boost::system::error_code& ec){

                std::cout<<m_tick<<std::endl;
                if (m_updateHandler != nullptr)
                {
                    m_updateHandler(std::to_string(m_tick));
                }

                Tick();
            });
        }



        inline void SetUpdateHandler(const std::function<void(const std::string&)>& handler)
        {
            m_updateHandler = handler;
        }
    private:
        boost::asio::io_context& m_ioc;

        int m_tick;
        boost::asio::steady_timer m_timer;

        std::function<void(const std::string&)> m_updateHandler;
    };
}
