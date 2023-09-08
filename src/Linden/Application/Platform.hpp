
#pragma once

#include <boost/asio.hpp>
#include <filesystem>
#include <functional>
#include <string>

#include "Ping.hpp"

namespace Linden::Application
{
    class Platform
    {
    public:
        inline explicit Platform(boost::asio::io_context& ioc) :
            m_ioc(ioc),
            m_ping(ioc)
        {
            Initialize();
        }

        void Initialize();

        void Read(const std::string& data, const std::function<void(const std::string&)>& responseHandler);

        inline void SetUpdateHandler(const std::function<void(const std::string&)>& handler)
        {
            m_updateHandler = handler;
        }

    private:

        std::string RequestInit() const;

        boost::asio::io_context& m_ioc;
        Ping m_ping;

        //! Working directory for the platform.
        std::filesystem::path m_wd;

        std::function<void(const std::string&)> m_updateHandler;

    };
}
