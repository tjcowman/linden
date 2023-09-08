
#pragma once

#include <boost/asio.hpp>
#include <boost/beast/core.hpp>
#include <boost/beast/websocket.hpp>
#include <boost/asio/dispatch.hpp>
#include <boost/asio/strand.hpp>
#include <iostream>
#include <list>
#include <memory>

#include "Listener.hpp"
#include "Session.hpp"

namespace beast = boost::beast;         // from <boost/beast.hpp>
namespace http = beast::http;           // from <boost/beast/http.hpp>
namespace websocket = beast::websocket; // from <boost/beast/websocket.hpp>
namespace net = boost::asio;            // from <boost/asio.hpp>
using tcp = boost::asio::ip::tcp;       // from <boost/asio/ip/tcp.hpp>

namespace Linden::Application
{
    // Accepts incoming connections and launches the sessions
    class Listener : public std::enable_shared_from_this<Listener>
    {
        net::io_context& ioc_;
        tcp::acceptor acceptor_;

    public:
        Listener(
            net::io_context& ioc,
            tcp::endpoint endpoint);

        // Start accepting incoming connections
        void run();

        void Update(const std::string& data);

    private:
        // Report a failure
        void fail(beast::error_code ec, char const* what);

        void do_accept();

        void on_accept(beast::error_code ec, tcp::socket socket);

        inline void SetReadHandler(std::function<void(const std::string&)>& handler)
        {
            m_readHandler = handler;
        }
        inline void SetWriteHandler(std::function<void(const std::string&)>& handler)
        {
            m_writeHandler = handler;
        }

        //! Called when the session reads from the socket.
        std::function<void(const std::string&)> m_readHandler;
        //! Called when the session writes from the socket.
        std::function<void(const std::string&)> m_writeHandler;

        // TODO does this need to be shared?
        std::list<std::shared_ptr<Session>> m_sessions;
    };
}
