
#pragma once

#include <boost/asio.hpp>
#include <boost/beast/core.hpp>
#include <boost/beast/websocket.hpp>
#include <boost/asio/dispatch.hpp>
#include <boost/asio/strand.hpp>
#include <iostream>

#include "Listener.hpp"
#include "Platform.hpp"

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
            tcp::endpoint endpoint,
            Platform& platform);

        // Start accepting incoming connections
        void run();

    private:
        // Report a failure
        void fail(beast::error_code ec, char const* what);

        void do_accept();

        void on_accept(beast::error_code ec, tcp::socket socket);

        Platform& m_platform;
    };
}
