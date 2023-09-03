
#pragma once

#include <boost/asio.hpp>
#include <boost/beast/core.hpp>
#include <boost/beast/websocket.hpp>
#include <boost/asio/dispatch.hpp>
#include <boost/asio/strand.hpp>
#include <iostream>

#include "Platform.hpp"
#include "Session.hpp"

namespace beast = boost::beast;         // from <boost/beast.hpp>
namespace http = beast::http;           // from <boost/beast/http.hpp>
namespace websocket = beast::websocket; // from <boost/beast/websocket.hpp>
namespace net = boost::asio;            // from <boost/asio.hpp>
using tcp = boost::asio::ip::tcp;       // from <boost/asio/ip/tcp.hpp>

namespace Linden::Application
{
    // Echoes back all received WebSocket messages
    class Session : public std::enable_shared_from_this<Session>
    {
        websocket::stream<beast::tcp_stream> ws_;
        beast::flat_buffer buffer_;

    public:
        // Take ownership of the socket
        explicit Session(tcp::socket&& socket, Platform& platform);

        // Get on the correct executor
        void run();

        // Start the asynchronous operation
        void on_run();

        void on_accept(beast::error_code ec);

        void do_read();

        void on_read(
            beast::error_code ec,
            std::size_t bytes_transferred);

        void on_write(
            beast::error_code ec,
            std::size_t bytes_transferred);

    private:
        // Report a failure
        void fail(beast::error_code ec, char const* what);

        Platform& m_platform;
    };
}
