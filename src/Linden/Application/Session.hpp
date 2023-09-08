
#pragma once

#include <boost/asio.hpp>
#include <boost/beast/core.hpp>
#include <boost/beast/websocket.hpp>
#include <boost/asio/dispatch.hpp>
#include <boost/asio/strand.hpp>
#include <functional>
#include <iostream>
#include <string>

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
        explicit Session(tcp::socket&& socket);

        // Get on the correct executor
        void run();

        void Write(const std::string& data);

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

        inline void SetReadHandler(std::function<void(const std::string&)>& handler)
        {
            m_readHandler = handler;
        }
        inline void SetWriteHandler(std::function<void(const std::string&)>& handler)
        {
            m_writeHandler = handler;
        }

    private:
        // Report a failure
        void fail(beast::error_code ec, char const* what);

        //! Called when the session reads from the socket.
        std::function<void(const std::string&)> m_readHandler;
        //! Called when the session writes from the socket.
        std::function<void(const std::string&)> m_writeHandler;
    };
}
