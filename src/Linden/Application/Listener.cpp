
#include <iostream>

#include "Listener.hpp"
#include "Session.hpp"

namespace beast = boost::beast;         // from <boost/beast.hpp>
namespace http = beast::http;           // from <boost/beast/http.hpp>
namespace websocket = beast::websocket; // from <boost/beast/websocket.hpp>
namespace net = boost::asio;            // from <boost/asio.hpp>
using tcp = boost::asio::ip::tcp;       // from <boost/asio/ip/tcp.hpp>

namespace Linden::Application
{
    Listener::Listener(
        net::io_context& ioc,
        tcp::endpoint endpoint,
        Platform& platform)
        : ioc_(ioc)
        , acceptor_(ioc)
        , m_platform(platform)
    {
        beast::error_code ec;

        // Open the acceptor
        acceptor_.open(endpoint.protocol(), ec);
        if(ec)
        {
            fail(ec, "open");
            return;
        }

        // Allow address reuse
        acceptor_.set_option(net::socket_base::reuse_address(true), ec);
        if(ec)
        {
            fail(ec, "set_option");
            return;
        }

        // Bind to the server address
        acceptor_.bind(endpoint, ec);
        if(ec)
        {
            fail(ec, "bind");
            return;
        }

        // Start listening for connections
        acceptor_.listen(
            net::socket_base::max_listen_connections, ec);
        if(ec)
        {
            fail(ec, "listen");
            return;
        }
    }

    // Report a failure
    void Listener::fail(beast::error_code ec, char const* what)
    {
        std::cerr << what << ": " << ec.message() << "\n";
    }

    // Start accepting incoming connections
    void Listener::run()
    {
        do_accept();
    }

    void Listener::do_accept()
    {
        // The new connection gets its own strand
        acceptor_.async_accept(
            net::make_strand(ioc_),
            beast::bind_front_handler(
                &Listener::on_accept,
                shared_from_this()));
    }

    void Listener::on_accept(beast::error_code ec, tcp::socket socket)
    {
        if(ec)
        {
            fail(ec, "accept");
        }
        else
        {
            // Create the session and run it
            std::make_shared<Session>(std::move(socket), m_platform)->run();
        }

        // Accept another connection
        do_accept();
    }
}