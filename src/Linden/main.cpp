
#include <boost/asio.hpp>
#include <iostream>

#include "Application/Session.hpp"
#include "Application/Listener.hpp"
#include "Platform.hpp"

int main(int argc, char *argv[])
{
    auto const address = boost::beast::net::ip::make_address(argv[1]);
    auto const port = static_cast<unsigned short>(std::atoi(argv[2]));
    auto const threads = std::max<int>(1, std::atoi(argv[3]));

    // The io_context is required for all I/O
    boost::beast::net::io_context ioc{threads};

    Linden::Application::Platform platform(ioc);

    // Create and launch a listening port
    auto listener = std::make_shared<Linden::Application::Listener>(ioc, tcp::endpoint{address, port});
    
    platform.SetUpdateHandler([&](const std::string& data)
    {
        listener->Update(data);
    });

    listener->run();

    // Run the I/O service on the requested number of threads
    std::vector<std::thread> v;
    v.reserve(threads - 1);
    for(auto i = threads - 1; i > 0; --i)
        v.emplace_back(
        [&ioc]
        {
            ioc.run();
        });
    ioc.run();

    return EXIT_SUCCESS;
}
