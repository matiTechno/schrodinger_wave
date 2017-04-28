#include <iostream>
#include <exception>
#include <MITS/app.hpp>
#include "game_scene.hpp"

int main()
{
    try
    {
        App app(600, 480, "test_scene", SDL_WINDOW_RESIZABLE | SDL_WINDOW_MAXIMIZED);
        app.start<Game_scene>();
    }
    catch(const std::exception& e)
    {
        std::cout << e.what() << std::endl;
    }
}
