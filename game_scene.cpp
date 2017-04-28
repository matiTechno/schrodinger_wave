#include "game_scene.hpp"
#include <fftw3.h>
#include <complex>

Game_scene::Game_scene()
{
    //    int num = 1;
    //    std::string filename = "data/psi_" + std::to_string(num) + ".txt";
    //    while(load_data(filename))
    //    {
    //        ++num;
    //        filename = "data/psi_" + std::to_string(num) + ".txt";
    //    }
    //    it = data.begin();

    data.dst_alpha = GL_ONE;
    update_wave_data();
    update_data();
}

void Game_scene::update()
{
    accumulator += dt;
    while(accumulator >= frametime)
    {
        accumulator -= frametime;
        //        ++it;
        //        if(it == data.end())
        //            it = data.begin();
        update_data();
    }
}

void Game_scene::update_coords()
{
    coords.size = App::get_fb_size();

    float proj_aspect = L / y_proj;
    float c_aspect = coords.size.x / float(coords.size.y);
    if(c_aspect > proj_aspect)
        coords.size.x = coords.size.y * proj_aspect;
    else
        coords.size.y = coords.size.x / proj_aspect;

    coords.pos = App::get_fb_size() / 2 - coords.size / 2;
}

#include <glm/common.hpp>

void Game_scene::render()
{
    //renderer.load_projection(glm::vec4(-1.f, -0.2f, 2.f, 0.4f));
    //renderer.rend_particles(*it);

    renderer.load_projection(glm::vec4(-L / 2.f, -y_proj, L, y_proj * 2));
    Sprite line;
    line.position = data.vbo_data[data.vbo_data.size() / 2].pos;
    line.position.y = -10000;
    line.size.x = 0.003;
    line.size.y = 100000;
    line.color = glm::vec4(1.f, 1.f, 0.f, 1.f);
    renderer.render(line);
    renderer.rend_particles(data);
}

void Game_scene::render_ImGui()
{
    ImGui::Begin("config");
    {
        if(ImGui::Button("reset"))
            update_wave_data();

        if(ImGui::SliderFloat("amplitude", &A, 1.f, 100.f));
        if(ImGui::SliderFloat("k0", &k_0, 50.f, 1000.f))
            update_wave_data();
        if(ImGui::InputInt("num particles", &N, 0, 0, ImGuiInputTextFlags_::ImGuiInputTextFlags_EnterReturnsTrue))
            update_wave_data();
        if(ImGui::SliderFloat("tau", &tau, 300.f, 200000.f))
            update_wave_data();
         if(ImGui::SliderFloat("p_size", &p_size, 0.0001f, 1.f, "%.6f", 4.f));
         if(ImGui::SliderFloat("sigma", &sigma, L / 100.f, L / 2.f))
             update_wave_data();
         ImGui::SliderFloat4("color", reinterpret_cast<float*>(&p_color), 0.f, 1.f);
    }
    ImGui::End();
    ImGui::ShowTestWindow();
}

void Game_scene::update_data()
{
    fftw_complex *ptr = reinterpret_cast<fftw_complex *>(&psi[0]);
    fftw_plan p = fftw_plan_dft_1d(N, ptr, ptr, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan q = fftw_plan_dft_1d(N, ptr, ptr, FFTW_BACKWARD, FFTW_ESTIMATE);

    psi = psi.cwiseProduct(H);
    fftw_execute(p);
    psi = psi.cwiseProduct(T);
    fftw_execute(q);
    psi = psi.cwiseProduct(H);

    for(std::size_t i = 0; i < data.vbo_data.size(); ++i)
    {
        data.vbo_data[i].pos.y =  psi[i].real() * A;
        data.vbo_data[i].size = glm::vec2(p_size);
        data.vbo_data[i].color = p_color;
    }
}

void Game_scene::update_wave_data()
{
    psi.resize(N);
    H.resize(N);
    T.resize(N);
    data.vbo_data.resize(N);
    data.num_to_render = N;

    std::complex<double> I(0, 1);
    fftw_plan p, q;
    Eigen::VectorXcd V(N);
    int plot_number;
    double k;
    std::vector<double> x(N);
    double x_0 = this->x_0, k_0 = this->k_0, tau = this->tau, sigma = this->sigma;
    for (int i = 0; i<N; i++)
    {
        x[i] = ((double)i / (N - 1) - 0.5) * 2;
        psi[i] = (exp(I*k_0*x[i] - pow((x[i] - x_0), 2) / (4 * pow(sigma, 2))));
        V[i] = i > N/2 ? 10 : 0;
        H[i] = exp(-I*V[i] * tau / 2.0);
        if (i <= (N / 2)) k = (double)i / N;
        else k = (double)(i - N) / N;
        T(i) = (1 / (double)N)*exp(-I*(2 * 3.1415*(k / L))*(2 * 3.1415*(k / L))*(tau / 2.0));
    }
    psi.normalize();
    for(std::size_t i = 0; i < data.vbo_data.size(); ++i)
    {
        data.vbo_data[i].pos.x =  x[i];
    }
}

//#include <fstream>
//#include <sstream>

//bool Game_scene::load_data(const std::string& filename)
//{
//    std::ifstream file(filename);
//    if(!file)
//        return false;

//    P_data temp;
//    temp.dst_alpha = GL_ONE;

//    std::string line;
//    while(std::getline(file, line))
//    {
//        std::stringstream ss(line);
//        Vbo_p p;
//        p.color = glm::vec4(1.f, 0.3f, 0.3f, 0.5f);
//        p.size = glm::vec2(0.005f, 0.005f);
//        ss >> p.pos.x;
//        float dum;
//        ss >> dum;
//        ss >> p.pos.y;
//        temp.vbo_data.push_back(std::move(p));
//    }

//    temp.num_to_render = temp.vbo_data.size();
//    data.push_back(std::move(temp));
//    return true;
//}
