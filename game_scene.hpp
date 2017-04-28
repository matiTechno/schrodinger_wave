#include <MITS/scene_common.hpp>
#include <vector>
#include <eigen3/Eigen/Dense>

// aspect ratio

class Game_scene: public Scene
{
public:
    Game_scene();
    void update() override;
    void render() override;
    void update_coords() override;
    void render_ImGui() override;

private:
    float accumulator = 0;
    float frametime = 0.01666666f;
    //std::vector<P_data>::iterator it;
    //std::vector<P_data> data;

    //bool load_data(const std::string& filename);
    P_data data;
    int N = 2000;
    float L = 2, x_0 = -0.5f, k_0 = 80, tau = 200, sigma = L / 30.f;
    float y_proj = 1.f;
    float A = 1.0;
    float p_size = 0.005f;
    Eigen::VectorXcd psi, H, T;
    glm::vec4 p_color{1.f, 0.3f, 0.3f, 0.3f};

    void update_data();
    void update_wave_data();
};
