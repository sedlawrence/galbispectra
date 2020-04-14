#include <cmath>
using namespace std;
int main(int argc, char** argv)
{
    float Nx = -1.3787706641;
    float Sx = 25.0;
    double r = Nx + sqrt(Sx);
    if (abs(r - 3.621229) > 0.01)
    {
        return -1;
    }
    return 0;
}
