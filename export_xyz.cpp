#include "export_xyz.h"
#include <fstream>

void export_xyz(std::ofstream& xyz_stream,  double t,
    const vector<Vector3d>& R,
    const vector<Vector3d>& V)
{
    int N = R.size();
    xyz_stream << N << "\n";
    xyz_stream << "t=" << t << "\n";
    for (int i = 0; i < N; i++)
    {
        Vector3d r = R[i];
        xyz_stream << "Al " << r(0) << " " << r(1) << " " << r(2) << "\n"; 
    }
}
