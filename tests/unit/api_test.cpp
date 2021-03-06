#include <cstdlib> // system()
#include <vector>
#include <string>
#include <sstream> // std::ostringstream

#include "../common.hpp"

using std::ostringstream;
using std::vector;
using std::string;

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting one argument - CMAKE_BINARY_DIR");

  string outdir;
  string opts_common = 
    "--outfreq=1000 --nt=2 --spinup=1 --dt=1 --serial=true --prs_tol=1e-3"; 
  vector<string> opts_dim({
    "--nx=4 --nz=4",
    "--nx=4 --ny=4 --nz=4"
  });
  string opts_dim_lgrngn_3d =  "--nx=40 --ny=40 --nz=40"; // = 4 caused multiplicity overflows in the Lagrangian 3D
  vector<string> opts_micro({
    "--micro=blk_1m --outdir=out"  ,
    "--async=false --micro=lgrngn --outdir=out --backend=serial --sd_conc=8 --z_rlx_sclr=100"
  });
  vector<string> opts_case({
    "--case=moist_thermal",
    "--case=dry_thermal --cond=0 --coal=0",
    "--case=dycoms"
  });
  vector<string> opts_piggy({
    "--piggy=0",
    "--piggy=0 --save_vel=1",
    "--piggy=1 --vel_in=out/velocity_out.dat"  // take vel file from blk, cause it's ran first
  });
  string opts_pig_from_lgrngn = "--piggy=1 --vel_in=out/velocity_out.dat";

  for (auto &opts_d : opts_dim)
    for (auto &opts_m : opts_micro)
      for (auto &opts_c : opts_case)
        for (auto &opts_p : opts_piggy) // piggy has to be last to prevent overwriting of vel_out
        {
          if((opts_c == opts_case[1]) && opts_d == opts_dim[1])
          {
            std::cout << "skipping 3d dry thermal tests" << std::endl;
            continue; 
          }

          // more cells in 3d lgrngn to avoid n overflow
          if(opts_d == opts_dim[1] && opts_m == opts_micro[1])
          {
            opts_d = opts_dim_lgrngn_3d;
            if(opts_p == opts_piggy[2])
              opts_p = opts_pig_from_lgrngn;
          }
          ostringstream cmd;
          cmd << av[1] << "/src/bicycles " << opts_common << " " << opts_m << " " << opts_d << " " << opts_c << " " << opts_p;
          notice_macro("about to call: " << cmd.str())
  
          if (EXIT_SUCCESS != system(cmd.str().c_str()))
            error_macro("model run failed: " << cmd.str())
        }
}
