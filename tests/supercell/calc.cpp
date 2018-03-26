#include <cstdlib> // system()
#include <set>
#include <string>
#include <sstream> // std::ostringstream

#include "../common.hpp"

using std::ostringstream;
using std::set;
using std::string;

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting one argument - CMAKE_BINARY_DIR");

  string outdir;

  /*  outdir=std::getenv("STORAGE");
  string tmp=std::getenv("NODE_CONF");
  outdir+="/"+tmp+"/out_lgrngn";
  std::cout << "output directory: " << outdir << std::endl;
*/
  string opts_common = 
    //"--case=supercell --outfreq=150 --nt=3600 --spinup=0 --nx=65 --ny=65 --nz=52 --dt=5 --relax_th_rv=false";
    "--case=supercell --outfreq=60 --nt=1440 --spinup=0 --nx=65 --ny=65 --nz=52 --dt=5 --relax_th_rv=false";
    //"--case=supercell --outfreq=30 --nt=2880 --spinup=0 --nx=65 --ny=65 --nz=52 --dt=2.5 --relax_th_rv=false";
  set<string> opts_micro({
    "--micro=blk_1m --outdir=out_blk_1m"
    //"--micro=blk_1m --outdir=out_blk_1m --cond=false --cevp=false --revp=false --conv=false --accr=false --sedi=false"
    //"--micro=blk_1m --outdir=out_blk_1m --revp=false --conv=false --accr=false --sedi=false"
  });

  for (auto &opts_m : opts_micro)
  {
    ostringstream cmd;
    cmd << av[1] << "/src/bicycles " << opts_common << " " << opts_m;  
    notice_macro("about to call: " << cmd.str())

    if (EXIT_SUCCESS != system(cmd.str().c_str()))
      error_macro("model run failed: " << cmd.str())
  }
}
