#pragma once

const std::set<std::string> series({
//"wvarmax", "clfrac", "lwp", "er",
// "surf_precip"/*, "mass_dry"*/, "acc_precip",
// "cl_nc", 
"ract_com", "ract_avg", 
//"th_com", "tot_water",
"com_mom0","com_mom1","com_mom2", // higher moments need lower ones enabled!!
"ract_std_dev"
});

std::set<std::string> profs({
"00rtot", "rliq", "thl", "wvar", 
"w3rd", "prflux", "act_conc", 
"clfrac", "N_c", "non_gccn_rw_up", 
"gccn_rw_up", "non_gccn_rw_down", "gccn_rw_down", "sat_RH"
}); // rtot has to be first

