#pragma once
#include <random>
#include "CasesCommon.hpp"

namespace setup 
{
  namespace supercell
  {
    namespace hydrostatic = libcloudphxx::common::hydrostatic;
    namespace theta_std = libcloudphxx::common::theta_std;
    namespace theta_dry = libcloudphxx::common::theta_dry;
    namespace lognormal = libcloudphxx::common::lognormal;

    template<long N, long D>
    using rat = boost::units::static_rational<N, D>;
  
    const quantity<si::pressure, real_t> p_0 = real_t(1e5) * si::pascals;
    const quantity<si::length, real_t> 
      Z    = real_t(17850) * si::metres,
      X    = real_t(128e3) * si::metres,
      Y    = real_t(128e3) * si::metres,
      z_tr = real_t(12e3) * si::metres;
    
    const quantity<si::temperature, real_t> tht_0  = 300 * si::kelvins;

    struct wk_tht_fctr
    {
      quantity<si::temperature, real_t> tht_tr = 343 * si::kelvins;
      quantity<si::temperature, real_t> T_tr   = 213 * si::kelvins;
      real_t operator()(const real_t &zz) const
      {
        using libcloudphxx::common::moist_air::c_pd;
        using libcloudphxx::common::earth::g;
        const quantity<si::length, real_t> z = zz * si::metres;
        const auto tht_below = static_cast<quantity<si::temperature, real_t>>(tht_0 + (tht_tr - tht_0) * pow<rat<5, 4>>(z / z_tr));
        //const quantity<si::dimensionless, real_t> zfac = std::pow(z.value() / z_tr.value(), 5. / 4); 
        //const auto tht_below = tht_0 + (tht_tr - tht_0) * zfac;
        const quantity<si::temperature, real_t> tht_above = tht_tr * std::exp(g<real_t>() / (c_pd<real_t>() * T_tr) * (z - z_tr));
        return z < z_tr ? tht_below / si::kelvins : tht_above / si::kelvins;
        //return tht_0.value();
      }
      BZ_DECLARE_FUNCTOR(wk_tht_fctr);
    };
    
    struct wk_RH_fctr
    {
      real_t operator()(const real_t &zz) const
      {
        const quantity<si::length, real_t> z = zz * si::metres;
        const auto RH_below = static_cast<quantity<si::dimensionless, real_t>>(1 - 0.75 * pow<rat<5, 4>>(z / z_tr));
        const quantity<si::dimensionless, real_t> RH_above = 0.25;
        return z < z_tr ? RH_below : RH_above;
        //return 0.0;
      }
      BZ_DECLARE_FUNCTOR(wk_RH_fctr);
    };
    
    struct u_fctr
    {
      real_t operator()(const real_t &zz) const
      {
        const quantity<si::length, real_t> z = zz * si::metres;
        return 25 * std::tanh(zz / 3e3) - 0.64 * 25;
        //return 0;
      }
      BZ_DECLARE_FUNCTOR(u_fctr);
    };

    template<class concurr_t>
    class supercell : public CasesCommon<concurr_t>
    {

      protected:
  
      template <class T, class U>
      void setopts_hlpr(T &params, const U &user_params)
      {
        params.outdir = user_params.outdir;
        params.outfreq = user_params.outfreq;
        params.spinup = user_params.spinup;
        params.w_src = user_params.w_src;
        params.uv_src = false;
        params.rv_src = false;
        params.th_src = false;
        params.dt = user_params.dt;
        params.nt = user_params.nt;
        params.buoyancy_wet = true;
        params.subsidence = false;
        params.friction = false;
      }

      template <class index_t>
      void intcond_hlpr(concurr_t &solver, arr_1D_t &rhod, arr_1D_t &rv_e, arr_1D_t &th_e, int rng_seed, index_t index)
      {
        using ix = typename concurr_t::solver_t::ix;
        int nz = solver.advectee().extent(ix::w);

        real_t dz = (Z / si::metres) / (nz - 1); 
      
        //arr_1D_t th_e(nz), RH_e(nz), p_e(nz), T_e(nz);
        //wk_1982_prof(th_e, RH_e, p_e, T_e, nz);

  
        solver.advectee(ix::rv) = rv_e(index); 
        solver.advectee(ix::u)  = u_fctr{}(index * dz);
        solver.advectee(ix::w)  = 0;  
       
        // absorbers
        const real_t z_abs = 15000;
        solver.vab_coefficient() = where(index * dz >= z_abs,
                                         1. / 100 * (index * dz - z_abs)/ (Z / si::metres - z_abs),
                                         0);
        //solver.vab_coefficient() = 0;

        solver.vab_relaxed_state(0) = 0;
        solver.vab_relaxed_state(ix::w) = 0;
  
        // density profile
        solver.g_factor() = rhod(index);
  
        // initial potential temperature
        solver.advectee(ix::th) = th_e(index); 
      }

      // profiles absed on Weisman & Klemp 1982
      void wk_1982_prof(arr_1D_t &th_e, arr_1D_t &RH_e, arr_1D_t &p_e, arr_1D_t &T_e, int nz)
      {
        const real_t dz = (Z / si::metres) / (nz - 1);
        const auto k = blitz::firstIndex{};

        th_e = wk_tht_fctr{}(k * dz);
        RH_e = wk_RH_fctr{}(k * dz);

        // obtaining pressure by integrating the hydrostatic equation
        using libcloudphxx::common::earth::g;
        using libcloudphxx::common::theta_std::p_1000;
        using libcloudphxx::common::moist_air::R_d_over_c_pd;
        using libcloudphxx::common::moist_air::c_pd;
        // 
        const auto p_1000_v = p_1000<real_t>().value();
        const auto g_v = g<real_t>().value();
        const auto c_pd_v = c_pd<real_t>().value();
        const auto R_d_over_c_pd_v = R_d_over_c_pd<real_t>().value();

        p_e(0) = p_0.value();
        for (int k = 1; k < nz; ++k)
        {
          auto del_th = th_e(k) - th_e(k - 1);
          // copied from eulag, why ?
          auto inv_th = del_th > 1e-4 ? std::log(th_e(k) / th_e(k - 1)) / del_th : 1. / th_e(k);
          p_e(k) = pow(pow(p_e(k - 1), R_d_over_c_pd_v) - g_v * std::pow(p_1000_v, R_d_over_c_pd_v) * inv_th * dz / c_pd_v
                       ,
                       1./ R_d_over_c_pd_v);
        }

        T_e = th_e(k) * pow(p_e(k) / p_1000_v, R_d_over_c_pd_v);
      }
  
      // calculate the initial environmental theta and rv profiles
      void env_prof(arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &th_ref, arr_1D_t &pre_ref, arr_1D_t &rhod, arr_1D_t &w_LS, arr_1D_t &hgt_fctr_vctr, arr_1D_t &hgt_fctr_sclr, int nz, const user_params_t &user_params)
      {
        blitz::firstIndex k;
        const real_t dz = (Z / si::metres) / (nz - 1);

        arr_1D_t RH_e(nz), T_e(nz);
        wk_1982_prof(th_e, RH_e, pre_ref, T_e, nz);
        
        arr_1D_t test_T(nz), test_p(nz), test_rvs(nz);

        const real_t max_rv_e = 0.014;
        using libcloudphxx::common::const_cp::r_vs;
        for (int k = 0; k < nz; ++k)
        {
          rv_e(k) = std::min(max_rv_e, static_cast<real_t>(RH_e(k) * r_vs(T_e(k) * si::kelvins, pre_ref(k) * si::pascals)));
          //rv_e(k) = 0;
        }

        using libcloudphxx::common::earth::g;
        using libcloudphxx::common::theta_std::p_1000;
        using libcloudphxx::common::moist_air::R_d_over_c_pd;
        using libcloudphxx::common::moist_air::c_pd;
        const auto g_v = g<real_t>().value();
        const auto c_pd_v = c_pd<real_t>().value();
        const auto R_d_over_c_pd_v = R_d_over_c_pd<real_t>().value();

        //// compute reference anelastic profiles
        //const real_t stab = 1.235e-5;
        //th_ref = tht_0.value() * exp(stab * k * dz);

        //auto cs_v = g_v / (c_pd_v * tht_0.value() * stab);
        //const real_t rho_0_v = 1.11;
        //rhod = rho_0_v * exp(-stab * k * dz) * pow(1 - cs_v * (1 - exp(-stab * k * dz)), 1. / R_d_over_c_pd_v - 1);
        ////rhod = rho_0_v;
        
        // theta_std env prof to theta_dry_e
        for(int k = 0; k < nz; ++k)
        {
          th_e(k) = theta_dry::std2dry<real_t>(th_e(k) * si::kelvins, quantity<si::dimensionless, real_t>(rv_e(k))) / si::kelvins;
        }

        using libcloudphxx::common::moist_air::R_d;
        const auto R_d_v = R_d<real_t>().value();
        const auto p_1000_v = p_1000<real_t>().value();
        // take reference profiles from initial wk soundings
        for (int k = 0; k < nz; ++k)
        {
          rhod(k) = std::pow(th_e(k) * std::pow(T_e(k), R_d_v / c_pd_v - 1), -c_pd_v / R_d_v) * p_1000_v / R_d_v;
          th_ref(k) = th_e(k);
        }
        
        for (int k = 0; k < nz; ++k)
        {
          test_T(k) = libcloudphxx::common::theta_dry::T<real_t>(th_e(k) * si::kelvins,
                                                                 rhod(k) * si::kilograms / si::cubic_meters) / si::kelvins;
          test_p(k) = libcloudphxx::common::theta_dry::p<real_t>(rhod(k) * si::kilograms / si::cubic_metres,
                                                                 rv_e(k), test_T(k) * si::kelvins) / si::pascals;
          test_rvs(k) = std::min(max_rv_e, static_cast<real_t>(RH_e(k) * r_vs(test_T(k) * si::kelvins, test_p(k) * si::pascals)));
        }
        
        std::cout << "Weisman & Klemp profiles" << std::endl;
        for (int k = 0; k < nz; ++k)
        {
          std::cout << std::setw(10) << k * dz << ' '
                    << std::setw(10) << th_e(k) << ' ' 
                    << std::setw(10) << RH_e(k) << ' ' 
                    << std::setw(10) << pre_ref(k) << ' ' 
                    << std::setw(10) << T_e(k) << std::endl;
        }
        
        std::cout << "env profiles" << std::endl;
        for (int k = 0; k < nz; ++k)
        {
          std::cout << std::setw(10) << k * dz << ' '
                    << std::setw(10) << th_e(k) << ' ' 
                    << std::setw(10) << rhod(k) << ' ' 
                    << std::setw(10) << rv_e(k) << std::endl;
        }
        
        std::cout << "test profiles" << std::endl;
        for (int k = 0; k < nz; ++k)
        {
          std::cout << std::setw(10) << k * dz << ' '
                    << std::setw(10) << test_T(k) << ' ' 
                    << std::setw(10) << T_e(k) << ' ' 
                    << std::setw(10) << test_p(k) << ' ' 
                    << std::setw(10) << pre_ref(k) << ' ' 
                    << std::setw(10) << test_rvs(k) << ' ' 
                    << std::setw(10) << rv_e(k)
                    << std::endl;
        }
        
        //for (int k = 0; k < nz; ++k)
        //{
        //  T_e(k)  = test_T(k);
        //  pre_ref(k)  = test_p(k);
        //  rv_e(k)  = test_rvs(k);
        //}
      }

      // ctor
      supercell()
      {
      }
    };

    template<class concurr_t>
    class supercell_3d : public supercell<concurr_t>
    {
      void setopts(typename concurr_t::solver_t::rt_params_t &params, int nx, int ny, int nz, const user_params_t &user_params)
      {
        this->setopts_hlpr(params, user_params);
        params.di = (X / si::metres) / (nx - 1); 
        params.dj = (Y / si::metres) / (ny - 1);
        params.dk = (Z / si::metres) / (nz - 1);
        params.dz = params.dk;
      }

      void intcond(concurr_t &solver, arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, int rng_seed)
      {

        blitz::firstIndex i;
        blitz::secondIndex j;
        blitz::thirdIndex k;

        this->intcond_hlpr(solver, rhod, rv_e, th_e, rng_seed, k);
        using ix = typename concurr_t::solver_t::ix;
        //this->make_cyclic(solver.advectee(ix::th));
  
        int nx = solver.advectee().extent(ix::u);
        int ny = solver.advectee().extent(ix::v);
        int nz = solver.advectee().extent(ix::w);
        real_t dx = (X / si::metres) / (nx - 1); 
        real_t dy = (Y / si::metres) / (ny - 1); 
        real_t dz = (Z / si::metres) / (nz - 1); 

        const real_t r0x = 10e3;
        const real_t r0y = r0x;
        const real_t r0z = 1.4e3;
        
        const real_t x0 = 0.5 * (X / si::metres);
        const real_t y0 = 0.5 * (Y / si::metres);
        const real_t z0 = r0z;

        const real_t pi = std::acos(-1.0);
  
        decltype(solver.advectee()) rr(solver.advectee().shape()), delta(solver.advectee().shape());
        rr = sqrt(pow2((i * dx - x0) / r0x) + pow2((j * dy - y0) / r0y) + pow2((k * dz - z0) / r0z)); 
        

        delta = where(rr(i, j, k) <= 1.0,
                      2 * pow2(cos(0.5 * pi * rr(i, j, k))),
                      0.0);
        
        arr_1D_t wk_th(nz), wk_RH(nz), wk_p(nz), wk_T(nz);
        this->wk_1982_prof(wk_th, wk_RH, wk_p, wk_T, nz);

        for (int k = 0; k < nz; ++k)
        {
          const auto i_r = blitz::Range(0, nx - 1);
          const auto j_r = blitz::Range(0, ny - 1);
          delta(i_r, j_r, k) *= th_e(k) / wk_th(k);
        }
        
        solver.advectee(ix::th) += delta(i, j, k),

        solver.advectee(ix::v) = 0;
        solver.vab_relaxed_state(1) = 0;
      }
    };
  };
};
