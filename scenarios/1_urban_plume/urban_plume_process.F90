! Copyright (C) 2009-2013, 2016, 2017 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The process program.

!> Read NetCDF output files and process them.
program process

  use pmc_output
  use pmc_stats
  use pmc_aero_state

  character(len=PMC_MAX_FILENAME_LEN), parameter :: prefix &
       = "out/urban_plume"

  character(len=PMC_MAX_FILENAME_LEN) :: in_filename, out_filename
  type(bin_grid_t) :: diam_grid, bc_grid, sc_grid, gamma_grid, avg_bin_grid, &
      gamma_coat_grid
  type(aero_data_t) :: aero_data
  type(aero_state_t) :: aero_state, aero_state_averaged
  type(env_state_t) :: env_state
  integer :: ncid, index, repeat, i_index, i_repeat, n_index, n_repeat
  real(kind=dp) :: time, del_t, tot_num_conc, tot_mass_conc, &
                   tot_so4_masses_conc, tot_no3_masses_conc, &
                   tot_bc_masses_conc
  real(kind=dp) :: d_alpha, d_gamma, chi
  real(kind=dp) :: d_alpha_n2o5, d_gamma_n2o5, chi_n2o5
  real(kind=dp) :: aero_state_n2o5_uptake, gamma_pop, gamma_core_pop, &
       gamma_coat_pop
  real(kind=dp), allocatable :: gamma_part(:), gamma_core(:), gamma_coat(:) 
  integer       :: n2o5_type 
  character(len=PMC_UUID_LEN) :: uuid
  real(kind=dp), allocatable :: times(:), dry_diameters(:), num_concs(:), &
       dry_masses(:), dry_masses_avg(:), wet_masses(:), wet_masses_avg(:), &
       bc_masses(:),bc_masses_avg(:), bc_fracs(:), &
       so4_masses(:), no3_masses(:), nh4_masses(:), wet_diameters(:), &
       crit_rhs(:), scs(:), num_dist(:), &
       mass_bc_dist(:), mass_so4_dist(:), mass_no3_dist(:), mass_nh4_dist(:), &
       diam_bc_dist_pr(:,:),diam_bc_dist_avg(:,:), diam_sc_dist(:,:), diam_gamma_dist_pr(:,:), &
       diam_gamma_dist_avg(:,:), gamma_size_pr(:), gamma_size_avg(:), &
       surf_area_pr(:), surf_area_avg(:), surf_area_dist_pr(:), &
       surf_area_dist_avg(:), h2o_masses(:), h2o_masses_avg(:), &
       diam_h2o_2d_avg(:,:), diam_h2o_2d_pr(:,:), &
       volumes_avg(:), volumes_pr(:), volumes_core_avg(:), volumes_core_pr(:), &
       diam_coating_2d_pr(:,:), diam_coating_2d_avg(:,:), gamma_coating_dist_avg(:), &
       gamma_coating_dist_pr(:), diam_coating_2d_pr_hasgamma(:, :), &
       diam_gamma_core_2d_avg(:,:), diam_gamma_coat_2d_avg(:,:), &
       diam_gamma_core_2d_pr(:,:), diam_gamma_coat_2d_pr(:,:), &
       core_masses_pr(:), core_masses_avg(:), core_mass_dist_pr(:), &
       core_mass_dist_avg(:)
  type(stats_1d_t) :: stats_num_dist, stats_d_alpha, stats_tot_num_conc, &
       stats_tot_mass_conc, stats_d_gamma, stats_chi, &
       stats_mass_bc_dist, stats_mass_so4_dist, stats_mass_no3_dist, stats_mass_nh4_dist, &
       stats_n2o5_uptake_pr, stats_n2o5_uptake_avg, stats_gamma_pop_pr, &
       stats_gamma_pop_avg, stats_gamma_size_pr, stats_gamma_size_avg, &
       stats_surf_area_dist_pr, stats_surf_area_dist_avg, &
       stats_tot_h2o_conc, stats_tot_h2o_conc_avg, &
       stats_tot_dry_conc, stats_tot_dry_conc_avg, &
       stats_tot_surf_area_pr, stats_tot_surf_area_avg, &
       stats_gamma_coating_dist_avg, stats_gamma_coating_dist_pr, &
       stats_d_alpha_n2o5, stats_d_gamma_n2o5, stats_chi_n2o5, &
       stats_gamma_core_pop_avg, stats_gamma_coat_pop_avg, &
       stats_gamma_core_pop_pr, stats_gamma_coat_pop_pr, &
       stats_tot_core_mass_conc, stats_core_mass_dist_pr, &
       stats_core_mass_dist_avg, stats_tot_so4_masses_conc, &
       stats_tot_no3_masses_conc, stats_tot_bc_masses_conc
  type(stats_2d_t) :: stats_diam_bc_dist_pr, stats_diam_bc_dist_avg, stats_diam_sc_dist, &
       stats_diam_gamma_dist_pr, stats_diam_gamma_dist_avg, &
       stats_diam_h2o_2d_avg, stats_diam_h2o_2d_pr, stats_diam_coating_2d_pr, &
       stats_diam_coating_2d_avg, stats_diam_coating_2d_pr_hasgamma, &
       stats_diam_gamma_core_2d_avg, stats_diam_gamma_coat_2d_avg, & 
       stats_diam_gamma_core_2d_pr, stats_diam_gamma_coat_2d_pr 
  logical, allocatable :: has_gamma(:)

  call pmc_mpi_init()

  call input_n_files(prefix, n_repeat, n_index)

  call bin_grid_make(diam_grid, BIN_GRID_TYPE_LOG, 180, 1d-9, 1d-3)
  call bin_grid_make(avg_bin_grid, BIN_GRID_TYPE_LOG, 1, 1d-30, 1d10)
  call bin_grid_make(bc_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(sc_grid, BIN_GRID_TYPE_LOG, 50, 1d-4, 1d0)
  call bin_grid_make(gamma_grid, BIN_GRID_TYPE_LOG, 50, 1d-3, 1d-1)
  call bin_grid_make(gamma_coat_grid, BIN_GRID_TYPE_LOG, 50, 1d-5, 1d-2)

  allocate(times(n_index))

  scs = [ real(kind=dp) :: ] ! silence avgiler warnings
  bc_fracs = [ real(kind=dp) :: ]

  do i_index = 1,n_index
     do i_repeat = 1,n_repeat
        call make_filename(in_filename, prefix, ".nc", i_index, i_repeat)
        write(*,*) "Processing " // trim(in_filename)
        ! FIXME: add UUID check into input_state(), keyed off of index or
        ! time or something? Or init to "" and check if not this.
        call input_state(in_filename, index, time, del_t, repeat, &
             uuid, aero_data=aero_data, aero_state=aero_state, &
             env_state=env_state)

        times(i_index) = time

        bc_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"BC"/))
        num_concs = aero_state_num_concs(aero_state, aero_data)
        tot_bc_masses_conc = sum(bc_masses * num_concs)
        call stats_1d_add_entry(stats_tot_bc_masses_conc, &
             tot_bc_masses_conc,i_index)

        ! Remove dry particles
        h2o_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"H2O"/))
        do i_part = aero_state_n_part(aero_state),1,-1
                if (h2o_masses(i_part) == 0.0d0) then
                 call aero_state_remove_particle_no_info(aero_state, &
                      i_part)
                end if
        end do

        !!!!Create an averaged aero_state!!!!
        aero_state_averaged = aero_state
        call aero_state_bin_average_avg(aero_state_averaged, avg_bin_grid, &
            aero_data)

        dry_diameters = aero_state_dry_diameters(aero_state, aero_data)
        wet_diameters = aero_state_diameters(aero_state, aero_data)
        num_concs = aero_state_num_concs(aero_state, aero_data)
        num_dist = bin_grid_histogram_1d(diam_grid, dry_diameters, num_concs)
        call stats_1d_add(stats_num_dist, num_dist)

        tot_num_conc = sum(num_concs)
        call stats_1d_add_entry(stats_tot_num_conc, tot_num_conc, i_index)

        wet_masses = aero_state_masses(aero_state, aero_data)
        tot_mass_conc = sum(wet_masses * num_concs)
        call stats_1d_add_entry(stats_tot_mass_conc, tot_mass_conc, i_index)

        wet_masses_avg = aero_state_masses(aero_state_averaged, aero_data)

        dry_masses = aero_state_masses(aero_state, aero_data, &
             exclude=(/"H2O"/))
        dry_masses_avg = aero_state_masses(aero_state_averaged, aero_data, &
             exclude=(/"H2O"/))

        call stats_1d_add_entry(stats_tot_dry_conc,  sum(dry_masses * num_concs), i_index)
        call stats_1d_add_entry(stats_tot_dry_conc_avg,  sum(dry_masses_avg *num_concs), i_index)

        ! Calcluate masses of different species
        !bc_masses = aero_state_masses(aero_state, aero_data, &
        !     include=(/"BC"/))
        so4_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"SO4"/))
        no3_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"NO3"/))
        nh4_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"NH4"/))
        h2o_masses_avg = aero_state_masses(aero_state_averaged, aero_data, &
             include=(/"H2O"/))
        bc_masses_avg = aero_state_masses(aero_state_averaged, aero_data, &
             include=(/"BC"/))

        tot_so4_masses_conc = sum(so4_masses * num_concs)
        tot_no3_masses_conc = sum(no3_masses * num_concs)
        !tot_bc_masses_conc = sum(bc_masses * num_concs)
        call stats_1d_add_entry(stats_tot_so4_masses_conc, & 
             tot_so4_masses_conc,i_index)
        call stats_1d_add_entry(stats_tot_no3_masses_conc, &
             tot_no3_masses_conc,i_index)
        !call stats_1d_add_entry(stats_tot_bc_masses_conc, &
        !     tot_bc_masses_conc,i_index)
        
       core_masses_pr  = aero_state_masses(aero_state, aero_data, &
             include=(/"SO4", "NO3", "Cl ", "NH4", "CO3", &
             "Na ", "Ca ", "OIN", "BC ", "H2O"/))
       core_masses_avg  = aero_state_masses(aero_state_averaged, aero_data, &
             include=(/"SO4", "NO3", "Cl ", "NH4", "CO3", &
             "Na ", "Ca ", "OIN", "BC ", "H2O"/))
       call stats_1d_add_entry(stats_tot_core_mass_conc, sum(core_masses_pr*num_concs), i_index)
       core_mass_dist_pr   = bin_grid_histogram_1d(diam_grid, wet_diameters, core_masses_pr * num_concs)
       core_mass_dist_avg  = bin_grid_histogram_1d(diam_grid, wet_diameters, core_masses_avg * num_concs)
        call stats_1d_add(stats_core_mass_dist_pr,  core_mass_dist_pr)
        call stats_1d_add(stats_core_mass_dist_avg,  core_mass_dist_avg)

        call stats_1d_add_entry(stats_tot_h2o_conc,  sum(h2o_masses * num_concs), i_index)
        call stats_1d_add_entry(stats_tot_h2o_conc_avg,  sum(h2o_masses_avg * num_concs), i_index)


        ! Make distribution for different species
        mass_bc_dist  = bin_grid_histogram_1d(diam_grid, dry_diameters, bc_masses * num_concs)
        mass_so4_dist = bin_grid_histogram_1d(diam_grid, dry_diameters, so4_masses * num_concs)
        mass_no3_dist = bin_grid_histogram_1d(diam_grid, dry_diameters, no3_masses * num_concs)
        mass_nh4_dist = bin_grid_histogram_1d(diam_grid, dry_diameters, nh4_masses * num_concs)
        
        call stats_1d_add(stats_mass_bc_dist,  mass_bc_dist)
        call stats_1d_add(stats_mass_so4_dist, mass_so4_dist)
        call stats_1d_add(stats_mass_no3_dist, mass_no3_dist)
        call stats_1d_add(stats_mass_nh4_dist, mass_nh4_dist)
        
        crit_rhs = aero_state_crit_rel_humids(aero_state, aero_data, &
             env_state)
        scs = crit_rhs - 1d0
        diam_sc_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             sc_grid, scs, num_concs)
        call stats_2d_add(stats_diam_sc_dist, diam_sc_dist)

        call aero_state_mixing_state_metrics(aero_state, aero_data, &
             d_alpha, d_gamma, chi, exclude=(/"H2O"/))

        call stats_1d_add_entry(stats_d_alpha, d_alpha, i_index)
        call stats_1d_add_entry(stats_d_gamma, d_gamma, i_index)
        call stats_1d_add_entry(stats_chi, chi, i_index)

        call aero_state_mixing_state_metrics(aero_state, aero_data, &
             d_alpha_n2o5, d_gamma_n2o5, chi_n2o5, include=(/"SO4", &
        "NO3", "Cl ", "NH4", "CO3", "Na ", "Ca ", "OIN", "BC "/), exclude=(/"H2O"/))

        call stats_1d_add_entry(stats_d_alpha_n2o5, d_alpha_n2o5, i_index)
        call stats_1d_add_entry(stats_d_gamma_n2o5, d_gamma_n2o5, i_index)
        call stats_1d_add_entry(stats_chi_n2o5, chi_n2o5, i_index)
    
        !!!!**************************************!!!!
        !!!!!!!!!!n2o5_type = N2O5_HYDR_PR!!!!!!!!!!!
        !!!!**************************************!!!!
        n2o5_type = N2O5_HYDR_PR 
        call aero_n2o5_uptake(aero_state, aero_data, &
        env_state, n2o5_type, gamma_part, gamma_core, gamma_coat, &
        aero_state_n2o5_uptake, gamma_pop, gamma_core_pop, gamma_coat_pop)
        
        has_gamma = gamma_part > 0.0   

        surf_area_pr = aero_state_surf_area_concs(aero_state, aero_data)
        call stats_1d_add_entry(stats_tot_surf_area_pr,  &
             sum(pack(surf_area_pr, has_gamma)), i_index)
        volumes_pr = sphere_vol2rad(aero_state_volumes(aero_state, aero_data))
        volumes_core_pr = sphere_vol2rad(aero_state_volumes(aero_state, aero_data, include=(/"SO4", &
        "NO3", "Cl ", "NH4", "CO3", "Na ", "Ca ", "OIN", "BC ", "H2O"/)))
        diam_coating_2d_pr = bin_grid_histogram_2d(diam_grid, wet_diameters, &
             bc_grid, 1-volumes_core_pr/volumes_pr, num_concs)
        call stats_2d_add(stats_diam_coating_2d_pr, diam_coating_2d_pr)

        diam_coating_2d_pr_hasgamma = bin_grid_histogram_2d(diam_grid, &
             pack(wet_diameters, has_gamma), bc_grid, &
             pack(1-volumes_core_pr/volumes_pr, has_gamma), &
             pack(num_concs, has_gamma))
        call stats_2d_add(stats_diam_coating_2d_pr_hasgamma, diam_coating_2d_pr_hasgamma)
       
       gamma_coating_dist_pr = bin_grid_histogram_1d(diam_grid, &
             wet_diameters, &
             (volumes_pr - volumes_core_pr)*num_concs)
       call stats_1d_add(stats_gamma_coating_dist_pr, gamma_coating_dist_pr)

        surf_area_dist_pr  = bin_grid_histogram_1d(diam_grid, &
             pack(wet_diameters, has_gamma), &
             pack(surf_area_pr,  has_gamma))
        call stats_1d_add(stats_surf_area_dist_pr,  surf_area_dist_pr)

        gamma_size_pr  = bin_grid_histogram_1d(diam_grid, &
             pack(wet_diameters,has_gamma), &
             pack(gamma_part * surf_area_pr,has_gamma ))
        call stats_1d_add(stats_gamma_size_pr,  gamma_size_pr)


        diam_gamma_dist_pr = bin_grid_histogram_2d(diam_grid, wet_diameters, &
             gamma_grid, gamma_part, num_concs)
        call stats_2d_add(stats_diam_gamma_dist_pr, diam_gamma_dist_pr)

   
        diam_gamma_core_2d_pr = bin_grid_histogram_2d(diam_grid, &
             pack(wet_diameters,has_gamma), &
             gamma_grid, &
             pack(gamma_core, has_gamma), &
             pack(num_concs, has_gamma))
        call stats_2d_add(stats_diam_gamma_core_2d_pr, diam_gamma_core_2d_pr)

        diam_gamma_coat_2d_pr = bin_grid_histogram_2d(diam_grid, &
             pack(wet_diameters,has_gamma), &
             gamma_coat_grid, &
             pack(gamma_coat, has_gamma), & 
             pack(num_concs, has_gamma))
        call stats_2d_add(stats_diam_gamma_coat_2d_pr, diam_gamma_coat_2d_pr)         


        call stats_1d_add_entry(stats_n2o5_uptake_pr, aero_state_n2o5_uptake, i_index)
        call stats_1d_add_entry(stats_gamma_pop_pr, gamma_pop, i_index)
        call stats_1d_add_entry(stats_gamma_core_pop_pr, gamma_core_pop, i_index)
        call stats_1d_add_entry(stats_gamma_coat_pop_pr, gamma_coat_pop, i_index)

        diam_h2o_2d_pr = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             bc_grid, h2o_masses/wet_masses, num_concs)
        call stats_2d_add(stats_diam_h2o_2d_pr, diam_h2o_2d_pr)

        bc_fracs = bc_masses / dry_masses
        diam_bc_dist_pr = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             bc_grid, bc_masses/dry_masses, num_concs)
        call stats_2d_add(stats_diam_bc_dist_pr, diam_bc_dist_pr)

        !!!!**************************************!!!!
        !!!!!!!!!!n2o5_type = N2O5_HYDR_avg!!!!!!!!!!!
        !!!!**************************************!!!!
        n2o5_type = N2O5_HYDR_avg
        call aero_n2o5_uptake(aero_state, aero_data, &
        env_state, n2o5_type, gamma_part, gamma_core, gamma_coat, &
        aero_state_n2o5_uptake, gamma_pop, gamma_core_pop, gamma_coat_pop)

        surf_area_avg = aero_state_surf_area_concs(aero_state_averaged, aero_data)
        call stats_1d_add_entry(stats_tot_surf_area_avg,  sum(surf_area_avg), i_index)

        surf_area_dist_avg  = bin_grid_histogram_1d(diam_grid, wet_diameters, &
             surf_area_avg)
        call stats_1d_add(stats_surf_area_dist_avg,  surf_area_dist_avg)

        gamma_size_avg  = bin_grid_histogram_1d(diam_grid, wet_diameters, &
             gamma_part * surf_area_avg)
        call stats_1d_add(stats_gamma_size_avg, gamma_size_avg)

        diam_gamma_core_2d_avg = bin_grid_histogram_2d(diam_grid, wet_diameters, &
             gamma_grid, gamma_core, num_concs)
        call stats_2d_add(stats_diam_gamma_core_2d_avg, diam_gamma_core_2d_avg)

        diam_gamma_coat_2d_avg = bin_grid_histogram_2d(diam_grid, wet_diameters, &
             gamma_coat_grid, gamma_coat, num_concs)
        call stats_2d_add(stats_diam_gamma_coat_2d_avg, diam_gamma_coat_2d_avg)

       diam_gamma_dist_avg = bin_grid_histogram_2d(diam_grid, wet_diameters, &
             gamma_grid, gamma_part, num_concs)
        call stats_2d_add(stats_diam_gamma_dist_avg, diam_gamma_dist_avg)

        call stats_1d_add_entry(stats_n2o5_uptake_avg, aero_state_n2o5_uptake, i_index)
        call stats_1d_add_entry(stats_gamma_pop_avg, gamma_pop, i_index)
          call stats_1d_add_entry(stats_gamma_core_pop_avg, gamma_core_pop, i_index)
        call stats_1d_add_entry(stats_gamma_coat_pop_avg, gamma_coat_pop, i_index)

        diam_h2o_2d_avg = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             bc_grid, h2o_masses_avg/wet_masses_avg, num_concs)
        call stats_2d_add(stats_diam_h2o_2d_avg, diam_h2o_2d_avg)

        diam_bc_dist_avg = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             bc_grid, bc_masses_avg/dry_masses_avg, num_concs)
        call stats_2d_add(stats_diam_bc_dist_avg, diam_bc_dist_avg)

        volumes_avg = sphere_vol2rad(aero_state_volumes(aero_state_averaged, aero_data))
        volumes_core_avg = sphere_vol2rad(aero_state_volumes(aero_state_averaged, aero_data, include=(/"SO4", &
        "NO3", "Cl ", "NH4", "CO3", "Na ", "Ca ", "OIN", "BC ", "H2O"/)))
        diam_coating_2d_avg = bin_grid_histogram_2d(diam_grid, wet_diameters, &
             bc_grid, 1 - volumes_core_avg/volumes_avg, num_concs)
        call stats_2d_add(stats_diam_coating_2d_avg, diam_coating_2d_avg)

       gamma_coating_dist_avg = bin_grid_histogram_1d(diam_grid, wet_diameters, &
             (volumes_avg - volumes_core_avg)*num_concs)
       call stats_1d_add(stats_gamma_coating_dist_avg, gamma_coating_dist_avg)

     end do

     call make_filename(out_filename, prefix, "_process.nc", index)
     write(*,*) "Writing " // trim(out_filename)
     call pmc_nc_open_write(out_filename, ncid)
     call pmc_nc_write_info(ncid, uuid, "1_urban_plume process")
     call bin_grid_output_netcdf(diam_grid, ncid, "diam", unit="m")
     call bin_grid_output_netcdf(bc_grid, ncid, "bc_frac", unit="1")
     call bin_grid_output_netcdf(sc_grid, ncid, "sc", unit="1")
     call bin_grid_output_netcdf(gamma_grid, ncid, "gamma", unit="1")

     call stats_1d_output_netcdf(stats_num_dist, ncid, "num_dist", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_num_dist)

     call stats_1d_output_netcdf(stats_core_mass_dist_pr, ncid,"core_mass_dist_pr", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_core_mass_dist_pr)

      call stats_1d_output_netcdf(stats_core_mass_dist_avg, ncid,"core_mass_dist_avg", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_core_mass_dist_avg)

     call stats_1d_output_netcdf(stats_mass_bc_dist, ncid, "mass_bc_dist", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_mass_bc_dist)
     call stats_1d_output_netcdf(stats_mass_so4_dist, ncid, "mass_so4_dist", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_mass_so4_dist)
     call stats_1d_output_netcdf(stats_mass_no3_dist, ncid, "mass_no3_dist", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_mass_no3_dist)
     call stats_1d_output_netcdf(stats_mass_nh4_dist, ncid, "mass_nh4_dist", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_mass_nh4_dist)

     call stats_1d_output_netcdf(stats_gamma_size_pr, ncid, "gamma_size_pr", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_gamma_size_pr)

     call stats_1d_output_netcdf(stats_gamma_size_avg, ncid, "gamma_size_avg", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_gamma_size_avg)
 
     call stats_1d_output_netcdf(stats_surf_area_dist_pr, ncid, "surf_area_dist_pr", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_surf_area_dist_pr)

     call stats_1d_output_netcdf(stats_surf_area_dist_avg, ncid, "surf_area_dist_avg", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_surf_area_dist_avg)

     call stats_2d_output_netcdf(stats_diam_bc_dist_pr, ncid, "diam_bc_dist_pr", &
          dim_name_1="diam", dim_name_2="bc_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_bc_dist_pr)

     call stats_2d_output_netcdf(stats_diam_bc_dist_avg, ncid, "diam_bc_dist_avg", &
          dim_name_1="diam", dim_name_2="bc_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_bc_dist_avg)

     call stats_2d_output_netcdf(stats_diam_sc_dist, ncid, "diam_sc_dist", &
          dim_name_1="diam", dim_name_2="sc", unit="m^{-3}")
     call stats_2d_clear(stats_diam_sc_dist)

     call stats_2d_output_netcdf(stats_diam_gamma_dist_pr, ncid, "diam_gamma_dist_pr", &
          dim_name_1="diam", dim_name_2="gamma", unit="m^{-3}")
     call stats_2d_clear(stats_diam_gamma_dist_pr)

     call stats_2d_output_netcdf(stats_diam_gamma_core_2d_pr, ncid, &
          "diam_gamma_core_2d_pr", dim_name_1="diam", dim_name_2="gamma",unit="m^{-3}")
     call stats_2d_clear(stats_diam_gamma_core_2d_pr)

     call stats_2d_output_netcdf(stats_diam_gamma_coat_2d_pr, ncid, &
          "diam_gamma_coat_2d_pr", dim_name_1="diam", dim_name_2="gamma",unit="m^{-3}")
     call stats_2d_clear(stats_diam_gamma_coat_2d_pr)

     call stats_2d_output_netcdf(stats_diam_gamma_dist_avg, ncid, &
          "diam_gamma_dist_avg", dim_name_1="diam", dim_name_2="gamma", unit="m^{-3}")
     call stats_2d_clear(stats_diam_gamma_dist_avg)

     call stats_2d_output_netcdf(stats_diam_gamma_core_2d_avg, ncid, &
          "diam_gamma_core_2d_avg", dim_name_1="diam", dim_name_2="gamma",unit="m^{-3}")
     call stats_2d_clear(stats_diam_gamma_core_2d_avg)

      call stats_2d_output_netcdf(stats_diam_gamma_coat_2d_avg, ncid, &
          "diam_gamma_coat_2d_avg", dim_name_1="diam", dim_name_2="gamma", unit="m^{-3}")
     call stats_2d_clear(stats_diam_gamma_coat_2d_avg)
     

     call stats_2d_output_netcdf(stats_diam_h2o_2d_avg, ncid, &
          "diam_h2o_2d_avg", dim_name_1="diam", dim_name_2="fracs", unit="m^{-3}")
     call stats_2d_clear(stats_diam_h2o_2d_avg)

     call stats_2d_output_netcdf(stats_diam_h2o_2d_pr, ncid, &
          "diam_h2o_2d_pr", dim_name_1="diam", dim_name_2="fracs", unit="m^{-3}")
     call stats_2d_clear(stats_diam_h2o_2d_pr)

     call stats_2d_output_netcdf(stats_diam_coating_2d_pr, ncid, &
          "diam_coating_2d_pr", dim_name_1="diam", dim_name_2="fracs", unit="m^{-3}")
     call stats_2d_clear(stats_diam_coating_2d_pr)

      call stats_2d_output_netcdf(stats_diam_coating_2d_pr_hasgamma, ncid, &
          "diam_coating_2d_pr_hasgamma", dim_name_1="diam", dim_name_2="fracs", unit="m^{-3}")
     call stats_2d_clear(stats_diam_coating_2d_pr_hasgamma)
     
     call stats_2d_output_netcdf(stats_diam_coating_2d_avg, ncid, &
          "diam_coating_2d_avg", dim_name_1="diam", dim_name_2="fracs", unit="m^{-3}")
     call stats_2d_clear(stats_diam_coating_2d_avg)

     call stats_1d_output_netcdf(stats_gamma_coating_dist_avg, ncid, & 
          "gamma_coating_dist_avg", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_gamma_coating_dist_avg)

     call stats_1d_output_netcdf(stats_gamma_coating_dist_pr, ncid, &
          "gamma_coating_dist_pr", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_gamma_coating_dist_pr)

     call pmc_nc_close(ncid)
  end do

  call make_filename(out_filename, prefix, "_process.nc")
  write(*,*) "Writing " // trim(out_filename)
  call pmc_nc_open_write(out_filename, ncid)
  call pmc_nc_write_info(ncid, uuid, "1_urban_plume process")
  call pmc_nc_write_real_1d(ncid, times, "time", dim_name="time", unit="s")
  call stats_1d_output_netcdf(stats_tot_num_conc, ncid, "tot_num_conc", &
       dim_name="time", unit="m^{-3}")
  call stats_1d_output_netcdf(stats_tot_mass_conc, ncid, "tot_mass_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_core_mass_conc, ncid, "tot_core_mass_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_h2o_conc, ncid, "tot_h2o_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_h2o_conc_avg, ncid, "tot_h2o_conc_avg", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_dry_conc, ncid, "tot_dry_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_dry_conc_avg, ncid, "tot_dry_conc_avg",&
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_so4_masses_conc, ncid, &
       "tot_so4_masses_conc", dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_no3_masses_conc, ncid, &
       "tot_no3_masses_conc", dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_bc_masses_conc, ncid, &
       "tot_bc_masses_conc", dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_d_alpha, ncid, "d_alpha", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_d_gamma, ncid, &
       "d_gamma", dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_chi, ncid, "chi", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_d_alpha_n2o5, ncid, "d_alpha_n2o5", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_d_gamma_n2o5, ncid, &
       "d_gamma_n2o5", dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_chi_n2o5, ncid, "chi_n2o5", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_n2o5_uptake_pr, ncid, "n2o5_uptake_pr", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_n2o5_uptake_avg, ncid, "n2o5_uptake_avg", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_gamma_pop_pr, ncid,"gamma_pop_pr", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_gamma_core_pop_pr, ncid,"gamma_core_pop_pr", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_gamma_coat_pop_pr, ncid,"gamma_coat_pop_pr", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_gamma_pop_avg, ncid,"gamma_pop_avg", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_gamma_core_pop_avg, ncid,"gamma_core_pop_avg", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_gamma_coat_pop_avg, ncid,"gamma_coat_pop_avg", &
       dim_name="time", unit="1") 
  call stats_1d_output_netcdf(stats_tot_surf_area_avg, ncid,"tot_surf_area_avg", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_tot_surf_area_pr, ncid,"tot_surf_area_pr", &
       dim_name="time", unit="1")
  call pmc_nc_close(ncid)

  call pmc_mpi_finalize()

end program process
