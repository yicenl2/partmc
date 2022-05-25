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
       coating_grid
  type(aero_data_t) :: aero_data
  type(aero_state_t) :: aero_state, aero_state_averaged
  type(env_state_t) :: env_state
  integer :: ncid, index, repeat, i_index, i_repeat, n_index, n_repeat, i_part
  real(kind=dp) :: time, del_t, tot_num_conc, tot_num_conc_avg, &
                   npart, nwet, npart_avg, nwet_avg, &
                   tot_wetmass_conc, tot_drymass_conc, &
                   tot_wetmass_conc_avg, tot_drymass_conc_avg, &
                   bulk_bc_masses, bulk_bc_masses_avg, & 
                   bulk_oc_masses, bulk_oc_masses_avg, &
                   bulk_so4_masses, bulk_so4_masses_avg, &
                   bulk_no3_masses, bulk_no3_masses_avg, &
                   bulk_nh4_masses, bulk_nh4_masses_avg, &
                   bulk_cl_masses, bulk_cl_masses_avg, &
                   bulk_h2o_masses, bulk_h2o_masses_avg, &
                   bulk_oin_masses, bulk_oin_masses_avg, &
                   bulk_na_masses, bulk_na_masses_avg, &
                   bulk_ca_masses, bulk_ca_masses_avg, &
                   bulk_co3_masses, bulk_co3_masses_avg, &
                   bulk_soa_masses, bulk_soa_masses_avg, &
                   tot_surf_area_avg, tot_surf_area_pr
  real(kind=dp) :: d_alpha_comp, d_alpha_pr, d_gamma_comp, d_gamma_pr, &
                   chi_comp, chi_pr, chi_h2o_comp, chi_h2o_pr, &
                   chi_no3_comp, chi_no3_pr, chi_so4_comp, chi_so4_pr, &
                   chi_org_comp, chi_org_pr, chi_dust_comp, chi_dust_pr
  real(kind=dp) :: aero_state_n2o5_uptake, gamma_pop
  real(kind=dp), allocatable :: gamma_part(:), gamma_core(:), gamma_coat(:)
  integer       :: n2o5_type
  character(len=PMC_UUID_LEN) :: uuid
  real(kind=dp), allocatable :: times(:), dry_diameters(:), dry_diameters_avg(:), &
       num_concs(:), num_concs_avg(:), &
       dry_masses(:), dry_masses_avg(:), &
       wet_masses(:), wet_masses_avg(:), &
       bc_masses(:),bc_masses_avg(:), bc_fracs(:), &
       so4_masses(:), so4_masses_avg(:), no3_masses(:), no3_masses_avg(:), &
       nh4_masses(:), nh4_masses_avg(:), wet_diameters(:), wet_diameters_avg(:), &
       cl_masses(:), cl_masses_avg(:), oin_masses(:), oin_masses_avg(:), &
       oc_masses(:), oc_masses_avg(:), &
       na_masses(:), na_masses_avg(:), ca_masses(:), ca_masses_avg(:), &
       co3_masses(:), co3_masses_avg(:), &
       soa_masses(:), soa_masses_avg(:), &
       wi_pr(:), wi_avg(:), &
       crit_rhs(:), scs(:), &
       crit_rhs_avg(:), scs_avg(:), &
       num_dist(:), num_dist_avg(:), &
       mass_so4_dist(:), mass_so4_dist_avg(:), mass_no3_dist(:), mass_no3_dist_avg(:), &
       diam_bc_dist_pr(:,:),diam_bc_dist_avg(:,:), &
       diam_sc_dist_pr(:,:), diam_sc_dist_avg(:,:), &
       gamma_surf_pr(:), gamma_surf_avg(:), &
       surf_area_pr(:), surf_area_avg(:), &
       surf_area_dist_pr(:), surf_area_dist_avg(:), &
       mass_dist_pr(:), mass_dist_avg(:), &
       h2o_masses(:), h2o_masses_new(:), h2o_masses_avg(:), &
       diam_gamma_dist_avg(:,:), diam_gamma_dist_pr(:,:), &
       diam_h2o_dist_avg(:,:), diam_h2o_dist_pr(:,:), &
       diam_no3_dist_pr(:,:),diam_no3_dist_avg(:,:), diam_so4_dist_pr(:,:),diam_so4_dist_avg(:,:), &
       diam_nh4_dist_pr(:,:),diam_nh4_dist_avg(:,:), diam_oin_dist_pr(:,:),diam_oin_dist_avg(:,:), &
       diam_oc_dist_pr(:,:),diam_oc_dist_avg(:,:), diam_cl_dist_pr(:,:),diam_cl_dist_avg(:,:), &
       diam_na_dist_pr(:,:),diam_na_dist_avg(:,:), diam_ca_dist_pr(:,:),diam_ca_dist_avg(:,:), &
       diam_co3_dist_pr(:,:),diam_co3_dist_avg(:,:), diam_soa_dist_pr(:,:),diam_soa_dist_avg(:,:), &
       diam_wi_dist_pr(:,:), diam_wi_dist_avg(:,:)
  type(stats_1d_t) :: stats_num_dist, stats_num_dist_avg, &
       stats_d_alpha_comp, stats_d_alpha_pr, & 
       stats_tot_num_conc, stats_tot_num_conc_avg, &
       stats_d_gamma_comp, stats_d_gamma_pr, & 
       stats_chi_pr, stats_chi_h2o_pr, &
       stats_chi_no3_pr, stats_chi_so4_pr, &
       stats_chi_org_pr, stats_chi_dust_pr, &
       stats_chi_comp, stats_chi_h2o_comp, &
       stats_chi_no3_comp, stats_chi_so4_comp, &
       stats_chi_org_comp, stats_chi_dust_comp, &
       stats_npart, stats_nwet, &
       stats_npart_avg, stats_nwet_avg, &
       stats_mass_so4_dist, stats_mass_so4_dist_avg, & 
       stats_mass_no3_dist, stats_mass_no3_dist_avg, &
       stats_h2o_masses, stats_h2o_masses_new, stats_h2o_masses_avg, &
       stats_n2o5_uptake_pr, stats_n2o5_uptake_comp, &
       stats_gamma_pop_pr, stats_gamma_pop_comp, &
       stats_surf_area_dist_pr, stats_surf_area_dist_avg, &
       stats_gamma_surf_pr, stats_gamma_surf_avg, &
       stats_mass_dist_pr, stats_mass_dist_avg, &
       stats_tot_drymass_conc, stats_tot_wetmass_conc, &
       stats_tot_drymass_conc_avg, stats_tot_wetmass_conc_avg, &
       stats_bulk_bc_masses, stats_bulk_bc_masses_avg, &
       stats_bulk_oc_masses, stats_bulk_oc_masses_avg, &
       stats_bulk_so4_masses, stats_bulk_so4_masses_avg, &
       stats_bulk_no3_masses, stats_bulk_no3_masses_avg, &
       stats_bulk_nh4_masses, stats_bulk_nh4_masses_avg, &
       stats_bulk_cl_masses, stats_bulk_cl_masses_avg, &
       stats_bulk_h2o_masses, stats_bulk_h2o_masses_avg, &
       stats_bulk_oin_masses, stats_bulk_oin_masses_avg, &
       stats_bulk_na_masses, stats_bulk_na_masses_avg, &
       stats_bulk_ca_masses, stats_bulk_ca_masses_avg, &
       stats_bulk_co3_masses, stats_bulk_co3_masses_avg, &
       stats_bulk_soa_masses, stats_bulk_soa_masses_avg, &
       stats_tot_surf_area_avg, stats_tot_surf_area_pr
  type(stats_2d_t) :: stats_diam_bc_dist_pr, stats_diam_bc_dist_avg, &
       stats_diam_sc_dist_pr, stats_diam_sc_dist_avg, &
       stats_diam_gamma_dist_avg, stats_diam_gamma_dist_pr, &
       stats_diam_h2o_dist_avg, stats_diam_h2o_dist_pr, &
       stats_diam_no3_dist_pr, stats_diam_no3_dist_avg, stats_diam_so4_dist_pr, &
       stats_diam_so4_dist_avg, stats_diam_nh4_dist_pr, stats_diam_nh4_dist_avg, &
       stats_diam_oin_dist_pr, stats_diam_oin_dist_avg, stats_diam_oc_dist_pr, &
       stats_diam_oc_dist_avg, stats_diam_cl_dist_pr, stats_diam_cl_dist_avg, &
       stats_diam_na_dist_pr, stats_diam_na_dist_avg, stats_diam_ca_dist_pr, &
       stats_diam_ca_dist_avg, stats_diam_co3_dist_pr, &
       stats_diam_co3_dist_avg, stats_diam_wi_dist_pr, &
       stats_diam_wi_dist_avg, stats_diam_soa_dist_pr, stats_diam_soa_dist_avg
  logical, allocatable :: has_gamma_comp(:), has_gamma(:)

  call pmc_mpi_init()

  call input_n_files(prefix, n_repeat, n_index)

  call bin_grid_make(diam_grid, BIN_GRID_TYPE_LOG, 180, 1d-9, 1d-3)
  call bin_grid_make(avg_bin_grid, BIN_GRID_TYPE_LOG, 1, 1d-30, 1d10)
  call bin_grid_make(bc_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(sc_grid, BIN_GRID_TYPE_LOG, 50, 1d-4, 1d0)
  call bin_grid_make(gamma_grid, BIN_GRID_TYPE_LOG, 50, 1d-3, 1d-1)
  call bin_grid_make(coating_grid, BIN_GRID_TYPE_LOG, 50, 1d-9, 1d-5)

  allocate(times(n_index))

  scs = [ real(kind=dp) :: ] ! silence compiler warnings
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

        !!!!**************************************!!!!
        !!!!!!!!!!Create an averaged aero_state!!!!!!!!!!!
        !!!!**************************************!!!!
        aero_state_averaged = aero_state
        call aero_state_bin_average_comp(aero_state_averaged, avg_bin_grid, &
            aero_data)
        !write(*,*) 'test 157'
        !!!!**************************************!!!!
        !!!!!!!!!!n2o5_type = N2O5_HYDR_COMP!!!!!!!!!!!
        !!!!**************************************!!!!
        n2o5_type = N2O5_HYDR_COMP
        call aero_n2o5_uptake(aero_state_averaged, aero_data, &
        env_state, n2o5_type, gamma_part, aero_state_n2o5_uptake, gamma_pop)

        has_gamma_comp = gamma_part > 0.0

        dry_diameters_avg = aero_state_dry_diameters(aero_state_averaged, aero_data)
        wet_diameters_avg = aero_state_diameters(aero_state_averaged, aero_data)
        num_concs_avg = aero_state_num_concs(aero_state_averaged, aero_data)
        num_dist_avg = bin_grid_histogram_1d(diam_grid, wet_diameters_avg, num_concs_avg)
        call stats_1d_add(stats_num_dist_avg, num_dist_avg)
        
        tot_num_conc_avg = sum(num_concs_avg)
        call stats_1d_add_entry(stats_tot_num_conc_avg, tot_num_conc_avg, i_index)

        wet_masses_avg = aero_state_masses(aero_state_averaged, aero_data)
        dry_masses_avg = aero_state_masses(aero_state_averaged, aero_data, &
             exclude=(/"H2O"/))
        tot_wetmass_conc_avg = sum(wet_masses_avg * num_concs_avg)
        tot_drymass_conc_avg = sum(dry_masses_avg * num_concs_avg)
        call stats_1d_add_entry(stats_tot_wetmass_conc_avg, tot_wetmass_conc_avg, i_index)
        call stats_1d_add_entry(stats_tot_drymass_conc_avg, tot_drymass_conc_avg, i_index)

        !==========Masses of different species========== 
        bc_masses_avg = aero_state_masses(aero_state_averaged, aero_data, &
             include=(/"BC"/))
        oc_masses_avg = aero_state_masses(aero_state_averaged, aero_data, &
             include=(/"OC"/))
        so4_masses_avg = aero_state_masses(aero_state_averaged, aero_data, &
             include=(/"SO4"/))
        no3_masses_avg = aero_state_masses(aero_state_averaged, aero_data, &
             include=(/"NO3"/))
        nh4_masses_avg = aero_state_masses(aero_state_averaged, aero_data, &
             include=(/"NH4"/))
        cl_masses_avg = aero_state_masses(aero_state_averaged, aero_data, &
             include=(/"Cl"/))
        h2o_masses_avg = aero_state_masses(aero_state_averaged, aero_data, &
             include=(/"H2O"/))
        oin_masses_avg = aero_state_masses(aero_state_averaged, aero_data, &
             include=(/"OIN"/))
        na_masses_avg = aero_state_masses(aero_state_averaged, aero_data, &
             include=(/"Na"/))
        ca_masses_avg = aero_state_masses(aero_state_averaged, aero_data, &
             include=(/"Ca"/))
        co3_masses_avg = aero_state_masses(aero_state_averaged, aero_data, &
             include=(/"CO3"/))
        soa_masses_avg = aero_state_masses(aero_state_averaged, aero_data, &
             include=(/"ARO1", "ARO2", "ALK1", "OLE1", "API1", "API2", "LIM1", "LIM2"/))

        !==========Bulk masses========== 
        bulk_bc_masses_avg = sum(bc_masses_avg * num_concs_avg)
        bulk_oc_masses_avg = sum(oc_masses_avg * num_concs_avg)
        bulk_so4_masses_avg = sum(so4_masses_avg * num_concs_avg)
        bulk_no3_masses_avg = sum(no3_masses_avg * num_concs_avg)
        bulk_nh4_masses_avg = sum(nh4_masses_avg * num_concs_avg)
        bulk_cl_masses_avg = sum(cl_masses_avg * num_concs_avg)
        bulk_oin_masses_avg = sum(oin_masses_avg * num_concs_avg)
        bulk_na_masses_avg = sum(na_masses_avg * num_concs_avg)
        bulk_ca_masses_avg = sum(ca_masses_avg * num_concs_avg)
        bulk_co3_masses_avg = sum(co3_masses_avg * num_concs_avg)
        bulk_soa_masses_avg = sum(soa_masses_avg * num_concs_avg)
        bulk_h2o_masses_avg = sum(h2o_masses_avg * num_concs_avg)

        call stats_1d_add_entry(stats_bulk_bc_masses_avg, bulk_bc_masses_avg,i_index)
        call stats_1d_add_entry(stats_bulk_oc_masses_avg, bulk_oc_masses_avg,i_index)
        call stats_1d_add_entry(stats_bulk_so4_masses_avg, bulk_so4_masses_avg,i_index)
        call stats_1d_add_entry(stats_bulk_no3_masses_avg, bulk_no3_masses_avg,i_index)
        call stats_1d_add_entry(stats_bulk_nh4_masses_avg, bulk_nh4_masses_avg,i_index)
        call stats_1d_add_entry(stats_bulk_cl_masses_avg, bulk_cl_masses_avg,i_index)
        call stats_1d_add_entry(stats_bulk_oin_masses_avg, bulk_oin_masses_avg,i_index)
        call stats_1d_add_entry(stats_bulk_na_masses_avg, bulk_na_masses_avg,i_index)
        call stats_1d_add_entry(stats_bulk_ca_masses_avg, bulk_ca_masses_avg,i_index)
        call stats_1d_add_entry(stats_bulk_co3_masses_avg, bulk_co3_masses_avg,i_index)
        call stats_1d_add_entry(stats_bulk_soa_masses_avg, bulk_soa_masses_avg,i_index)
        call stats_1d_add_entry(stats_bulk_h2o_masses_avg, bulk_h2o_masses_avg,i_index)

        ! Make distribution for different species
        mass_so4_dist_avg = bin_grid_histogram_1d(diam_grid, wet_diameters_avg, so4_masses_avg * num_concs_avg)
        mass_no3_dist_avg = bin_grid_histogram_1d(diam_grid, wet_diameters_avg, no3_masses_avg * num_concs_avg)

        call stats_1d_add(stats_mass_so4_dist_avg, mass_so4_dist_avg)
        call stats_1d_add(stats_mass_no3_dist_avg, mass_no3_dist_avg)
        
        !write(*,*) 'test 244'
        crit_rhs_avg = aero_state_crit_rel_humids(aero_state_averaged, aero_data, &
             env_state)
        scs_avg = crit_rhs_avg - 1d0
        diam_sc_dist_avg = bin_grid_histogram_2d(diam_grid, wet_diameters_avg, &
             sc_grid, scs_avg, num_concs_avg)
        call stats_2d_add(stats_diam_sc_dist_avg, diam_sc_dist_avg)

        call aero_state_mixing_state_metrics(aero_state_averaged, aero_data, &
             d_alpha_comp, d_gamma_comp, chi_comp, exclude=(/"H2O"/))

        call aero_state_mixing_state_metrics(aero_state_averaged, aero_data, &
             d_alpha_comp, d_gamma_comp, chi_h2o_comp, group=(/"H2O"/))

        call aero_state_mixing_state_metrics(aero_state_averaged, aero_data, &
             d_alpha_comp, d_gamma_comp, chi_no3_comp, group=(/"NO3"/))

        call aero_state_mixing_state_metrics(aero_state_averaged, aero_data, &
             d_alpha_comp, d_gamma_comp, chi_so4_comp, group=(/"SO4"/))

        call aero_state_mixing_state_metrics(aero_state_averaged, aero_data, &
             d_alpha_comp, d_gamma_comp, chi_org_comp, group=(/"OC  ","ARO1", &
             "ARO2","ALK1","OLE1","API1","API2","LIM1","LIM2"/))

        call aero_state_mixing_state_metrics(aero_state_averaged, aero_data, &
             d_alpha_comp, d_gamma_comp, chi_dust_comp, group=(/"OIN"/))

        call stats_1d_add_entry(stats_d_alpha_comp, d_alpha_comp, i_index)
        call stats_1d_add_entry(stats_d_gamma_comp, d_gamma_comp, i_index)
        call stats_1d_add_entry(stats_chi_comp, chi_comp, i_index)
        call stats_1d_add_entry(stats_chi_h2o_comp, chi_h2o_comp, i_index)
        call stats_1d_add_entry(stats_chi_no3_comp, chi_no3_comp, i_index)
        call stats_1d_add_entry(stats_chi_so4_comp, chi_so4_comp, i_index)
        call stats_1d_add_entry(stats_chi_org_comp, chi_org_comp, i_index)
        call stats_1d_add_entry(stats_chi_dust_comp, chi_dust_comp, i_index)

        !write(*,*) 'test 263'
        !==========1D histogram==========
        surf_area_avg = aero_state_surf_area_concs(aero_state_averaged, aero_data)

        tot_surf_area_avg = sum(pack(surf_area_avg, has_gamma_comp))
        call stats_1d_add_entry(stats_tot_surf_area_avg, tot_surf_area_avg, i_index)

        surf_area_dist_avg  = bin_grid_histogram_1d(diam_grid, wet_diameters_avg, &
             surf_area_avg)
        call stats_1d_add(stats_surf_area_dist_avg,  surf_area_dist_avg)

        gamma_surf_avg  = bin_grid_histogram_1d(diam_grid, wet_diameters_avg, &
             gamma_part * surf_area_avg)
        call stats_1d_add(stats_gamma_surf_avg, gamma_surf_avg)

        mass_dist_avg = bin_grid_histogram_1d(diam_grid, wet_diameters_avg, &
             wet_masses_avg * num_concs_avg)
        call stats_1d_add(stats_mass_dist_avg, mass_dist_avg)

        call stats_1d_add_entry(stats_n2o5_uptake_comp, aero_state_n2o5_uptake, i_index)
        
        call stats_1d_add_entry(stats_gamma_pop_comp, gamma_pop, i_index)

        !==========2D histogram==========
        diam_gamma_dist_avg = bin_grid_histogram_2d(diam_grid, wet_diameters_avg, &
             gamma_grid, gamma_part, num_concs_avg)
        call stats_2d_add(stats_diam_gamma_dist_avg, diam_gamma_dist_avg)        

        diam_h2o_dist_avg = bin_grid_histogram_2d(diam_grid, wet_diameters_avg, &
             bc_grid, h2o_masses_avg/wet_masses_avg, num_concs_avg)
        call stats_2d_add(stats_diam_h2o_dist_avg, diam_h2o_dist_avg)

        diam_bc_dist_avg = bin_grid_histogram_2d(diam_grid, wet_diameters_avg, &
             bc_grid, bc_masses_avg/wet_masses_avg, num_concs_avg)
        call stats_2d_add(stats_diam_bc_dist_avg, diam_bc_dist_avg)

        diam_oc_dist_avg = bin_grid_histogram_2d(diam_grid, wet_diameters_avg, &
             bc_grid, oc_masses_avg/wet_masses_avg, num_concs_avg)
        call stats_2d_add(stats_diam_oc_dist_avg, diam_oc_dist_avg)

        diam_no3_dist_avg = bin_grid_histogram_2d(diam_grid, wet_diameters_avg, &
             bc_grid, no3_masses_avg/wet_masses_avg, num_concs_avg)
        call stats_2d_add(stats_diam_no3_dist_avg, diam_no3_dist_avg)

        diam_so4_dist_avg = bin_grid_histogram_2d(diam_grid, wet_diameters_avg, &
             bc_grid, so4_masses_avg/wet_masses_avg, num_concs_avg)
        call stats_2d_add(stats_diam_so4_dist_avg, diam_so4_dist_avg)

        diam_nh4_dist_avg = bin_grid_histogram_2d(diam_grid, wet_diameters_avg, &
             bc_grid, nh4_masses_avg/wet_masses_avg, num_concs_avg)
        call stats_2d_add(stats_diam_nh4_dist_avg, diam_nh4_dist_avg)

        diam_cl_dist_avg = bin_grid_histogram_2d(diam_grid, wet_diameters_avg, &
             bc_grid, cl_masses_avg/wet_masses_avg, num_concs_avg)
        call stats_2d_add(stats_diam_cl_dist_avg, diam_cl_dist_avg)

        diam_oin_dist_avg = bin_grid_histogram_2d(diam_grid, wet_diameters_avg, &
             bc_grid, oin_masses_avg/wet_masses_avg, num_concs_avg)
        call stats_2d_add(stats_diam_oin_dist_avg, diam_oin_dist_avg)

        diam_ca_dist_avg = bin_grid_histogram_2d(diam_grid, wet_diameters_avg, &
             bc_grid, ca_masses_avg/wet_masses_avg, num_concs_avg)
        call stats_2d_add(stats_diam_ca_dist_avg, diam_ca_dist_avg)

        diam_na_dist_avg = bin_grid_histogram_2d(diam_grid, wet_diameters_avg, &
             bc_grid, na_masses_avg/wet_masses_avg, num_concs_avg)
        call stats_2d_add(stats_diam_na_dist_avg, diam_na_dist_avg)

        diam_co3_dist_avg = bin_grid_histogram_2d(diam_grid, wet_diameters_avg, &
             bc_grid, co3_masses_avg/wet_masses_avg, num_concs_avg)
        call stats_2d_add(stats_diam_co3_dist_avg, diam_co3_dist_avg)

        diam_soa_dist_avg = bin_grid_histogram_2d(diam_grid, wet_diameters_avg, &
             bc_grid, soa_masses_avg/wet_masses_avg, num_concs_avg)
        call stats_2d_add(stats_diam_soa_dist_avg, diam_soa_dist_avg)

        wi_avg = so4_masses_avg/(so4_masses_avg + no3_masses_avg)
        diam_wi_dist_avg = bin_grid_histogram_2d(diam_grid, wet_diameters_avg, &
             bc_grid, wi_avg, num_concs_avg)
        call stats_2d_add(stats_diam_wi_dist_avg, diam_wi_dist_avg)

        npart_avg = aero_state_n_part(aero_state_averaged)
        !!!!**************************************!!!!
        !!!!!!!!!!Remove dry particles!!!!!!!!!!!
        !!!!**************************************!!!!
        do i_part = aero_state_n_part(aero_state),1,-1
           if (h2o_masses_avg(i_part) == 0.0d0) then
             call aero_state_remove_particle_no_info(aero_state, &
                  i_part)
           end if
        end do

        nwet_avg = aero_state_n_part(aero_state_averaged)
        call stats_1d_add_entry(stats_npart_avg, npart_avg, i_index)
        call stats_1d_add_entry(stats_nwet_avg, nwet_avg, i_index)
 
        !write(*,*) 'test 344'
        !!!!**************************************!!!!
        !!!!!!!!!!n2o5_type = N2O5_HYDR_PR!!!!!!!!!!!
        !!!!**************************************!!!!
        !write(*,*) "test 274"
        n2o5_type = N2O5_HYDR_PR
        call aero_n2o5_uptake(aero_state, aero_data, &
        env_state, n2o5_type, gamma_part, aero_state_n2o5_uptake, gamma_pop)

        has_gamma = gamma_part > 0.0

        dry_diameters = aero_state_dry_diameters(aero_state, aero_data)
        wet_diameters = aero_state_diameters(aero_state, aero_data)
        num_concs = aero_state_num_concs(aero_state, aero_data)
        num_dist = bin_grid_histogram_1d(diam_grid, wet_diameters, num_concs)
        call stats_1d_add(stats_num_dist, num_dist)

        tot_num_conc = sum(num_concs)
        call stats_1d_add_entry(stats_tot_num_conc, tot_num_conc, i_index)

        wet_masses = aero_state_masses(aero_state, aero_data)
        dry_masses = aero_state_masses(aero_state, aero_data, &
             exclude=(/"H2O"/))
        tot_wetmass_conc = sum(wet_masses * num_concs)
        tot_drymass_conc = sum(dry_masses * num_concs)
        call stats_1d_add_entry(stats_tot_wetmass_conc, tot_wetmass_conc, i_index)
        call stats_1d_add_entry(stats_tot_drymass_conc, tot_drymass_conc, i_index)

        !==========Masses of different species==========
        bc_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"BC"/))
        oc_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"OC"/))
        so4_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"SO4"/))
        no3_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"NO3"/))
        nh4_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"NH4"/))
        cl_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"Cl"/))
        oin_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"OIN"/))
        na_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"Na"/))
        ca_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"Ca"/))
        co3_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"CO3"/))
        soa_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"ARO1", "ARO2", "ALK1", "OLE1", "API1", "API2", "LIM1", "LIM2"/))
        h2o_masses_new = aero_state_masses(aero_state, aero_data, &
             include=(/"H2O"/))

        !==========Bulk masses========== 
        bulk_bc_masses = sum(bc_masses * num_concs)
        bulk_oc_masses = sum(oc_masses * num_concs)
        bulk_so4_masses = sum(so4_masses * num_concs)
        bulk_no3_masses = sum(no3_masses * num_concs)
        bulk_nh4_masses = sum(nh4_masses * num_concs)
        bulk_cl_masses = sum(cl_masses * num_concs)
        bulk_h2o_masses = sum(h2o_masses_new * num_concs)
        bulk_oin_masses = sum(oin_masses * num_concs)
        bulk_na_masses = sum(na_masses * num_concs)
        bulk_ca_masses = sum(ca_masses * num_concs)
        bulk_co3_masses = sum(co3_masses * num_concs)
        bulk_soa_masses = sum(soa_masses * num_concs)

        call stats_1d_add_entry(stats_bulk_bc_masses, bulk_bc_masses,i_index)
        call stats_1d_add_entry(stats_bulk_oc_masses, bulk_oc_masses,i_index)
        call stats_1d_add_entry(stats_bulk_so4_masses, bulk_so4_masses,i_index)
        call stats_1d_add_entry(stats_bulk_no3_masses, bulk_no3_masses,i_index)
        call stats_1d_add_entry(stats_bulk_nh4_masses, bulk_nh4_masses,i_index)
        call stats_1d_add_entry(stats_bulk_cl_masses, bulk_cl_masses,i_index)
        call stats_1d_add_entry(stats_bulk_h2o_masses, bulk_h2o_masses,i_index)
        call stats_1d_add_entry(stats_bulk_oin_masses, bulk_oin_masses,i_index)
        call stats_1d_add_entry(stats_bulk_na_masses, bulk_na_masses,i_index)
        call stats_1d_add_entry(stats_bulk_ca_masses, bulk_ca_masses,i_index)
        call stats_1d_add_entry(stats_bulk_co3_masses, bulk_co3_masses,i_index)
        call stats_1d_add_entry(stats_bulk_soa_masses, bulk_soa_masses,i_index)        

        ! Make distribution for different species
        mass_so4_dist = bin_grid_histogram_1d(diam_grid, &
             pack(wet_diameters,has_gamma), &
             pack(so4_masses * num_concs,has_gamma))
        mass_no3_dist = bin_grid_histogram_1d(diam_grid, &
             pack(wet_diameters,has_gamma), &
             pack(no3_masses * num_concs,has_gamma))

        call stats_1d_add(stats_mass_so4_dist, mass_so4_dist)
        call stats_1d_add(stats_mass_no3_dist, mass_no3_dist)

        !write(*,*) 'test 436'
        crit_rhs = aero_state_crit_rel_humids(aero_state, aero_data, &
             env_state)
        scs = crit_rhs - 1d0
        diam_sc_dist_pr = bin_grid_histogram_2d(diam_grid, wet_diameters, &
             sc_grid, scs, num_concs)
        call stats_2d_add(stats_diam_sc_dist_pr, diam_sc_dist_pr)

        call aero_state_mixing_state_metrics(aero_state, aero_data, &
             d_alpha_pr, d_gamma_pr, chi_pr, exclude=(/"H2O"/))

        call aero_state_mixing_state_metrics(aero_state, aero_data, &
             d_alpha_pr, d_gamma_pr, chi_h2o_pr, group=(/"H2O"/))

        call aero_state_mixing_state_metrics(aero_state, aero_data, &
             d_alpha_pr, d_gamma_pr, chi_no3_pr, group=(/"NO3"/))
        
        call aero_state_mixing_state_metrics(aero_state, aero_data, &
             d_alpha_pr, d_gamma_pr, chi_so4_pr, group=(/"SO4"/))
        
        call aero_state_mixing_state_metrics(aero_state, aero_data, &
             d_alpha_pr, d_gamma_pr, chi_org_pr, group=(/"OC  ","ARO1", &
             "ARO2","ALK1","OLE1","API1","API2","LIM1","LIM2"/))
        
        call aero_state_mixing_state_metrics(aero_state, aero_data, &
             d_alpha_pr, d_gamma_pr, chi_dust_pr, group=(/"OIN"/))
        
        call stats_1d_add_entry(stats_d_alpha_pr, d_alpha_pr, i_index)
        call stats_1d_add_entry(stats_d_gamma_pr, d_gamma_pr, i_index)
        call stats_1d_add_entry(stats_chi_pr, chi_pr, i_index)
        call stats_1d_add_entry(stats_chi_h2o_pr, chi_h2o_pr, i_index)
        call stats_1d_add_entry(stats_chi_no3_pr, chi_no3_pr, i_index)
        call stats_1d_add_entry(stats_chi_so4_pr, chi_so4_pr, i_index)
        call stats_1d_add_entry(stats_chi_org_pr, chi_org_pr, i_index)
        call stats_1d_add_entry(stats_chi_dust_pr, chi_dust_pr, i_index)

        !write(*,*) 'test 455'
        !==========1D histogram========== 
        surf_area_pr = aero_state_surf_area_concs(aero_state, aero_data)

        tot_surf_area_pr = sum(pack(surf_area_pr,has_gamma))
        call stats_1d_add_entry(stats_tot_surf_area_pr, tot_surf_area_pr, i_index)

        surf_area_dist_pr = bin_grid_histogram_1d(diam_grid, &
             pack(wet_diameters,has_gamma), &
             pack(surf_area_pr,has_gamma))
        call stats_1d_add(stats_surf_area_dist_pr,  surf_area_dist_pr)

        gamma_surf_pr = bin_grid_histogram_1d(diam_grid, &
              pack(wet_diameters,has_gamma), &
              pack(gamma_part * surf_area_pr,has_gamma))
        call stats_1d_add(stats_gamma_surf_pr, gamma_surf_pr)

        mass_dist_pr = bin_grid_histogram_1d(diam_grid, &
             pack(wet_diameters,has_gamma), &
             pack(wet_masses * num_concs,has_gamma))
        call stats_1d_add(stats_mass_dist_pr, mass_dist_pr)

        call stats_1d_add_entry(stats_n2o5_uptake_pr, aero_state_n2o5_uptake, i_index)

        call stats_1d_add_entry(stats_gamma_pop_pr, gamma_pop, i_index)

        !==========2D histogram========== 
        diam_gamma_dist_pr = bin_grid_histogram_2d(diam_grid, &
             pack(wet_diameters,has_gamma), &
             gamma_grid, pack(gamma_part,has_gamma), pack(num_concs,has_gamma))
        call stats_2d_add(stats_diam_gamma_dist_pr, diam_gamma_dist_pr)

        diam_h2o_dist_pr = bin_grid_histogram_2d(diam_grid, wet_diameters, &
             bc_grid, h2o_masses_new/wet_masses, num_concs)
        call stats_2d_add(stats_diam_h2o_dist_pr, diam_h2o_dist_pr)

        diam_bc_dist_pr = bin_grid_histogram_2d(diam_grid, wet_diameters, &
             bc_grid, bc_masses/wet_masses, num_concs)
        call stats_2d_add(stats_diam_bc_dist_pr, diam_bc_dist_pr)
        
        diam_oc_dist_pr = bin_grid_histogram_2d(diam_grid, wet_diameters, &
             bc_grid, oc_masses/wet_masses, num_concs)
        call stats_2d_add(stats_diam_oc_dist_pr, diam_oc_dist_pr)

        diam_no3_dist_pr = bin_grid_histogram_2d(diam_grid, wet_diameters, &
             bc_grid, no3_masses/wet_masses, num_concs)
        call stats_2d_add(stats_diam_no3_dist_pr, diam_no3_dist_pr)

        diam_so4_dist_pr = bin_grid_histogram_2d(diam_grid, wet_diameters, &
             bc_grid, so4_masses/wet_masses, num_concs)
        call stats_2d_add(stats_diam_so4_dist_pr, diam_so4_dist_pr)

        diam_nh4_dist_pr = bin_grid_histogram_2d(diam_grid, wet_diameters, &
             bc_grid, nh4_masses/wet_masses, num_concs)
        call stats_2d_add(stats_diam_nh4_dist_pr, diam_nh4_dist_pr)

        diam_cl_dist_pr = bin_grid_histogram_2d(diam_grid, wet_diameters, &
             bc_grid, cl_masses/wet_masses, num_concs)
        call stats_2d_add(stats_diam_cl_dist_pr, diam_cl_dist_pr)

        diam_oin_dist_pr = bin_grid_histogram_2d(diam_grid, wet_diameters, &
             bc_grid, oin_masses/wet_masses, num_concs)
        call stats_2d_add(stats_diam_oin_dist_pr, diam_oin_dist_pr)

        diam_ca_dist_pr = bin_grid_histogram_2d(diam_grid, wet_diameters, &
             bc_grid, ca_masses/wet_masses, num_concs)
        call stats_2d_add(stats_diam_ca_dist_pr, diam_ca_dist_pr)
        
        diam_na_dist_pr = bin_grid_histogram_2d(diam_grid, wet_diameters, &
             bc_grid, na_masses/wet_masses, num_concs)
        call stats_2d_add(stats_diam_na_dist_pr, diam_na_dist_pr)

        diam_co3_dist_pr = bin_grid_histogram_2d(diam_grid, wet_diameters, &
             bc_grid, co3_masses/wet_masses, num_concs)
        call stats_2d_add(stats_diam_co3_dist_pr, diam_co3_dist_pr)
      
        diam_soa_dist_pr = bin_grid_histogram_2d(diam_grid, wet_diameters, &
             bc_grid, soa_masses/wet_masses, num_concs)
        call stats_2d_add(stats_diam_soa_dist_pr, diam_soa_dist_pr)
 
        wi_pr = so4_masses/(so4_masses + no3_masses)
        diam_wi_dist_pr = bin_grid_histogram_2d(diam_grid, wet_diameters, &
             bc_grid, wi_pr, num_concs)
        call stats_2d_add(stats_diam_wi_dist_pr, diam_wi_dist_pr)        
       
        npart = aero_state_n_part(aero_state)
       
        !!!!**************************************!!!!
        !!!!!!!!!!Remove dry particles!!!!!!!!!!!
        !!!!**************************************!!!!
        h2o_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"H2O"/))
        do i_part = aero_state_n_part(aero_state),1,-1
           if (h2o_masses(i_part) == 0.0d0) then
             call aero_state_remove_particle_no_info(aero_state, &
                  i_part)
           end if
        end do
 
        nwet = aero_state_n_part(aero_state)
        call stats_1d_add_entry(stats_npart, npart, i_index)
        call stats_1d_add_entry(stats_nwet, nwet, i_index)

     end do

     ! Output h2o_masses
     !open(1, file = 'h2o_masses.dat', position = 'append', status = 'unknown')
     !do i=1,aero_state_n_part(aero_state)
     !     write(1,*) wet_diameters(i), h2o_masses(i)
     !end do
     !close(1)

     !do i = 1,aero_state_n_part(aero_state)  
     !     write(*,*) wet_diameters(i), h2o_masses(i)
     !end do

     call make_filename(out_filename, prefix, "_process.nc", index)
     write(*,*) "Writing " // trim(out_filename)
     call pmc_nc_open_write(out_filename, ncid)
     call pmc_nc_write_info(ncid, uuid, "1_urban_plume process")
     call bin_grid_output_netcdf(diam_grid, ncid, "diam", unit="m")
     call bin_grid_output_netcdf(bc_grid, ncid, "bc_frac", unit="1")
     call bin_grid_output_netcdf(sc_grid, ncid, "sc", unit="1")
     call bin_grid_output_netcdf(gamma_grid, ncid, "gamma", unit="1")
     call bin_grid_output_netcdf(coating_grid, ncid, "coating", unit="1")

     call stats_1d_output_netcdf(stats_num_dist, ncid, "num_dist", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_num_dist)

     call stats_1d_output_netcdf(stats_num_dist_avg, ncid, "num_dist_avg", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_num_dist_avg)

     call stats_1d_output_netcdf(stats_mass_so4_dist, ncid, "mass_so4_dist", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_mass_so4_dist)

     call stats_1d_output_netcdf(stats_mass_so4_dist_avg, ncid, "mass_so4_dist_avg", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_mass_so4_dist_avg)

     call stats_1d_output_netcdf(stats_mass_no3_dist, ncid, "mass_no3_dist", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_mass_no3_dist)

     call stats_1d_output_netcdf(stats_mass_no3_dist_avg, ncid, "mass_no3_dist_avg", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_mass_no3_dist_avg)

     call stats_1d_output_netcdf(stats_gamma_surf_pr, ncid, "gamma_surf_pr", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_gamma_surf_pr)

     call stats_1d_output_netcdf(stats_gamma_surf_avg, ncid, "gamma_surf_avg", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_gamma_surf_avg)

     call stats_1d_output_netcdf(stats_surf_area_dist_pr, ncid, "surf_area_dist_pr", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_surf_area_dist_pr)

     call stats_1d_output_netcdf(stats_surf_area_dist_avg, ncid, "surf_area_dist_avg", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_surf_area_dist_avg)

     call stats_1d_output_netcdf(stats_mass_dist_pr, ncid, "mass_dist_pr", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_mass_dist_pr)

     call stats_1d_output_netcdf(stats_mass_dist_avg, ncid, "mass_dist_avg", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_mass_dist_avg)

     call stats_2d_output_netcdf(stats_diam_bc_dist_pr, ncid, "diam_bc_dist_pr", &
          dim_name_1="diam", dim_name_2="bc_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_bc_dist_pr)

     call stats_2d_output_netcdf(stats_diam_bc_dist_avg, ncid, "diam_bc_dist_avg", &
          dim_name_1="diam", dim_name_2="bc_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_bc_dist_avg)

     call stats_2d_output_netcdf(stats_diam_no3_dist_pr, ncid, "diam_no3_dist_pr", &
          dim_name_1="diam", dim_name_2="no3_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_no3_dist_pr)

     call stats_2d_output_netcdf(stats_diam_no3_dist_avg, ncid, "diam_no3_dist_avg", &
          dim_name_1="diam", dim_name_2="no3_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_no3_dist_avg)

     call stats_2d_output_netcdf(stats_diam_so4_dist_pr, ncid, "diam_so4_dist_pr", &
          dim_name_1="diam", dim_name_2="so4_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_so4_dist_pr)

     call stats_2d_output_netcdf(stats_diam_so4_dist_avg, ncid, "diam_so4_dist_avg", &
          dim_name_1="diam", dim_name_2="so4_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_so4_dist_avg)

     call stats_2d_output_netcdf(stats_diam_nh4_dist_pr, ncid, "diam_nh4_dist_pr", &
          dim_name_1="diam", dim_name_2="nh4_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_nh4_dist_pr)

     call stats_2d_output_netcdf(stats_diam_nh4_dist_avg, ncid, "diam_nh4_dist_avg", &
          dim_name_1="diam", dim_name_2="nh4_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_nh4_dist_avg)

     call stats_2d_output_netcdf(stats_diam_oin_dist_pr, ncid, "diam_oin_dist_pr", &
          dim_name_1="diam", dim_name_2="oin_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_oin_dist_pr)

     call stats_2d_output_netcdf(stats_diam_oin_dist_avg, ncid, "diam_oin_dist_avg", &
          dim_name_1="diam", dim_name_2="oin_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_oin_dist_avg)

     call stats_2d_output_netcdf(stats_diam_oc_dist_pr, ncid, "diam_oc_dist_pr", &
          dim_name_1="diam", dim_name_2="oc_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_oc_dist_pr)

     call stats_2d_output_netcdf(stats_diam_oc_dist_avg, ncid, "diam_oc_dist_avg", &
          dim_name_1="diam", dim_name_2="oc_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_oc_dist_avg)

     call stats_2d_output_netcdf(stats_diam_cl_dist_pr, ncid, "diam_cl_dist_pr", &
          dim_name_1="diam", dim_name_2="cl_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_cl_dist_pr)

     call stats_2d_output_netcdf(stats_diam_cl_dist_avg, ncid, "diam_cl_dist_avg", &
          dim_name_1="diam", dim_name_2="cl_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_cl_dist_avg)

     call stats_2d_output_netcdf(stats_diam_ca_dist_pr, ncid, "diam_ca_dist_pr", &
          dim_name_1="diam", dim_name_2="ca_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_ca_dist_pr)

     call stats_2d_output_netcdf(stats_diam_ca_dist_avg, ncid, "diam_ca_dist_avg", &
          dim_name_1="diam", dim_name_2="ca_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_ca_dist_avg)

     call stats_2d_output_netcdf(stats_diam_na_dist_pr, ncid, "diam_na_dist_pr", &
          dim_name_1="diam", dim_name_2="na_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_na_dist_pr)

     call stats_2d_output_netcdf(stats_diam_na_dist_avg, ncid, "diam_na_dist_avg", &
          dim_name_1="diam", dim_name_2="na_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_na_dist_avg)

     call stats_2d_output_netcdf(stats_diam_co3_dist_pr, ncid, "diam_co3_dist_pr", &
          dim_name_1="diam", dim_name_2="co3_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_co3_dist_pr)

     call stats_2d_output_netcdf(stats_diam_co3_dist_avg, ncid, "diam_co3_dist_avg", &
          dim_name_1="diam", dim_name_2="co3_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_co3_dist_avg)

     call stats_2d_output_netcdf(stats_diam_soa_dist_pr, ncid, "diam_soa_dist_pr", &
          dim_name_1="diam", dim_name_2="soa_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_soa_dist_pr)

     call stats_2d_output_netcdf(stats_diam_soa_dist_avg, ncid, "diam_soa_dist_avg", &
          dim_name_1="diam", dim_name_2="soa_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_soa_dist_avg)

     call stats_2d_output_netcdf(stats_diam_sc_dist_avg, ncid, "diam_sc_dist_avg", &
          dim_name_1="diam", dim_name_2="sc", unit="m^{-3}")
     call stats_2d_clear(stats_diam_sc_dist_avg)

     call stats_2d_output_netcdf(stats_diam_sc_dist_pr, ncid, "diam_sc_dist_pr", &
          dim_name_1="diam", dim_name_2="sc", unit="m^{-3}")
     call stats_2d_clear(stats_diam_sc_dist_pr)

     call stats_2d_output_netcdf(stats_diam_h2o_dist_avg, ncid, &
          "diam_h2o_dist_avg", dim_name_1="diam", dim_name_2="fracs", unit="m^{-3}")
     call stats_2d_clear(stats_diam_h2o_dist_avg)

     call stats_2d_output_netcdf(stats_diam_h2o_dist_pr, ncid, &
          "diam_h2o_dist_pr", dim_name_1="diam", dim_name_2="fracs", unit="m^{-3}")
     call stats_2d_clear(stats_diam_h2o_dist_pr)

     call stats_2d_output_netcdf(stats_diam_wi_dist_pr, ncid, "diam_wi_dist_pr", &
          dim_name_1="diam", dim_name_2="wi", unit="m^{-3}")
     call stats_2d_clear(stats_diam_wi_dist_pr)

     call stats_2d_output_netcdf(stats_diam_wi_dist_avg, ncid, "diam_wi_dist_avg", &
          dim_name_1="diam", dim_name_2="wi", unit="m^{-3}")
     call stats_2d_clear(stats_diam_wi_dist_avg)
       
     call stats_2d_output_netcdf(stats_diam_gamma_dist_pr, ncid, "diam_gamma_dist_pr", &
          dim_name_1="diam", dim_name_2="gamma", unit="m^{-3}")
     call stats_2d_clear(stats_diam_gamma_dist_pr)

     call stats_2d_output_netcdf(stats_diam_gamma_dist_avg, ncid, "diam_gamma_dist_avg", &
          dim_name_1="diam", dim_name_2="gamma", unit="m^{-3}")
     call stats_2d_clear(stats_diam_gamma_dist_avg)
 
     call pmc_nc_close(ncid)
  end do

  call make_filename(out_filename, prefix, "_process.nc")
  write(*,*) "Writing " // trim(out_filename)
  call pmc_nc_open_write(out_filename, ncid)
  call pmc_nc_write_info(ncid, uuid, "1_urban_plume process")
  call pmc_nc_write_real_1d(ncid, times, "time", dim_name="time", unit="s")
  call stats_1d_output_netcdf(stats_tot_num_conc, ncid, "tot_num_conc", &
       dim_name="time", unit="m^{-3}")
  call stats_1d_output_netcdf(stats_tot_num_conc_avg, ncid, "tot_num_conc_avg", &
       dim_name="time", unit="m^{-3}")
  call stats_1d_output_netcdf(stats_tot_wetmass_conc, ncid, "tot_wetmass_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_wetmass_conc_avg, ncid, "tot_wetmass_conc_avg", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_drymass_conc, ncid, "tot_drymass_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_drymass_conc_avg, ncid, "tot_drymass_conc_avg", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_bulk_bc_masses, ncid, "bulk_bc_masses", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_bulk_oc_masses, ncid, "bulk_oc_masses", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_bulk_so4_masses, ncid, "bulk_so4_masses", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_bulk_no3_masses, ncid, "bulk_no3_masses", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_bulk_nh4_masses, ncid, "bulk_nh4_masses", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_bulk_cl_masses, ncid, "bulk_cl_masses", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_bulk_h2o_masses, ncid, "bulk_h2o_masses", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_bulk_oin_masses, ncid, "bulk_oin_masses", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_bulk_na_masses, ncid, "bulk_na_masses", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_bulk_ca_masses, ncid, "bulk_ca_masses", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_bulk_co3_masses, ncid, "bulk_co3_masses", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_bulk_soa_masses, ncid, "bulk_soa_masses", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_d_alpha_comp, ncid, "d_alpha_comp", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_d_alpha_pr, ncid, "d_alpha_pr", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_d_gamma_comp, ncid, &
       "d_gamma_comp", dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_d_gamma_pr, ncid, &
       "d_gamma_pr", dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_chi_comp, ncid, "chi_comp", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_chi_pr, ncid, "chi_pr", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_chi_h2o_comp, ncid, "chi_h2o_comp", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_chi_h2o_pr, ncid, "chi_h2o_pr", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_chi_no3_comp, ncid, "chi_no3_comp", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_chi_no3_pr, ncid, "chi_no3_pr", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_chi_so4_comp, ncid, "chi_so4_comp", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_chi_so4_pr, ncid, "chi_so4_pr", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_chi_org_comp, ncid, "chi_org_comp", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_chi_org_pr, ncid, "chi_org_pr", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_chi_dust_comp, ncid, "chi_dust_comp", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_chi_dust_pr, ncid, "chi_dust_pr", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_n2o5_uptake_pr, ncid, "n2o5_uptake_pr", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_n2o5_uptake_comp, ncid, "n2o5_uptake_comp", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_gamma_pop_pr, ncid,"gamma_pop_pr", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_gamma_pop_comp, ncid,"gamma_pop_comp", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_npart, ncid, "npart", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_npart_avg, ncid, "npart_avg", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_nwet, ncid, "nwet", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_nwet_avg, ncid, "nwet_avg", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_tot_surf_area_avg, ncid, "tot_surf_area_avg", &
       dim_name="time", unit="m^{-3}")
  call stats_1d_output_netcdf(stats_tot_surf_area_pr, ncid, "tot_surf_area_pr", &
       dim_name="time", unit="m^{-3}")
  call pmc_nc_close(ncid)

  call pmc_mpi_finalize()

end program process
