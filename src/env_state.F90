! Copyright (C) 2005-2012 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_env_state module.

!> The env_state_t structure and associated subroutines.
module pmc_env_state

  use pmc_constants
  use pmc_util
  use pmc_spec_file
  use pmc_mpi
  use pmc_netcdf
#ifdef PMC_USE_CAMP
  use camp_camp_state
#endif
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> Current environment state.
  !!
  !! All quantities are instantaneous, describing the state at a
  !! particular instant of time. Constant data and other data not
  !! associated with the current environment state is stored in
  !! scenario_t.
  type env_state_t
     !> Temperature (K).
     real(kind=dp) :: temp
     !> Relative humidity (1).
     real(kind=dp) :: rel_humid
     !> Ambient pressure (Pa).
     real(kind=dp) :: pressure
     !> Longitude (degrees).
     real(kind=dp) :: longitude
     !> Latitude (degrees).
     real(kind=dp) :: latitude
     !> Altitude (m).
     real(kind=dp) :: altitude
#ifdef PMC_USE_WRF
     !> Height of lower edge of grid cell (m).
     real(kind=dp) :: z_min
     !> Height of upper edge of grid cell (m).
     real(kind=dp) :: z_max
     !> Inverse density (m^3 kg^{-1}).
     real(kind=dp) :: inverse_density
     !> Grid cell volume.
     real(kind=dp) :: cell_volume
     !> East-West index for grid cell.
     integer :: cell_ix
     !> North-South index for grid cell.
     integer :: cell_iy
     !> Top-Bottom index for grid cell.
     integer :: cell_iz
     !> Eddy diffusivity coefficient (m^2 s^{-2}).
     real(kind=dp) :: diff_coef
     !> Transfer probability in all directions due to advection.
     real(kind=dp), allocatable :: prob_advection(:,:,:,:)
     !> Transfer probability in all directions due to diffusion.
     real(kind=dp), allocatable :: prob_diffusion(:,:,:,:)
     !> Transfer probability in vertical direction due to turbulent diffusion.
     real(kind=dp), allocatable :: prob_vert_diffusion(:,:)
#endif
     !> Start time (s since 00:00 UTC on \c start_day).
     real(kind=dp) :: start_time
     !> Start day of year (UTC).
     integer :: start_day
     !> Time since \c start_time (s).
     real(kind=dp) :: elapsed_time
     !> Solar zenith angle (radians from zenith).
     real(kind=dp) :: solar_zenith_angle
     !> Box height (m).
     real(kind=dp) :: height
     !> Scaling coefficient for additive coagulation kernel.
     real(kind=dp) :: additive_kernel_coefficient 
  end type env_state_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> env_state += env_state_delta
  subroutine env_state_add(env_state, env_state_delta)

    !> Environment.
    type(env_state_t), intent(inout) :: env_state
    !> Increment.
    type(env_state_t), intent(in) :: env_state_delta

    env_state%temp = env_state%temp + env_state_delta%temp
    env_state%rel_humid = env_state%rel_humid + env_state_delta%rel_humid
    env_state%pressure = env_state%pressure + env_state_delta%pressure
    env_state%longitude = env_state%longitude + env_state_delta%longitude
    env_state%latitude = env_state%latitude + env_state_delta%latitude
    env_state%altitude = env_state%altitude + env_state_delta%altitude
    env_state%start_time = env_state%start_time + env_state_delta%start_time
    env_state%start_day = env_state%start_day + env_state_delta%start_day
    env_state%elapsed_time = env_state%elapsed_time &
         + env_state_delta%elapsed_time
    env_state%solar_zenith_angle = env_state%solar_zenith_angle &
         + env_state_delta%solar_zenith_angle
    env_state%height = env_state%height + env_state_delta%height

  end subroutine env_state_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> env_state *= alpha
  subroutine env_state_scale(env_state, alpha)

    !> Environment.
    type(env_state_t), intent(inout) :: env_state
    !> Scale factor.
    real(kind=dp), intent(in) :: alpha

    env_state%temp = env_state%temp * alpha
    env_state%rel_humid = env_state%rel_humid * alpha
    env_state%pressure = env_state%pressure * alpha
    env_state%longitude = env_state%longitude * alpha
    env_state%latitude = env_state%latitude * alpha
    env_state%altitude = env_state%altitude * alpha
    env_state%start_time = env_state%start_time * alpha
    env_state%start_day = nint(real(env_state%start_day, kind=dp) * alpha)
    env_state%elapsed_time = env_state%elapsed_time * alpha
    env_state%solar_zenith_angle = env_state%solar_zenith_angle * alpha
    env_state%height = env_state%height * alpha

  end subroutine env_state_scale

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Adds the given water volume to the water vapor and updates all
  !> environment quantities.
  subroutine env_state_change_water_volume(env_state, dv)

    !> Environment state to update.
    type(env_state_t), intent(inout) :: env_state
    !> Volume concentration of water added (m^3/m^3).
    real(kind=dp), intent(in) :: dv

    real(kind=dp) pmv     ! ambient water vapor pressure (Pa)
    real(kind=dp) mv      ! ambient water vapor density (kg m^{-3})
                   ! pmv and mv are related by the factor molec_weight/(R*T)
    real(kind=dp) dmv     ! change of water density (kg m^{-3})

    dmv = dv * const%water_density
    pmv = env_state_sat_vapor_pressure(env_state) * env_state%rel_humid
    mv = const%water_molec_weight / (const%univ_gas_const*env_state%temp) * pmv
    mv = mv - dmv
    if (mv < 0d0) then
       call warn_msg(980320483, "relative humidity tried to go negative")
       mv = 0d0
    end if
    env_state%rel_humid = const%univ_gas_const * env_state%temp &
         / const%water_molec_weight * mv &
         / env_state_sat_vapor_pressure(env_state)

  end subroutine env_state_change_water_volume

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the current saturation vapor pressure (Pa).
  real(kind=dp) function env_state_sat_vapor_pressure(env_state)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    env_state_sat_vapor_pressure = const%water_eq_vap_press &
         * 10d0**(7.45d0 * (env_state%temp - const%water_freeze_temp) &
         / (env_state%temp - 38d0))

  end function env_state_sat_vapor_pressure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Air density (kg m^{-3}).
  real(kind=dp) function env_state_air_den(env_state)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    env_state_air_den = const%air_molec_weight &
         * env_state_air_molar_den(env_state)

  end function env_state_air_den

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Air molar density (mol m^{-3}).
  real(kind=dp) function env_state_air_molar_den(env_state)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    env_state_air_molar_den = env_state%pressure &
         / (const%univ_gas_const * env_state%temp)

  end function env_state_air_molar_den

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Converts relative humidity (1) to water vapor mixing ratio (ppb).
  !! 
  !! Uses equation (1.10) of Seinfeld and Pandis Atmospheric Chemistry and
  !! Physics From Air Pollution to Climate Change Second Edition, 2006.
  real(kind=dp) function env_state_rel_humid_to_mix_rat(env_state)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    real(kind=dp), parameter :: t_steam = 373.15 ! steam temperature (K)
    real(kind=dp) :: a, water_vp

    a = 1.0 - t_steam / env_state%temp
    a = (((-0.1299 * a - 0.6445) * a - 1.976) * a + 13.3185) * a
    water_vp = 101325.0 * exp(a)  ! (Pa)
    env_state_rel_humid_to_mix_rat = env_state%rel_humid * water_vp * 1.0e9 &
         / env_state%pressure ! (ppb)

  end function env_state_rel_humid_to_mix_rat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Condensation \f$A\f$ parameter.
  real(kind=dp) function env_state_A(env_state)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    env_state_A = 4d0 * const%water_surf_eng * const%water_molec_weight &
         / (const%univ_gas_const * env_state%temp * const%water_density)

  end function env_state_A

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert (ppb) to (molecules m^{-3}).
  real(kind=dp) function env_state_ppb_to_conc(env_state, ppb)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Mixing ratio (ppb).
    real(kind=dp), intent(in) :: ppb

    env_state_ppb_to_conc = ppb / 1d9 * env_state_air_molar_den(env_state) &
         * const%avagadro

  end function env_state_ppb_to_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert (molecules m^{-3}) to (ppb).
  real(kind=dp) function env_state_conc_to_ppb(env_state, conc)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Concentration (molecules m^{-3}).
    real(kind=dp), intent(in) :: conc

    env_state_conc_to_ppb = conc * 1d9 / env_state_air_molar_den(env_state) &
         / const%avagadro

  end function env_state_conc_to_ppb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PMC_USE_CAMP
  !> Sets CAMP environmental variables from PartMC environmental variables.
  subroutine env_state_set_camp_env_state(env_state, camp_state)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Environment state.
    type(camp_state_t), intent(in) :: camp_state


    camp_state%env_states(1)%val%rel_humid = env_state%rel_humid
    camp_state%env_states(1)%val%temp = env_state%temp


  end subroutine env_state_set_camp_env_state
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read environment specification from a spec file.
  subroutine spec_file_read_env_state(file, env_state)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Environment data.
    type(env_state_t), intent(inout) :: env_state

    !> \page input_format_env_state Input File Format: Environment State
    !!
    !! The environment parameters are divided into those specified at
    !! the start of the simulation and then either held constant or
    !! computed for the rest of the simulation, and those parameters
    !! given as prescribed profiles for the entire simulation
    !! duration. The variables below are for the first type --- for
    !! the prescribed profiles see \ref input_format_scenario.
    !!
    !! The environment state is specified by the parameters:
    !! - \b rel_humidity (real, dimensionless): the relative humidity
    !!   (0 is completely unsaturated and 1 is fully saturated)
    !! - \b latitude (real, unit degrees_north): the latitude of the
    !!   simulation location
    !! - \b longitude (real, unit degrees_east): the longitude of the
    !!   simulation location
    !! - \b altitude (real, unit m): the altitude of the simulation
    !!   location
    !! - \b start_time (real, unit s): the time-of-day of the start of
    !!   the simulation (in seconds past midnight)
    !! - \b start_day (integer): the day-of-year of the start of the
    !!   simulation (starting from 1 on the first day of the year)
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format
    !!   - \ref output_format_env_state --- the corresponding output
    !!     format
    !!   - \ref input_format_scenario --- the prescribed profiles of
    !!     other environment data

    call spec_file_read_real(file, 'rel_humidity', env_state%rel_humid)
    call spec_file_read_real(file, 'latitude', env_state%latitude)
    call spec_file_read_real(file, 'longitude', env_state%longitude)
    call spec_file_read_real(file, 'altitude', env_state%altitude)
    call spec_file_read_real(file, 'start_time', env_state%start_time)
    call spec_file_read_integer(file, 'start_day', env_state%start_day)

  end subroutine spec_file_read_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Average val over all processes.
  subroutine env_state_mix(val)

    !> Value to average.
    type(env_state_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    type(env_state_t) :: val_avg

    call pmc_mpi_allreduce_average_real(val%temp, val_avg%temp)
    call pmc_mpi_allreduce_average_real(val%rel_humid, val_avg%rel_humid)
    call pmc_mpi_allreduce_average_real(val%pressure, val_avg%pressure)
    val%temp = val_avg%temp
    val%rel_humid = val_avg%rel_humid
    val%pressure = val_avg%pressure
#endif

  end subroutine env_state_mix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Average val over all processes, with the result only on the root
  !> process.
  subroutine env_state_reduce_avg(val)

    !> Value to average.
    type(env_state_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    type(env_state_t) :: val_avg

    call pmc_mpi_reduce_avg_real(val%temp, val_avg%temp)
    call pmc_mpi_reduce_avg_real(val%rel_humid, val_avg%rel_humid)
    call pmc_mpi_reduce_avg_real(val%pressure, val_avg%pressure)
    if (pmc_mpi_rank() == 0) then
       val%temp = val_avg%temp
       val%rel_humid = val_avg%rel_humid
       val%pressure = val_avg%pressure
    end if
#endif

  end subroutine env_state_reduce_avg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_env_state(val)

    !> Value to pack.
    type(env_state_t), intent(in) :: val

    pmc_mpi_pack_size_env_state = &
         pmc_mpi_pack_size_real(val%temp) &
         + pmc_mpi_pack_size_real(val%rel_humid) &
         + pmc_mpi_pack_size_real(val%pressure) &
         + pmc_mpi_pack_size_real(val%longitude) &
         + pmc_mpi_pack_size_real(val%latitude) &
         + pmc_mpi_pack_size_real(val%altitude) &
#ifdef PMC_USE_WRF
         + pmc_mpi_pack_size_real(val%z_min) &
         + pmc_mpi_pack_size_real(val%z_max) &
         + pmc_mpi_pack_size_real(val%inverse_density) &
         + pmc_mpi_pack_size_real(val%cell_volume) &
         + pmc_mpi_pack_size_integer(val%cell_ix) &
         + pmc_mpi_pack_size_integer(val%cell_iy) &
         + pmc_mpi_pack_size_integer(val%cell_iz) &
         + pmc_mpi_pack_size_real(val%diff_coef) &
         + pmc_mpi_pack_size_real_array_4d(val%prob_advection) &
         + pmc_mpi_pack_size_real_array_4d(val%prob_diffusion) &
         + pmc_mpi_pack_size_real_array_2d(val%prob_vert_diffusion) &
#endif
         + pmc_mpi_pack_size_real(val%start_time) &
         + pmc_mpi_pack_size_integer(val%start_day) &
         + pmc_mpi_pack_size_real(val%elapsed_time) &
         + pmc_mpi_pack_size_real(val%solar_zenith_angle) &
         + pmc_mpi_pack_size_real(val%height) &
         + pmc_mpi_pack_size_real(val%additive_kernel_coefficient)

  end function pmc_mpi_pack_size_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_env_state(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(env_state_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_real(buffer, position, val%temp)
    call pmc_mpi_pack_real(buffer, position, val%rel_humid)
    call pmc_mpi_pack_real(buffer, position, val%pressure)
    call pmc_mpi_pack_real(buffer, position, val%longitude)
    call pmc_mpi_pack_real(buffer, position, val%latitude)
    call pmc_mpi_pack_real(buffer, position, val%altitude)
#ifdef PMC_USE_WRF
    call pmc_mpi_pack_real(buffer, position, val%z_min)
    call pmc_mpi_pack_real(buffer, position, val%z_max)
    call pmc_mpi_pack_real(buffer, position, val%inverse_density)
    call pmc_mpi_pack_real(buffer, position, val%cell_volume)
    call pmc_mpi_pack_integer(buffer, position, val%cell_ix)
    call pmc_mpi_pack_integer(buffer, position, val%cell_iy)
    call pmc_mpi_pack_integer(buffer, position, val%cell_iz)
    call pmc_mpi_pack_real(buffer, position, val%diff_coef)
    call pmc_mpi_pack_real_array_4d(buffer, position, val%prob_advection)
    call pmc_mpi_pack_real_array_4d(buffer, position, val%prob_diffusion)
    call pmc_mpi_pack_real_array_2d(buffer, position, val%prob_vert_diffusion)
#endif
    call pmc_mpi_pack_real(buffer, position, val%start_time)
    call pmc_mpi_pack_integer(buffer, position, val%start_day)
    call pmc_mpi_pack_real(buffer, position, val%elapsed_time)
    call pmc_mpi_pack_real(buffer, position, val%solar_zenith_angle)
    call pmc_mpi_pack_real(buffer, position, val%height)
    call pmc_mpi_pack_real(buffer, position, val%additive_kernel_coefficient)
    call assert(464101191, &
         position - prev_position <= pmc_mpi_pack_size_env_state(val))
#endif

  end subroutine pmc_mpi_pack_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_env_state(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(env_state_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_real(buffer, position, val%temp)
    call pmc_mpi_unpack_real(buffer, position, val%rel_humid)
    call pmc_mpi_unpack_real(buffer, position, val%pressure)
    call pmc_mpi_unpack_real(buffer, position, val%longitude)
    call pmc_mpi_unpack_real(buffer, position, val%latitude)
    call pmc_mpi_unpack_real(buffer, position, val%altitude)
#ifdef PMC_USE_WRF
    call pmc_mpi_unpack_real(buffer, position, val%z_min)
    call pmc_mpi_unpack_real(buffer, position, val%z_max)
    call pmc_mpi_unpack_real(buffer, position, val%inverse_density)
    call pmc_mpi_unpack_real(buffer, position, val%cell_volume)
    call pmc_mpi_unpack_integer(buffer, position, val%cell_ix)
    call pmc_mpi_unpack_integer(buffer, position, val%cell_iy)
    call pmc_mpi_unpack_integer(buffer, position, val%cell_iz)
    call pmc_mpi_unpack_real(buffer, position, val%diff_coef)
    call pmc_mpi_unpack_real_array_4d(buffer, position, val%prob_advection)
    call pmc_mpi_unpack_real_array_4d(buffer, position, val%prob_diffusion)
    call pmc_mpi_unpack_real_array_2d(buffer, position, &
         val%prob_vert_diffusion)
#endif
    call pmc_mpi_unpack_real(buffer, position, val%start_time)
    call pmc_mpi_unpack_integer(buffer, position, val%start_day)
    call pmc_mpi_unpack_real(buffer, position, val%elapsed_time)
    call pmc_mpi_unpack_real(buffer, position, val%solar_zenith_angle)
    call pmc_mpi_unpack_real(buffer, position, val%height)
    call pmc_mpi_unpack_real(buffer, position, val%additive_kernel_coefficient)
    call assert(205696745, &
         position - prev_position <= pmc_mpi_pack_size_env_state(val))
#endif

  end subroutine pmc_mpi_unpack_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the average of val across all processes, storing the
  !> result in val_avg on the root process.
  subroutine pmc_mpi_reduce_avg_env_state(val, val_avg)

    !> Value to average.
    type(env_state_t), intent(in) :: val
    !> Result.
    type(env_state_t), intent(inout) :: val_avg

    val_avg = val
    call pmc_mpi_reduce_avg_real(val%temp, val_avg%temp)
    call pmc_mpi_reduce_avg_real(val%rel_humid, val_avg%rel_humid)
    call pmc_mpi_reduce_avg_real(val%pressure, val_avg%pressure)

  end subroutine pmc_mpi_reduce_avg_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full state.
  subroutine env_state_output_netcdf(env_state, ncid)

    !> Environment state to write.
    type(env_state_t), intent(in) :: env_state
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    !> \page output_format_env_state Output File Format: Environment State
    !!
    !! The environment state NetCDF variables are:
    !!   - \b temperature (unit K): current air temperature
    !!   - \b relative_humidity (dimensionless): current air
    !!     relative humidity (value of 1 means completely saturated)
    !!   - \b pressure (unit Pa): current air pressure
    !!   - \b longitude (unit degrees_east): longitude of simulation location
    !!   - \b latitude (unit degrees_north): latitude of simulation location
    !!   - \b altitude (unit m): altitude of simulation location
    !!   - \b start_time_of_day (unit s): time-of-day of the
    !!     simulation start measured in seconds after midnight UTC
    !!   - \b start_day_of_year: day-in-year number of the simulation start
    !!     (starting from 1 on the first day of the year)
    !!   - \b elapsed_time (unit s): elapsed time since the simulation start
    !!   - \b solar_zenith_angle (unit radians): current angle from
    !!     the zenith to the sun
    !!   - \b height (unit m): current boundary layer mixing height
    !!
    !! See also:
    !!   - \ref input_format_env_state and \ref input_format_scenario
    !!     --- the corresponding input formats

    call pmc_nc_write_real(ncid, env_state%temp, "temperature", unit="K", &
         standard_name="air_temperature")
    call pmc_nc_write_real(ncid, env_state%rel_humid, &
         "relative_humidity", unit="1", standard_name="relative_humidity")
    call pmc_nc_write_real(ncid, env_state%pressure, "pressure", unit="Pa", &
         standard_name="air_pressure")
    call pmc_nc_write_real(ncid, env_state%longitude, "longitude", &
         unit="degree_east", standard_name="longitude")
    call pmc_nc_write_real(ncid, env_state%latitude, "latitude", &
         unit="degree_north", standard_name="latitude")
    call pmc_nc_write_real(ncid, env_state%altitude, "altitude", unit="m", &
         standard_name="altitude")
#ifdef PMC_USE_WRF
    call pmc_nc_write_real(ncid, env_state%z_min, "bottom_boundary_altitude", &
         unit="m", standard_name="bottom_altitude")
    call pmc_nc_write_real(ncid, env_state%z_max, "top_boundary_altitude", &
         unit="m", standard_name="top_altitude")
    call pmc_nc_write_real(ncid, env_state%inverse_density, &
         "inverse_density", unit="m3kg-1", standard_name="inverse_density")
    call pmc_nc_write_real(ncid, env_state%cell_volume, "cell_volume", &
         unit="m3", standard_name="cell_volume")
    call pmc_nc_write_integer(ncid,env_state%cell_ix,"x_index", &
         description="east-west index of WRF domain")
    call pmc_nc_write_integer(ncid,env_state%cell_iy,"y_index", &
         description="north-south index of WRF domain")
    call pmc_nc_write_integer(ncid,env_state%cell_iz,"z_index", &
         description="top-bottom index of WRF domain")
    call pmc_nc_write_real(ncid, env_state%diff_coef, "eddy_diff", &
         unit="m2s-1", description="eddy diffusion coefficient")
#endif
    call pmc_nc_write_real(ncid, env_state%start_time, &
         "start_time_of_day", unit="s", description="time-of-day of " &
         // "simulation start in seconds since midnight")
    call pmc_nc_write_integer(ncid, env_state%start_day, &
         "start_day_of_year", &
         description="day-of-year number of simulation start")
    call pmc_nc_write_real(ncid, env_state%elapsed_time, "elapsed_time", &
         unit="s", description="elapsed time since simulation start")
    call pmc_nc_write_real(ncid, env_state%solar_zenith_angle, &
         "solar_zenith_angle", unit="radian", &
         description="current angle from the zenith to the sun")
    call pmc_nc_write_real(ncid, env_state%height, "height", unit="m", &
         long_name="boundary layer mixing height")

  end subroutine env_state_output_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read full state.
  subroutine env_state_input_netcdf(env_state, ncid)

    !> Environment state to read.
    type(env_state_t), intent(inout) :: env_state
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    call pmc_nc_read_real(ncid, env_state%temp, "temperature")
    call pmc_nc_read_real(ncid, env_state%rel_humid, "relative_humidity")
    call pmc_nc_read_real(ncid, env_state%pressure, "pressure")
    call pmc_nc_read_real(ncid, env_state%longitude, "longitude")
    call pmc_nc_read_real(ncid, env_state%latitude, "latitude")
    call pmc_nc_read_real(ncid, env_state%altitude, "altitude")
#ifdef PMC_USE_WRF
    call pmc_nc_read_real(ncid, env_state%z_min, "bottom_boundary_altitude")
    call pmc_nc_read_real(ncid, env_state%z_max, "top_boundary_altitude")
    call pmc_nc_read_real(ncid, env_state%inverse_density, "inverse_density")
    call pmc_nc_read_real(ncid, env_state%cell_volume, "cell_volume")
    call pmc_nc_read_integer(ncid,env_state%cell_ix,"x_index")
    call pmc_nc_read_integer(ncid,env_state%cell_iy,"y_index")
    call pmc_nc_read_integer(ncid,env_state%cell_iz,"z_index")
    call pmc_nc_read_real(ncid, env_state%diff_coef, "eddy_diff")
#endif
    call pmc_nc_read_real(ncid, env_state%start_time, &
         "start_time_of_day")
    call pmc_nc_read_integer(ncid, env_state%start_day, &
         "start_day_of_year")
    call pmc_nc_read_real(ncid, env_state%elapsed_time, "elapsed_time")
    call pmc_nc_read_real(ncid, env_state%solar_zenith_angle, &
         "solar_zenith_angle")
    call pmc_nc_read_real(ncid, env_state%height, "height")

  end subroutine env_state_input_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_env_state
