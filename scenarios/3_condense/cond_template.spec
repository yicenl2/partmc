run_type particle               # Monte Carlo run
output_prefix %%OUTPUT_PREFIX%% # prefix of output files
n_loop 1                        # number of Monte Carlo loops
n_part 1                        # total number of particles
kernel brown                    # coagulation kernel
nucleate none                   # nucleation function
restart yes                     # whether to restart from saved state (yes/no)
restart_file %%RESTART_FILE%%   # saved state file to restart from

t_max 600                       # total simulation time (s)
del_t 1                         # timestep (s)
t_output 1                      # output interval (0 disables) (s)
t_progress 1                    # progress printing interval (0 disables) (s)

n_bin 160                       # number of bins
d_min 1e-10                     # minimum diameter (m)
d_max 1e-1                      # maximum diameter (m)

weight none                     # unweighted particles

gas_data gas_data.dat           # file containing gas data

aerosol_data aero_data.dat      # file containing aerosol data

temp_profile %%TEMP_PROFILE%%   # temperature profile file
height_profile height.dat       # height profile file
gas_emissions gas_emit.dat      # gas emissions file
gas_background gas_back.dat     # background gas concentrations file
aero_emissions aero_emit.dat    # aerosol emissions file
aero_background aero_back.dat   # aerosol background file

rel_humidity 0.95               # initial relative humidity (1)
pressure 1e5                    # initial pressure (Pa)
latitude 0                      # latitude (degrees, -90 to 90)
longitude 0                     # longitude (degrees, -180 to 180)
altitude 0                      # altitude (m)
start_time 21600                # start time (s since 00:00 UTC)
start_day 200                   # start day of year (UTC)

rand_init 7                     # random initialization (0 to use time)
do_coagulation no               # whether to do coagulation (yes/no)
allow_doubling no               # whether to allow doubling (yes/no)
allow_halving no                # whether to allow halving (yes/no)
do_condensation yes             # whether to do condensation (yes/no)
do_mosaic no                    # whether to do MOSAIC (yes/no)
record_removals no              # whether to record particle removals (yes/no)
do_parallel no                  # whether to run in parallel (yes/no)
