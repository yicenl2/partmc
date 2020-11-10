import numpy as np
import os
import shutil

if not os.path.exists("spec"):
    os.mkdir("spec")
if not os.path.exists("temp"):
    os.mkdir("temp")

values = np.loadtxt('lhs.txt')
n_scenarios = (values.shape[0])
n_modes = int((values.shape[1]) / 3)
print(n_scenarios,n_modes)
n_bin = 20

# Current in PartMC, the radius grid is:
#   call bin_grid_make(aero_state%bin_grid, BIN_GRID_TYPE_LOG, &
#               n_bin=20, min=.5d-9, max=25.0d-6)
#
bin_edges_diam = np.logspace(-9,np.log10(5e-5),n_bin+1)
np.set_printoptions(precision=3)
print(bin_edges_diam)
for counter in range(n_scenarios):
    print("counter: %d" % counter)
    filename_in = "aero_init_dist_template.dat"
    directory = "./inputs/scenario_%04i" % counter
    if not os.path.exists(directory):
        os.makedirs(directory)
    filename_out = "%s/aero_init_dist.dat" % (directory)
    print("filename_out: %s" % filename_out)
    f_in = open(filename_in, 'r')
    f_out = open(filename_out, 'w')
    d_min = np.zeros(3)
    d_max = np.zeros(3)
    for line in f_in:
        for i_mode in range(n_modes):
            line = line.replace('NUM_CONC_%1i' % (i_mode+1), 
                 "%e" %(values[counter,2*n_modes + i_mode]))
            line = line.replace('GEOM_MEAN_DIAM_%1i' % (i_mode+1),
                 "%e" %(values[counter,i_mode]))
            line = line.replace('LOG10_GEO_STD_DEV_%1i' % (i_mode+1),
                 "%f" %(np.log10(values[counter,n_modes+i_mode])))
            # Standard deviation diameter range
            num_dev = 3.0
            dev = np.log(num_dev*values[counter,n_modes+i_mode])
            d_min[i_mode] = np.exp(np.log(values[counter,i_mode]) - dev)
            d_max[i_mode] = np.exp(np.log(values[counter,i_mode]) + dev)
        f_out.write(line)

    f_in.close()
    f_out.close()
    # Make the distribution of emissions
    filename_in = "aero_emit_size_dist_template.dat"
    filename_out = "%s/aero_emit_size_dist.dat" % (directory)
    f_in = open(filename_in, 'r')
    f_out = open(filename_out, 'w')
    number_conc = np.zeros(n_bin)
    # These simulations have uniform emissions
    if (np.random.random() < 2.0/3.0):
       n_i = np.sum(values[counter,6:]) * 1e-2
       for i_bin in range(n_bin):
           number_conc[i_bin] = n_i / n_bin
    # These simulations have emissions for modes present
    else:
       all_bins = np.zeros(20)
       for i_mode in range(n_modes):
          all_bins[np.where(((bin_edges_diam > d_min[i_mode]) &  \
                            (bin_edges_diam < d_max[i_mode])))] += 1
       n_i = np.sum(values[counter,6:]) * 1e-2
       for i_bin in range(n_bin):
           if (all_bins[i_bin] > 0):
               number_conc[i_bin] = n_i / len(np.where(all_bins > 0))
    for line in f_in:
        line = line.replace('BIN_EDGES', str(bin_edges_diam).replace("\n", "")[1:-1])
        line = line.replace('NUMBER', str(number_conc).replace("\n", "")[1:-1])
        f_out.write(line)

    f_in.close()
    f_out.close()

    # Copy files
    input_files = ["gas_back.dat", "temp.dat", "aero_back.dat",
                   "aero_emit.dat", "height.dat", "aero_back_comp.dat",
                   "aero_emit_comp.dat", "aero_data.dat",
                   "aero_back_dist.dat",
                   "gas_emit.dat", "pres.dat",
                   "gas_init.dat", "aero_emit_dist.dat",
                   "aero_init_comp.dat", "gas_data.dat", "coag_brownian.spec", "run.sh"]
    for i_file, filename in enumerate(input_files):
        dest = shutil.copy("%s" %filename, "%s" % directory)
