# radial-drifter-comparison

get_closest_radial_to_drifter.m creates an array matching the size of the drifter hourly data for each variable you want to match from the radial files from a single site (radial variables are named by their 4-character ID). For each hourly drifter data point, the radial measurement closest to the drifter location is used. Drifter velocity is rotated by the bearing of the closest radial measurement to compare the same velocity component between each. As an option, only radial measurements within some distance of the drifter track are included in the comparison. Paths and some other information in the beginning of the script may need to be changed depending what you need. I'm not familiar with ellipticals so more may need adjusting but I don't know what. If you need to add variables not included in typical radial files, just add the 4-char ID to the 'savetypes' and 'adjustdist' variables.

plot_closest_radial_to_drifter.m uses the output from get_closest_radial_to_drifter.m to create images comparing radial velocity to drifter velocity, plus metrics including spatial and temporal quality. For them to be consistent with mine, you shouldn't really have to change much, but you do need to define printDir as the directory you want to print the images to. Again, I don't know if there are big differences with elliptical data that you'll need to change something in here for.

['USCG_' drifterid '_2016_05_10.mat'] have data for the drifters off NJ (may include only a portion of the drifter track...some may have more data going further into the summer, but the process would be the same with the full dataset if you get that).

Other files are functions called in either get_closest_radial_to_drifter.m or plot_closest_radial_to_drifter.m.
