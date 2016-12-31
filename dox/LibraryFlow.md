# The Flow of ForceManII                      {#flow}
========================

The call sequence of the ForceManII library will be very similar in almost all
cases with some exceptions when the user desires more customability.

In general the call sequence will be:

1. `parse_file` to generate a ForceField instance
2. `get_coords` to generate a CoordArray instance for the system
3. `get_params` to parameterize each coordinate of the system
4. `deriv` to compute the requested derivative

For convenience we define a function `run_forcemanii` that does this sequence
for you returning the requested derivative.


