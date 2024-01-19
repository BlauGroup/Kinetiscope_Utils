# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 10:18:42 2024

@author: jacob
"""

def write_time_and_power(f, initial_time):
    """
    Writes time and corresponding optical power values to a file.

    Parameters:
    - f (file object): The file object to write the data.
    - initial_time (float): The initial time value.

    Returns:
    float: The updated time value after writing the data.

    This function writes four lines to the specified file 'f', each representing a time
    and corresponding optical power value. The time values are rounded to 3 decimal places
    in scientific notation, and the optical power values alternate between 1 and 0.

    The function returns the updated time value after writing the data.
    """
    rounded_time = round(initial_time, 9)  # Round to 9 decimal places
    f.write(f'{rounded_time:.3e} 1\n')

    time = initial_time + 1e-6 - 1e-9
    rounded_time = round(time, 9)  # Round to 9 decimal places
    f.write(f'{rounded_time:.3e} 1\n')

    time = initial_time + 1e-6
    rounded_time = round(time, 9)  # Round to 9 decimal places
    f.write(f'{rounded_time:.3e} 0\n')

    off_interval_start = initial_time + 2.0e-5 - 1e-9
    rounded_off_interval_start = round(off_interval_start, 9)  # Round to 9 decimal places
    f.write(f'{rounded_off_interval_start:.3e} 0\n')

    return initial_time + 2.0e-5

total_exposure_time = 0.136  # Total exposure time in seconds
pulse_frequency = 50000  # Pulse frequency in Hz
pulse_duration = 1e-6  # Pulse duration in seconds

with open("corrected_parameters.txt", "w") as f:
    f.write('time (s) optical power (a.u)\n')  # Column titles

    time = 0.0
    num_pulses = 0

    # Calculate the number of pulses based on total exposure time and pulse frequency
    num_pulses_expected = int(total_exposure_time * pulse_frequency)

    while num_pulses < num_pulses_expected:
        num_pulses += 1
        time = write_time_and_power(f, time)  # Increment to just before the next pulse
    
    # Write the value at 1 second as 0
    f.write('2.000e-1 0\n')

    assert num_pulses == num_pulses_expected





        