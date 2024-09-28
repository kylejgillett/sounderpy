


"""
    SOUNDERPY SOUNDING PARAMETER TIMESERIES

    Purpose of module:

    retrieve sounding data, calculate sounding parameters and plot
    sounding parameters as a timeseries

    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2024
"""

####################################################
# CREATE FIGURE
####################################################
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(20, 15))
####################################################



"""
timeseries_data['times']

"""


# AXIS
####################################################
ax1.plot(timeseries_data['times'], timeseries_data['param'], color='cornflowerblue', ls='-', lw=3, label='SRV')
ax1.set_ylabel('0-500m Streamwise\nVorticity (/sec)', color='cornflowerblue', weight='bold', fontsize=18)
ax1.tick_params(axis='y', labelcolor='cornflowerblue', labelsize=15)
ax1.set_ylim(0, 0.045)
ax1.margins(0)

ax1b = ax1.twinx()
ax1b.plot(times, [val.m for val in srw], 'b-', lw=3, label='SRW')
ax1b.set_ylabel('0-500m Storm\nRelative Wind (RM) (kts)', color='b', weight='bold', fontsize=18)
ax1b.tick_params(axis='y', labelcolor='b', labelsize=15)
ax1b.set_ylim(10, 45)
ax1b.margins(0)
####################################################



# KINEMATICS AXIS
####################################################
ax1.plot(times, srv, color='cornflowerblue', ls='-', lw=3, label='SRV')
ax1.set_ylabel('0-500m Streamwise\nVorticity (/sec)', color='cornflowerblue', weight='bold', fontsize=18)
ax1.tick_params(axis='y', labelcolor='cornflowerblue', labelsize=15)
ax1.set_ylim(0, 0.045)
ax1.margins(0)

ax1b = ax1.twinx()
ax1b.plot(times, [val.m for val in srw], 'b-', lw=3, label='SRW')
ax1b.set_ylabel('0-500m Storm\nRelative Wind (RM) (kts)', color='b', weight='bold', fontsize=18)
ax1b.tick_params(axis='y', labelcolor='b', labelsize=15)
ax1b.set_ylim(10, 45)
ax1b.margins(0)
####################################################
















