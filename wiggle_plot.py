"""
Defines a function for plotting seismograms

Nathan Stoikopoulos
University of Toronto
2018-09-22

This script is a port to python (with some vectorization and other minor
modifications) of a MATLAB wiggleplot function,
CGB, 18 Aug 1999, modified 13 Oct 2000
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm


def prep_plot(traces, t_start, dt, geophones, trace_normed, scale, mode='wiggle'):
    n_samples = traces.shape[0]  # Samples per trace
    n_traces = traces.shape[1]  # Number of traces
    geophone_spacing = np.abs(geophones[1] - geophones[0])
    t_end = t_start + (n_samples - 1) * dt  # Max time on seismogram
    x, y = np.meshgrid(np.arange(0, n_samples, 1), geophones)
    x, y = x.T, y.T
    times = np.arange(t_start, t_end + dt, dt)
    if trace_normed:
        # Normalize each trace individually
        plot_traces = y + (np.max(np.abs(traces), 0).T)**-1 * traces \
                      * (geophone_spacing if mode == 'wiggle' else 1)
    else:
        # Normalize all traces as a whole
        normval = np.max(np.abs(traces))
        plot_traces = y + scale/normval * traces

    return times, plot_traces

def imageplot(traces, t_start, dt, geophones, trace_normed, scale, title='', reverse=False, **kwargs):

    times, plot_traces = prep_plot(traces, t_start, dt, geophones, trace_normed, scale)
    ax = plt.gca()
    ax.imshow(traces, aspect= )

def wiggleplot(traces, t_start, dt, geophones, trace_normed, scale, title='', reverse=False, **kwargs):
    """
    Produces a seismogram plot, with location (and amplitude) along the y-axis,
    and time on the x-axis
    :param ndarray traces: Trace data to be plotted
    :param float t_start: Time at which seismogram should begin
    :param float dt: Time-step between samples
    :param ndarray geophones: Row-vector of geophone locations
    :param bool trace_normed: Normalizes trace-by-trace if true, otherwise
                               normalizes relative to whole data set
    :param float scale: Scaling factor for traces
    :param dict kwargs: Keyword arguments to be passed to plot() function
    :return: None
    """
    #plt.figure(400)
    times, plot_traces = prep_plot(traces, t_start, dt,
                                   geophones, trace_normed, scale)

    ax = plt.gca()
    plt.plot(plot_traces, times, 'k-', **kwargs)
    geophone_ones = geophones * np.ones(plot_traces.shape)
    for gp, trace in enumerate(plot_traces.T):
        ax.fill_betweenx(times, geophone_ones[:, gp], trace, where=trace > geophone_ones[:, gp], color='k')
        #ax.fill_betweenx(times, geophone_ones[:, gp], trace, where=trace < geophone_ones[:, gp], color='k')
    #plt.plot(times, plot_traces[:, 18], 'r-')
    ax.invert_yaxis()
    plt.title(title)
    plt.xlabel("Position (m)")
    plt.ylabel("Time (s)")
    return ax, plot_traces
