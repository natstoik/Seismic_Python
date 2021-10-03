"""
Script for picking first arrivals in seismic data stored in SEG2 format. Picks
are saved in numpy.ndarray format as $(FILE_NAME)_picks.txt

Nathan Stoikopoulos
University of Toronto
2018-09-22

This script is based on MATLAB script seis_pick.m by Carl-Georg Bank, included
in SIGKIT software package
"""


import numpy as np
import pathlib as pth
import matplotlib
import matplotlib.pyplot as plt
from wiggle_plot import wiggleplot
import os
import scipy.signal as sig
import seg2load
matplotlib.use('TkAgg', force=True)

def butter_bandpass(lowcut, highcut, fs, order=5):
    """
    Taken from
    https://scipy-cookbook.readthedocs.io/items/ButterworthBandpass.html

    Creates a butterworth bandpass filter of order 'order', over frequency
    band [lowcut, highcut].
    :param lowcut: Lowcut frequency in Hz
    :param highcut: Highcut frequency in Hz
    :param fs: Sampling frequency in Hz
    :param order: width of influence in points
    :return filt: Filter to be passed to butter_bandpass_filter for evaluation
    """
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = sig.butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    """
    Taken from
    https://scipy-cookbook.readthedocs.io/items/ButterworthBandpass.html

    Creates a bandpass filter, and applies it to data
    :param data:
    :param lowcut: 
    :param highcut:
    :param fs:
    :param order:
    :return:
    """
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = sig.lfilter(b, a, data.T)
    return y.T

def binary_airwave(shot, data, geophones, dt, direct_zone = 4, speed = 340):
    """
    Binary filter which sets all signal after the airwave arrival (assuming
    velocity of 340 m/s) from a shot at location 'shot'
    :param shot: Coordinate of shot in m
    :param data: Geophone responses
    :param geophones: Locations of geophones in m
    :param dt: Sampling time of data
    :param direct_zone: Radius in m around shot in which traces are not filtered
    :return: filtered_data
    """
    filt = np.ones_like(data) * np.array([np.mgrid[0:data.shape[0]*dt:dt]]).T
    phone_dist = np.abs(geophones - shot)
    airwave = phone_dist / speed
    filt[np.logical_and(filt >= airwave, phone_dist >= direct_zone)] = 0
    print(filt[:, 5])
    filt[filt > 0] = 1
    filt[0, :] = 1
    return data * filt


def read_seg2(data_file):

    seis_data, seis_header = seg2load.seg2_load(str(data_file))
    dt = seis_header['tr']['SAMPLE_INTERVAL'][0]  # sampling interval
    header = seis_header
    return seis_data, dt, header


def ps_deconvolve(p_file, s_file, p_pick_file):
    p, dt, header = read_seg2(p_file)
    ps, dt, header = read_seg2(s_file)

    def fft_deconvolve(f, conv):

        conv_ft = np.fft.fft(conv)
        f_ft = np.fft.fft(f, len(conv_ft), 0)
        deconv_ft = f_ft /conv_ft
        deconv = np.fft.ifft(deconv_ft)

        return deconv

    s = fft_deconvolve(p, ps)
    return s, dt


def load_picks(pick_file):
    return np.loadtxt(pick_file, skiprows=1)


def pick_arrivals():
    picks = np.asarray(plt.ginput(n=-1, show_clicks=True, timeout=0))
    picks_x = picks[:, 0]
    picks_y = picks[:, 1]
    sorted_picks = picks_y[np.argsort(picks_x)]
    return sorted_picks


def seis_pick(data_file, t_start, t_end, geophones, trace_normed=True,
              p_pickfile=None, bandpass=None, remove_airwave=False, shot=False):

    # Where first-arrival pick coordinates are saved in ndarray format
    pick_file = '{}_picks'.format(data_file.name)
    pick_path = './picks/'
    if os.path.exists(pick_path):
        pass
    else:
        os.mkdir('./picks/')

    seis_view(data_file, t_start, t_end, geophones,
                   trace_normed=trace_normed, bandpass=bandpass, show=False)
    if p_pickfile:
        p_picks = np.loadtxt(pick_path + p_pickfile, skiprows=1)
        plt.scatter(p_picks[:, 1], p_picks[:, 0], c='k', marker='o')

    # picks = load_picks(pick_path+'373.dat_picks')
    sorted_picks = pick_arrivals()
    np.savetxt(pick_path + pick_file, np.array([geophones, sorted_picks]).T,
               header='Position (m)     Twt (s)', fmt='%2f %.18e')
    # plt.plot(picks_x, picks_y, c='b', marker='+')
    plt.show()


def flatten_arrivals(data_file, t_start, t_end, geophones, bandpass=None, title=None):
    ax, traces, dt = seis_view(data_file, t_start, t_end, geophones, bandpass=bandpass, show=False, trace_normed=True)
    picks = pick_arrivals()
    arrival_indices = np.asarray(picks // dt, dtype=int)
    if len(arrival_indices) != traces.shape[1]:
        raise ValueError('Wrong number of picks! ' + str(len(arrival_indices)))
    print(traces.shape)
    for trace, arrival in enumerate(arrival_indices):
        traces[:-arrival, trace] = traces[arrival:, trace]
        traces[-arrival:, trace] = geophones[trace]
    plt.figure()
    wiggleplot(traces, t_start, dt, geophones, trace_normed=True, scale=2, title=title)
    plt.ylim([(len(traces) - np.max(arrival_indices))*dt, 0])
    plt.show()


def seis_view(data_file, t_start, t_end, geophones, trace_normed=True,
              bandpass=None, show=True, title=None):

    seis_data, dt, header = read_seg2(data_file)

    # Trim and filter traces
    trace_bounds = np.array([t_start, (t_end+dt)]) // dt
    seis_data = seis_data[int(trace_bounds[0]):int(trace_bounds[1]), :]
    if not geophones:
        geophones = header['tr']['RECEIVER_LOCATION'].flatten()
    if bandpass:
        seis_data = butter_bandpass_filter(seis_data, bandpass[0], bandpass[1], 1/dt, order=5)


    #Plot seismogram
    plt.tight_layout()
    if not title:
        title = data_file.name
    ax, plot_traces = wiggleplot(seis_data, t_start, dt, geophones, trace_normed, 1,
               linewidth=0.5, title=title, reverse=True)

   # plt.xlim([t_start, t_end])
   # ax = plt.gca()
    if show:
        plt.show()
    return ax, seis_data, dt


def cmp_gather(file_info):
    gathers = {}
    for file in file_info:
        traces, dt = read_seg2(file[0])
        receivers = np.linspace(file[2], file[3], traces.shape[1])
        offsets = receivers - file[1]
        midpoints = receivers - offsets / 2
        gather_keys = gathers.keys()
        for col, mp in enumerate(midpoints):
            if mp in gather_keys:
                traces[-1, col] = offsets[col]
                print(traces.shape)
                try:
                    gathers[mp] = np.array([gathers[mp], traces[:, col]]).T
                except ValueError:
                    pass
            else:
                gathers[mp] = traces[:, col]
    for mp in gathers.keys():
        traces = gathers[mp]
        print(traces.shape)
        if traces.ndim ==2 and traces.shape[1]>=3:
            offsets = traces[-1, :]
            wiggleplot(traces, 0, dt, offsets, True, 1,
                       title='Gathers at '+str(mp)+'m')
            plt.show()


if __name__ == '__main__':
    # Path to folder containing .dat files
    data_path = pth.PurePosixPath('./Nash_Creek_Seismic/UNSORTED_TRIP_2/')
    data_file = data_path / '415.dat'  # Append file name to path
    t_min = 0.0
    t_max = 0.5  # Upper limit on time axis
    geop = np.linspace(0, 115, 24, dtype=np.float64)

    # file_info = [['503.dat', 465, 490, 605], ['504.dat', 482.5, 490, 605],
    #             ['505.dat', 542.5, 490, 605], ['507.dat', 587.5, 490, 605]]
    # for file in file_info:
    #    file[0] = str(data_path/file[0])
    # cmp_gather(file_info)
    picture_path = data_path / 'Pictures'
    for f in os.scandir(data_path):
        if f.is_file():
            if int(f.name.split('.')[0]) in range(550, 565):
                data_file = data_path / f.name
                seis_view(data_file, t_min, t_max, None, show=False)
                plt.savefig(picture_path / '{}.pdf'.format(f.name.split('.')[0]))
                plt.show()
                plt.close()
    # flatten_arrivals(data_file, t_min, t_max, geop, bandpass=None, title ='Area 2 Line 3 - Side Trail')
