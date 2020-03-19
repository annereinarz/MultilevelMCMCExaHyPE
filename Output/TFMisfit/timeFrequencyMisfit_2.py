#
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA

import numpy as np

from obspy.imaging.cm import obspy_sequential, obspy_divergent
from obspy.signal import util
#
#Single Valued Envelope Misfit
def em(st1, st2, dt=0.01, fmin=1., fmax=10., nf=100, w0=6, norm='global',
       st2_isref=True):
    """
    :return: Single Valued Envelope Misfit
    """
    if len(np.shape(st1)) == 1:
        w_1 = np.zeros((1, nf, (np.shape(st1))[0]), dtype=np.complex)
        w_2 = np.zeros((1, nf, (np.shape(st1))[0]), dtype=np.complex)

        w_1[0] = cwt(st1, dt, w0, fmin, fmax, nf)
        w_2[0] = cwt(st2, dt, w0, fmin, fmax, nf)
    else:
        w_1 = np.zeros(((np.shape(st1))[0], nf, (np.shape(st1))[1]), dtype=np.complex)
        w_2 = np.zeros(((np.shape(st2))[0], nf, (np.shape(st2))[1]), dtype=np.complex)

        for i in np.arange((np.shape(st1))[0]):
            w_1[i] = cwt(st1[i], dt, w0, fmin, fmax, nf)
            w_2[i] = cwt(st2[i], dt, w0, fmin, fmax, nf)

    if st2_isref:
        _ar = np.abs(w_2)
    else:
        if np.abs(w_1).max() > np.abs(w_2).max():
            _ar = np.abs(w_1)
        else:
            _ar = np.abs(w_2)

    _em = (np.sum(np.sum((np.abs(w_1) - np.abs(w_2)) ** 2, axis=2),
                  axis=1)) ** .5

    if norm == 'global':
        if len(np.shape(st1)) == 1:
            return _em[0] / (np.sum(_ar ** 2)) ** .5
        else:
            return _em / ((np.sum(np.sum(_ar ** 2, axis=2),
                                  axis=1)) ** .5).max()
    elif norm == 'local':
        if len(np.shape(st1)) == 1:
            return _em[0] / (np.sum(_ar ** 2)) ** .5
        else:
            return _em / (np.sum(np.sum(_ar ** 2, axis=2), axis=1)) ** .5
    else:
        raise ValueError('norm "' + norm + '" not defined!')

############################
#
# Continuous Wavelet Transformation in the Frequency Domain.
def cwt(st, dt, w0, fmin, fmax, nf=100, wl='morlet'):
    from obspy.signal import util
    npts = len(st) * 2
    tmax = (npts - 1) * dt
    t = np.linspace(0., tmax, npts)
    f = np.logspace(np.log10(fmin), np.log10(fmax), nf)

    cwt = np.zeros((npts // 2, nf), dtype=np.complex)

    if wl == 'morlet':

        def psi(t):
            return np.pi ** (-.25) * np.exp(1j * w0 * t) * \
                np.exp(-t ** 2 / 2.)

        def scale(f):
            return w0 / (2 * np.pi * f)
    else:
        raise ValueError('wavelet type "' + wl + '" not defined!')

    nfft = util.next_pow_2(npts) * 2
    sf = np.fft.fft(st, n=nfft)

    # Ignore underflows.
    with np.errstate(under="ignore"):
        for n, _f in enumerate(f):
            a = scale(_f)
            # time shift necessary, because wavelet is defined around t = 0
            psih = psi(-1 * (t - t[-1] / 2.) / a).conjugate() / np.abs(a) ** .5
            psihf = np.fft.fft(psih, n=nfft)
            tminin = int(t[-1] / 2. / (t[1] - t[0]))
            cwt[:, n] = np.fft.ifft(psihf * sf)[tminin:tminin + npts // 2] * \
                (t[1] - t[0])

    return cwt.T
################
#
# Single Valued Phase Misfit
def pm(st1, st2, dt=0.01, fmin=1., fmax=10., nf=100, w0=6, norm='global',
       st2_isref=True):
    """
    Single Valued Phase Misfit
    :return: Single Valued Phase Misfit
    """
    if len(np.shape(st1)) == 1:
        w_1 = np.zeros((1, nf, (np.shape(st1))[0]), dtype=np.complex)
        w_2 = np.zeros((1, nf, (np.shape(st1))[0]), dtype=np.complex)

        w_1[0] = cwt(st1, dt, w0, fmin, fmax, nf)
        w_2[0] = cwt(st2, dt, w0, fmin, fmax, nf)
    else:
        w_1 = np.zeros(((np.shape(st1))[0], nf, (np.shape(st1))[1]), dtype=np.complex)
        w_2 = np.zeros(((np.shape(st2))[0], nf, (np.shape(st2))[1]), dtype=np.complex)

        for i in np.arange((np.shape(st1))[0]):
            w_1[i] = cwt(st1[i], dt, w0, fmin, fmax, nf)
            w_2[i] = cwt(st2[i], dt, w0, fmin, fmax, nf)

    if st2_isref:
        _ar = np.abs(w_2)
    else:
        if np.abs(w_1).max() > np.abs(w_2).max():
            _ar = np.abs(w_1)
        else:
            _ar = np.abs(w_2)

    _pm = np.angle(w_1 / w_2) / np.pi

    _pm = (np.sum(np.sum((_ar * _pm) ** 2, axis=2), axis=1)) ** .5

    if norm == 'global':
        if len(np.shape(st1)) == 1:
            return _pm[0] / (np.sum(_ar ** 2)) ** .5
        else:
            return _pm / ((np.sum(np.sum(_ar ** 2, axis=2),
                                  axis=1)) ** .5).max()
    elif norm == 'local':
        if len(np.shape(st1)) == 1:
            return _pm[0] / (np.sum(_ar ** 2)) ** .5
        else:
            return _pm / (np.sum(np.sum(_ar ** 2, axis=2), axis=1)) ** .5
    else:
        raise ValueError('norm "' + norm + '" not defined!')
####################################################################
#
# plot misfits
def plot_tf_misfits(st1, st2,fout,output_dir, dt=0.01, t0=0.0, fmin=1.0, fmax=10.0, nf=100, w0=6,\
                    norm='global', st2_isref=True, left=0.1, bottom=0.1,\
                    h_1=0.2, h_2=0.125, h_3=0.2, w_1=0.2, w_2=0.6, w_cb=0.01,\
                    d_cb=0.0, show=True, plot_args=['k', 'r', 'b'], ylim=0.0,\
                    clim=0.0, cmap=obspy_divergent):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    npts = (np.shape(st1))[-1]
    tmax = (npts - 1) * dt
    t = np.linspace(0., tmax, npts) + t0
    f = np.logspace(np.log10(fmin), np.log10(fmax), nf)

    # compute time frequency misfits
    _tfem = tfem(st1, st2, dt=dt, fmin=fmin, fmax=fmax, nf=nf, w0=w0,
                 norm=norm, st2_isref=st2_isref)
    _tem = tem(st1, st2, dt=dt, fmin=fmin, fmax=fmax, nf=nf, w0=w0, norm=norm,
               st2_isref=st2_isref)
    _fem = fem(st1, st2, dt=dt, fmin=fmin, fmax=fmax, nf=nf, w0=w0, norm=norm,
               st2_isref=st2_isref)
    _em = em(st1, st2, dt=dt, fmin=fmin, fmax=fmax, nf=nf, w0=w0, norm=norm,
             st2_isref=st2_isref)
    _tfpm = tfpm(st1, st2, dt=dt, fmin=fmin, fmax=fmax, nf=nf, w0=w0,
                 norm=norm, st2_isref=st2_isref)
    _tpm = tpm(st1, st2, dt=dt, fmin=fmin, fmax=fmax, nf=nf, w0=w0, norm=norm,
               st2_isref=st2_isref)
    _fpm = fpm(st1, st2, dt=dt, fmin=fmin, fmax=fmax, nf=nf, w0=w0, norm=norm,
               st2_isref=st2_isref)
    _pm = pm(st1, st2, dt=dt, fmin=fmin, fmax=fmax, nf=nf, w0=w0, norm=norm,
             st2_isref=st2_isref)

    if len(np.shape(st1)) == 1:
        _tfem = _tfem.reshape((1, nf, npts))
        _tem = _tem.reshape((1, npts))
        _fem = _fem.reshape((1, nf))
        _em = _em.reshape((1, 1))
        _tfpm = _tfpm.reshape((1, nf, npts))
        _tpm = _tpm.reshape((1, npts))
        _fpm = _fpm.reshape((1, nf))
        _pm = _pm.reshape((1, 1))
        st1 = st1.reshape((1, npts))
        st2 = st2.reshape((1, npts))
        ntr = 1
    else:
        ntr = (np.shape(st1))[0]

    figs = []
    channel_num='xyz'

    for itr in np.arange(ntr):
        plt.rc('xtick', labelsize=8)
        plt.rc('ytick', labelsize=8)
        plt.figure(figsize=(16,12))
        fig = plt.figure()

        # plot signals
        ax_sig = fig.add_axes([left + w_1+0.05, bottom + h_2 + h_3, w_2-0.05, h_1])
        ax_sig.plot(t, st2[itr], plot_args[0])
        ax_sig.plot(t, st1[itr], plot_args[1])

        # plot TEM
        ax_tem = fig.add_axes([left + w_1+0.05, bottom + h_1 + h_2 + h_3, w_2-0.05, h_2])
        ax_tem.plot(t, _tem[itr], plot_args[2])

        # plot TFEM
        ax_tfem = fig.add_axes([left + w_1+0.05, bottom + h_1 + 2 * h_2 + h_3, w_2-0.05,
                                h_3])

        x, y = np.meshgrid(
            t, np.logspace(np.log10(fmin), np.log10(fmax),
                           _tfem[itr].shape[0]))
        img_tfem = ax_tfem.pcolormesh(x, y, _tfem[itr], cmap=cmap)
        img_tfem.set_rasterized(True)
        ax_tfem.set_yscale("log")
        ax_tfem.set_ylim(fmin, fmax)
        #ax_tfem.set_ylabel('frequency')

        # plot FEM
        ax_fem = fig.add_axes([left-0.01, bottom + h_1 + 2 * h_2 + h_3, w_1, h_3])
        ax_fem.semilogy(_fem[itr], f, plot_args[2])
        ax_fem.set_ylim(fmin, fmax)

        # plot TPM
        ax_tpm = fig.add_axes([left + w_1+0.05, bottom, w_2-0.05, h_2])
        ax_tpm.plot(t, _tpm[itr], plot_args[2])

        # plot TFPM
        ax_tfpm = fig.add_axes([left + w_1+0.05, bottom + h_2, w_2-0.05, h_3])

        x, y = np.meshgrid(t, f)
        img_tfpm = ax_tfpm.pcolormesh(x, y, _tfpm[itr], cmap=cmap)
        img_tfpm.set_rasterized(True)
        ax_tfpm.set_yscale("log")
        ax_tfpm.set_ylim(f[0], f[-1])
        #ax_tfpm.set_ylabel('frequency')

        # add colorbars
        ax_cb_tfpm = fig.add_axes([left + w_1 + w_2 + d_cb + w_cb, bottom,
                                   w_cb, h_2 + h_3])
        fig.colorbar(img_tfpm, cax=ax_cb_tfpm)

        # plot FPM
        ax_fpm = fig.add_axes([left-0.01, bottom + h_2, w_1, h_3])
        ax_fpm.semilogy(_fpm[itr], f, plot_args[2])
        ax_fpm.set_ylim(fmin, fmax)

        # set limits
        ylim_sig = np.max([np.abs(st1).max(), np.abs(st2).max()]) * 1.1
        ax_sig.set_ylim(-ylim_sig, ylim_sig)

        if ylim == 0.:
            ylim = np.max([np.abs(_tem).max(), np.abs(_tpm).max(),
                           np.abs(_fem).max(), np.abs(_fpm).max()]) * 1.1

        ax_tem.set_ylim(-ylim, ylim)
        ax_fem.set_xlim(-ylim, ylim)
        ax_tpm.set_ylim(-ylim, ylim)
        ax_fpm.set_xlim(-ylim, ylim)

        ax_sig.set_xlim(t[0], t[-1])
        ax_tem.set_xlim(t[0], t[-1])
        ax_tpm.set_xlim(t[0], t[-1])

        if clim == 0.:
            clim = np.max([np.abs(_tfem).max(), np.abs(_tfpm).max()])

        img_tfpm.set_clim(-clim, clim)
        img_tfem.set_clim(-clim, clim)

        # add text box for EM + PM
        textstr = 'EM = %.5f\nPM = %.5f' % (_em[itr], _pm[itr])
        print(_em[itr],',', _pm[itr])
        props = dict(boxstyle='round', facecolor='white')
        ax_sig.text(-0.5, 0.5, textstr, transform=ax_sig.transAxes,
                    verticalalignment='center', horizontalalignment='left',
                    bbox=props)

        ax_tpm.set_xlabel('time')
        ax_fem.set_ylabel('frequency')
        ax_fpm.set_ylabel('frequency')

        # add text boxes
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        ax_tfem.text(0.95, 0.85, 'TFEM', transform=ax_tfem.transAxes,
                     verticalalignment='top', horizontalalignment='right',
                     bbox=props, fontsize=10)
        ax_tfpm.text(0.95, 0.85, 'TFPM', transform=ax_tfpm.transAxes,
                     verticalalignment='top', horizontalalignment='right',
                     bbox=props, fontsize=10)
        ax_tem.text(0.95, 0.75, 'TEM', transform=ax_tem.transAxes,
                    verticalalignment='top', horizontalalignment='right',
                    bbox=props)
        ax_tpm.text(0.95, 0.75, 'TPM', transform=ax_tpm.transAxes,
                    verticalalignment='top', horizontalalignment='right',
                    bbox=props)
        ax_fem.text(0.9, 0.85, 'FEM', transform=ax_fem.transAxes,
                    verticalalignment='top', horizontalalignment='right',
                    bbox=props)
        ax_fpm.text(0.9, 0.85, 'FPM', transform=ax_fpm.transAxes,
                    verticalalignment='top', horizontalalignment='right',
                    bbox=props)

        # remove axis labels
        ax_tfpm.xaxis.set_major_formatter(NullFormatter())
        ax_tfem.xaxis.set_major_formatter(NullFormatter())
        ax_tem.xaxis.set_major_formatter(NullFormatter())
        ax_sig.xaxis.set_major_formatter(NullFormatter())
        #ax_tfpm.yaxis.set_major_formatter(NullFormatter())
        #ax_tfem.yaxis.set_major_formatter(NullFormatter())

        figs.append(fig)
        fig.savefig(output_dir+fout+'%s.pdf'%channel_num[itr])

    if show:
        plt.show()
    else:
        if ntr == 1:
            return figs[0]
        else:
            return figs
###################################################
#
#Time Frequency Envelope Misfit
def tfem(st1, st2, dt=0.01, fmin=1., fmax=10., nf=100, w0=6, norm='global',
         st2_isref=True):
    """
    Time Frequency Envelope Misfit

    """
    if len((np.shape(st1))) == 1:
        w_1 = np.zeros((1, nf, (np.shape(st1))[0]), dtype=np.complex)
        w_2 = np.zeros((1, nf, (np.shape(st1))[0]), dtype=np.complex)

        w_1[0] = cwt(st1, dt, w0, fmin, fmax, nf)
        w_2[0] = cwt(st2, dt, w0, fmin, fmax, nf)
    else:
        w_1 = np.zeros(((np.shape(st1))[0], nf, (np.shape(st1))[1]), dtype=np.complex)
        w_2 = np.zeros(((np.shape(st2))[0], nf, (np.shape(st2))[1]), dtype=np.complex)

        for i in np.arange((np.shape(st1))[0]):
            w_1[i] = cwt(st1[i], dt, w0, fmin, fmax, nf)
            w_2[i] = cwt(st2[i], dt, w0, fmin, fmax, nf)

    if st2_isref:
        ar = np.abs(w_2)
    else:
        if np.abs(w_1).max() > np.abs(w_2).max():
            ar = np.abs(w_1)
        else:
            ar = np.abs(w_2)

    _tfem = (np.abs(w_1) - np.abs(w_2))

    if norm == 'global':
        if len(np.shape(st1)) == 1:
            return _tfem[0] / np.max(ar)
        else:
            return _tfem / np.max(ar)
    elif norm == 'local':
        if len(np.shape(st1)) == 1:
            return _tfem[0] / ar[0]
        else:
            return _tfem / ar
    else:
        raise ValueError('norm "' + norm + '" not defined!')
############################################################
#
# Time-dependent Envelope Misfit
def tem(st1, st2, dt=0.01, fmin=1., fmax=10., nf=100, w0=6, norm='global',
        st2_isref=True):
    """
    Time-dependent Envelope Misfit
    """
    if len(np.shape(st1)) == 1:
        w_1 = np.zeros((1, nf, (np.shape(st1))[0]), dtype=np.complex)
        w_2 = np.zeros((1, nf, (np.shape(st1))[0]), dtype=np.complex)

        w_1[0] = cwt(st1, dt, w0, fmin, fmax, nf)
        w_2[0] = cwt(st2, dt, w0, fmin, fmax, nf)
    else:
        w_1 = np.zeros(((np.shape(st1))[0], nf, (np.shape(st1))[1]), dtype=np.complex)
        w_2 = np.zeros(((np.shape(st2))[0], nf, (np.shape(st2))[1]), dtype=np.complex)

        for i in np.arange((np.shape(st1))[0]):
            w_1[i] = cwt(st1[i], dt, w0, fmin, fmax, nf)
            w_2[i] = cwt(st2[i], dt, w0, fmin, fmax, nf)

    if st2_isref:
        _ar = np.abs(w_2)
    else:
        if np.abs(w_1).max() > np.abs(w_2).max():
            _ar = np.abs(w_1)
        else:
            _ar = np.abs(w_2)

    _tem = np.sum((np.abs(w_1) - np.abs(w_2)), axis=1)

    if norm == 'global':
        if len(np.shape(st1)) == 1:
            return _tem[0] / np.max(np.sum(_ar, axis=1))
        else:
            return _tem / np.max(np.sum(_ar, axis=1))
    elif norm == 'local':
        if len(np.shape(st1)) == 1:
            return _tem[0] / np.sum(_ar, axis=1)[0]
        else:
            return _tem / np.sum(_ar, axis=1)
    else:
        raise ValueError('norm "' + norm + '" not defined!')
########################################################
#
# Frequency-dependent Envelope Misfit
def fem(st1, st2, dt=0.01, fmin=1., fmax=10., nf=100, w0=6, norm='global',
        st2_isref=True):
    """
    Frequency-dependent Envelope Misfit
    """
    if len(np.shape(st1)) == 1:
        w_1 = np.zeros((1, nf, (np.shape(st1))[0]), dtype=np.complex)
        w_2 = np.zeros((1, nf, (np.shape(st1))[0]), dtype=np.complex)

        w_1[0] = cwt(st1, dt, w0, fmin, fmax, nf)
        w_2[0] = cwt(st2, dt, w0, fmin, fmax, nf)
    else:
        w_1 = np.zeros(((np.shape(st1))[0], nf, (np.shape(st1))[1]), dtype=np.complex)
        w_2 = np.zeros(((np.shape(st2))[0], nf, (np.shape(st2))[1]), dtype=np.complex)

        for i in np.arange((np.shape(st1))[0]):
            w_1[i] = cwt(st1[i], dt, w0, fmin, fmax, nf)
            w_2[i] = cwt(st2[i], dt, w0, fmin, fmax, nf)

    if st2_isref:
        _ar = np.abs(w_2)
    else:
        if np.abs(w_1).max() > np.abs(w_2).max():
            _ar = np.abs(w_1)
        else:
            _ar = np.abs(w_2)

    _tem = np.abs(w_1) - np.abs(w_2)
    _tem = np.sum(_tem, axis=2)

    if norm == 'global':
        if len(np.shape(st1)) == 1:
            return _tem[0] / np.max(np.sum(_ar, axis=2))
        else:
            return _tem / np.max(np.sum(_ar, axis=2))
    elif norm == 'local':
        if len(np.shape(st1)) == 1:
            return _tem[0] / np.sum(_ar, axis=2)[0]
        else:
            return _tem / np.sum(_ar, axis=2)
    else:
        raise ValueError('norm "' + norm + '" not defined!')
###################################################################
#
# Time Frequency Phase Misfit
def tfpm(st1, st2, dt=0.01, fmin=1., fmax=10., nf=100, w0=6, norm='global',
         st2_isref=True):
    """
    Time Frequency Phase Misfit
    """
    if len(np.shape(st1)) == 1:
        w_1 = np.zeros((1, nf, (np.shape(st1))[0]), dtype=np.complex)
        w_2 = np.zeros((1, nf, (np.shape(st1))[0]), dtype=np.complex)

        w_1[0] = cwt(st1, dt, w0, fmin, fmax, nf)
        w_2[0] = cwt(st2, dt, w0, fmin, fmax, nf)
    else:
        w_1 = np.zeros(((np.shape(st1))[0], nf, (np.shape(st1))[1]), dtype=np.complex)
        w_2 = np.zeros(((np.shape(st2))[0], nf, (np.shape(st2))[1]), dtype=np.complex)

        for i in np.arange((np.shape(st1))[0]):
            w_1[i] = cwt(st1[i], dt, w0, fmin, fmax, nf)
            w_2[i] = cwt(st2[i], dt, w0, fmin, fmax, nf)

    if st2_isref:
        _ar = np.abs(w_2)
    else:
        if np.abs(w_1).max() > np.abs(w_2).max():
            _ar = np.abs(w_1)
        else:
            _ar = np.abs(w_2)
    
    _tfpm = np.angle(w_1 / w_2) / np.pi

    if norm == 'global':
        if len(np.shape(st1)) == 1:
            return _ar[0] * _tfpm[0] / np.max(_ar)
        else:
            return _ar * _tfpm / np.max(_ar)
    elif norm == 'local':
        if len(np.shape(st1)) == 1:
            return _tfpm[0]
        else:
            return _tfpm
    else:
        raise ValueError('norm "' + norm + '" not defined!')
##############################################################
#
# Time-dependent Phase Misfit
def tpm(st1, st2, dt=0.01, fmin=1., fmax=10., nf=100, w0=6, norm='global',
        st2_isref=True):
    """
    Time-dependent Phase Misfit
    """
    if len(np.shape(st1)) == 1:
        w_1 = np.zeros((1, nf, (np.shape(st1))[0]), dtype=np.complex)
        w_2 = np.zeros((1, nf, (np.shape(st1))[0]), dtype=np.complex)

        w_1[0] = cwt(st1, dt, w0, fmin, fmax, nf)
        w_2[0] = cwt(st2, dt, w0, fmin, fmax, nf)
    else:
        w_1 = np.zeros(((np.shape(st1))[0], nf, (np.shape(st1))[1]), dtype=np.complex)
        w_2 = np.zeros(((np.shape(st2))[0], nf, (np.shape(st2))[1]), dtype=np.complex)

        for i in np.arange((np.shape(st1))[0]):
            w_1[i] = cwt(st1[i], dt, w0, fmin, fmax, nf)
            w_2[i] = cwt(st2[i], dt, w0, fmin, fmax, nf)

    if st2_isref:
        _ar = np.abs(w_2)
    else:
        if np.abs(w_1).max() > np.abs(w_2).max():
            _ar = np.abs(w_2)
        else:
            _ar = np.abs(w_1)

    _tpm = np.angle(w_1 / w_2) / np.pi
    _tpm = np.sum(_ar * _tpm, axis=1)

    if norm == 'global':
        if len(np.shape(st1)) == 1:
            return _tpm[0] / np.max(np.sum(_ar, axis=1))
        else:
            return _tpm / np.max(np.sum(_ar, axis=1))
    elif norm == 'local':
        if len(np.shape(st1)) == 1:
            return _tpm[0] / np.sum(_ar, axis=1)[0]
        else:
            return _tpm / np.sum(_ar, axis=1)
    else:
        raise ValueError('norm "' + norm + '" not defined!')
##################################################################
#
# Frequency-dependent Phase Misfit
def fpm(st1, st2, dt=0.01, fmin=1., fmax=10., nf=100, w0=6, norm='global',
        st2_isref=True):
    """
    Frequency-dependent Phase Misfit
    """
    if len(np.shape(st1)) == 1:
        w_1 = np.zeros((1, nf, (np.shape(st1))[0]), dtype=np.complex)
        w_2 = np.zeros((1, nf, (np.shape(st1))[0]), dtype=np.complex)

        w_1[0] = cwt(st1, dt, w0, fmin, fmax, nf)
        w_2[0] = cwt(st2, dt, w0, fmin, fmax, nf)
    else:
        w_1 = np.zeros(((np.shape(st1))[0], nf, (np.shape(st1))[1]), dtype=np.complex)
        w_2 = np.zeros(((np.shape(st2))[0], nf, (np.shape(st2))[1]), dtype=np.complex)

        for i in np.arange((np.shape(st1))[0]):
            w_1[i] = cwt(st1[i], dt, w0, fmin, fmax, nf)
            w_2[i] = cwt(st2[i], dt, w0, fmin, fmax, nf)

    if st2_isref:
        _ar = np.abs(w_2)
    else:
        if np.abs(w_1).max() > np.abs(w_2).max():
            _ar = np.abs(w_1)
        else:
            _ar = np.abs(w_2)

    _tpm = np.angle(w_1 / w_2) / np.pi
    _tpm = np.sum(_ar * _tpm, axis=2)

    if norm == 'global':
        if len(np.shape(st1)) == 1:
            return _tpm[0] / np.max(np.sum(_ar, axis=2))
        else:
            return _tpm / np.max(np.sum(_ar, axis=2))
    elif norm == 'local':
        if len(np.shape(st1)) == 1:
            return _tpm[0] / np.sum(_ar, axis=2)[0]
        else:
            return _tpm / np.sum(_ar, axis=2)
    else:
        raise ValueError('norm "' + norm + '" not defined!')
################################################################
