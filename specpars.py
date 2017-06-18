import numpy as np
import os
import logging
import matplotlib.pyplot as plt
from . import moog, errors
from .tools import linfit, wlinear_fit
from .star import Star
import datetime
from scipy import ma
from collections import OrderedDict
from bokeh.plotting import *
from bokeh.models import HoverTool
import copy

logger = logging.getLogger(__name__)


class SolvePars:
    def __init__(self, grid='odfnew'):
        self.step_teff = 32
        self.step_logg = 0.32
        self.step_vt = 0.32
        self.niter = 50
        self.grid = grid
        self.solar_afe = 7.50
        self.errors = False
        self.check_converged = True
        self.ignore = []

        # LM boundaries for parameter exploration
        self.teff_min = 4100.0
        self.teff_max = 7000.0

        self.vmic_min = 0.0
        self.vmic_max = 4.0

        self.logg_min = 1.5
        self.logg_max = 4.95

        self.gfeh_min = -2.5
        self.gfeh_max = 0.50

        # LM: new variables for alternative converging algortihm
        self.step_teff_alt = 50.0
        self.step_logg_alt = 0.10
        self.step_vmic_alt = 0.10
        self.stef_gfeh_alt = 0.10

        # Number of steps in the interpolation scheme of each variable
        self.inter_halfsteps = 3

        # critical limits for iteration stopping criteria
        self.teff_crit = [10.0, 50.0]
        self.vmic_crit = [0.01, 0.05]
        self.logg_crit = [0.01, 0.10]
        self.gfeh_crit = [0.01, 0.10]

        self.count_iters_limit = 3
        self.count_recurs_limit = 25
        self.iter_poly_order = 2
        self.QP = [0.25, 0.25, 1., 2.]
        self.QP_sigma_rejection = [5.0, 4.0, 3.0, 2.0]


class PlotPars:
    def __init__(self):
        self.afe = None
        self.wavelength_range = None
        self.make_figure = True
        self.figure_format = 'png'
        self.title = None
        self.title_inside = None


def iron_stats(Star, Ref=object, plot=None, PlotPars=object, silent=True):
    if hasattr(Ref, 'name'):
        if Star.name == Ref.name:
            x = {'afe': 0, 'err_afe': 0,
                 'afe1': 0, 'err_afe1': 0, 'nfe1': 0,
                 'afe2': 0, 'err_afe2': 0, 'nfe2': 0,
                 'slope_ep': 0,
                 'err_slope_ep': 0,
                 'slope_rew': 0,
                 'err_slope_rew': 0,
                 'rank': 0,
                 'reference': Ref.name}
            Star.iron_stats = x
            return None
    logger.info('Begin iron_stats for ' + Star.name)
    logger.info('Calculating abundances for ' + Star.name)
    fe1_done = moog.abfind(Star, 26.0, 'fe1')
    fe2_done = moog.abfind(Star, 26.1, 'fe2')
    if not fe1_done and not fe2_done:
        logger.warning('No fe1/fe2 attribute(s) added to ' + Star.name)
        return None

    # LM create a conditionnl array if not present already
    if 'flag' not in Star.fe1:
        Star.fe1['flag'] = np.ones(np.size(Star.fe1['ep']), dtype=bool)
    if 'flag' not in Star.fe2:
        Star.fe2['flag'] = np.ones(np.size(Star.fe2['ep']), dtype=bool)

    if hasattr(Ref, 'name'):
        logger.info('Differential analysis. Reference star is ' + Ref.name)
        if not (hasattr(Ref, 'fe1')):
            logger.info('Reference star does not have abundances as ' + \
                        'attributes')
            logger.info('Calculating abundances for reference star')
            moog.abfind(Ref, 26.0, 'fe1')
            moog.abfind(Ref, 26.1, 'fe2')
        ww1, ww2 = Star.fe1['ww'], Star.fe2['ww']
        ww1r, ww2r = Ref.fe1['ww'], Ref.fe2['ww']
        w1, w2 = np.intersect1d(ww1, ww1r), np.intersect1d(ww2, ww2r)
        k1 = [i for i, w in zip(range(len(ww1)), ww1) if w in w1]
        k1r = [i for i, w in zip(range(len(ww1r)), ww1r) if w in w1]
        k2 = [i for i, w in zip(range(len(ww2)), ww2) if w in w2]
        k2r = [i for i, w in zip(range(len(ww2r)), ww2r) if w in w2]
        afe1 = Star.fe1['ab'][k1] - Ref.fe1['ab'][k1r]
        afe2 = Star.fe2['ab'][k2] - Ref.fe2['ab'][k2r]
        afe1_e = Star.fe1['ab_e'][k1]  # LM Reference star is supposed error-less
        afe2_e = Star.fe2['ab_e'][k2]  # LM Reference star is supposed error-less
        rew1 = np.log10(1e-3 * Star.fe1['ew'][k1] / Star.fe1['ww'][k1])
        rew2 = np.log10(1e-3 * Star.fe2['ew'][k2] / Star.fe2['ww'][k2])
        ep1, ep2 = Star.fe1['ep'][k1], Star.fe2['ep'][k2]
        w1 = Star.fe1['ww'][k1]
        w2 = Star.fe2['ww'][k2]
        #
        Star.fe1['ww'], Star.fe2['ww'] = w1, w2
        Star.fe1['ep'], Star.fe2['ep'] = ep1, ep2
        Star.fe1['ew'], Star.fe2['ew'] = Star.fe1['ew'][k1], Star.fe2['ew'][k2]
        Star.fe1['rew'], Star.fe2['rew'] = rew1, rew2
        Star.fe1['ab'], Star.fe2['ab'] = Star.fe1['ab'][k1], Star.fe2['ab'][k2]
        Star.fe1['ab_e'], Star.fe2['ab_e'] = Star.fe1['ab_e'][k1], Star.fe2['ab_e'][k2]
        Star.fe1['difab'], Star.fe2['difab'] = afe1, afe2
        # Added by LM
        Star.fe1['flag'], Star.fe2['flag'] = Star.fe1['flag'][k1], Star.fe2['flag'][k2]

        if plot:
            # ylabel = '$\Delta$[Fe/H]'
            ylabel = '[Fe/H]'
    else:
        logger.info('Working with absolute abundances')
        w1, w2 = Star.fe1['ww'], Star.fe2['ww']
        afe1 = Star.fe1['ab']
        afe2 = Star.fe2['ab']
        afe1_e = Star.fe1['ab_e']
        afe2_e = Star.fe2['ab_e']
        rew1 = np.log10(1e-3 * Star.fe1['ew'] / w1)
        rew2 = np.log10(1e-3 * Star.fe2['ew'] / w2)
        ep1, ep2 = Star.fe1['ep'], Star.fe2['ep']
        if plot:
            ylabel = 'A(Fe)'

    mfe1, efe1 = np.mean(afe1), np.std(afe1, ddof=1)
    mfe2, efe2 = np.mean(afe2), np.std(afe2, ddof=1)
    mafe = np.mean(list(afe1) + list(afe2))
    eafe = np.std(list(afe1) + list(afe2))
    nfe1, nfe2 = len(afe1), len(afe2)

    # LM  Linear fit is substituted with weighted fit
    # afe1_v = 1./(afe1_e**2)
    afe1_v = (0.05 / afe1_e) ** 2
    afe1_v[np.asarray(afe1_e) < 0.05] = 1.00  # to avoid giving too much weight to unrealistically small error bars
    if np.sum(Star.fe1['flag']) < np.size(Star.fe1['flag']):
        afe1_v[not Star.fe1['flag']] = 0.00

    zero_ep, slope_ep, _ = wlinear_fit(ep1, afe1, afe1_v)
    zero_rew, slope_rew, _ = wlinear_fit(rew1, afe1, afe1_v)
    _, _, err_slope_ep = linfit(ep1, afe1)
    _, _, err_slope_rew = linfit(rew1, afe1)
    #zero_ep, slope_ep, err_slope_ep = linfit(ep1, afe1)
    #zero_rew, slope_rew, err_slope_rew = linfit(rew1, afe1)


    x_epfit = np.array([min(ep1), max(ep1)])
    y_epfit = zero_ep + slope_ep * x_epfit
    x_rewfit = np.array([min(rew1), max(rew1)])
    y_rewfit = zero_rew + slope_rew * x_rewfit

    if plot:
        logger.info('Making figure')
        plt.figure(figsize=(7, 9))
        title = Star.name + ' : ' + str(Star.teff) + ', ' + str(Star.logg) + ', ' \
                + str(round(Star.feh, 3)) + ', ' + str(Star.vt)
        if hasattr(Ref, 'name'):
            title += ' [' + Ref.name + ']'
        if hasattr(PlotPars, 'title'):
            if PlotPars.title != None:
                title = PlotPars.title
        plt.suptitle(title, fontsize=16)
        plt.subplots_adjust(hspace=0.35, top=0.93, left=0.2)
        plt.rc("axes", labelsize=15, titlesize=12)
        plt.rc("xtick", labelsize=14)
        plt.rc("ytick", labelsize=14)
        plt.rc("xtick.major", size=6, width=1)
        plt.rc("ytick.major", size=6, width=1)
        plt.rc("lines", markersize=10, markeredgewidth=2)
        plt.rc("lines", linewidth=3)

        try:
            if PlotPars.afe[0] != -1000:
                ylim = [PlotPars.afe[0], PlotPars.afe[1]]
            else:
                ylim = [mafe - abs(PlotPars.afe[1]),
                        mafe + abs(PlotPars.afe[1])]
        except:
            ylim = [mafe - 4 * eafe, mafe + 4 * eafe]

        panel_a = plt.subplot(311)
        plt.xlabel('EP = $\chi$ (eV)')
        plt.ylabel(ylabel)
        plt.xlim(-0.2, 5.2)
        plt.ylim(ylim)
        if hasattr(PlotPars, 'title_inside'):
            if PlotPars.title_inside != None:
                plt.text(plt.xlim()[0] + 0.50 * (plt.xlim()[1] - plt.xlim()[0]),
                         plt.ylim()[0] + 0.85 * (plt.ylim()[1] - plt.ylim()[0]),
                         PlotPars.title_inside,
                         horizontalalignment='center',
                         size=15)
        panel_b = plt.subplot(312)
        plt.xlabel('REW = log (EW/$\lambda$)')
        plt.ylabel(ylabel)
        plt.xlim(1.02 * min(list(rew1) + list(rew2)),
                 0.98 * max(list(rew1) + list(rew2)))
        plt.ylim(ylim)

        panel_c = plt.subplot(313)
        plt.xlabel('Wavelength ($\AA$)')
        plt.ylabel(ylabel)
        try:
            plt.xlim(PlotPars.wavelength_range[0], PlotPars.wavelength_range[1])
        except:
            plt.xlim(4100, 7900)
        plt.ylim(ylim)

        panel_a.plot(ep1, afe1, 'b+')
        panel_a.plot(ep2, afe2, 'go')
        panel_a.plot(x_epfit, y_epfit, 'b')

        panel_b.plot(rew1, afe1, 'b+')
        panel_b.plot(rew2, afe2, 'go')
        panel_b.plot(x_rewfit, y_rewfit, 'b')

        panel_c.plot(w1, afe1, 'b+')
        panel_c.plot(w2, afe2, 'go')
        panel_c.plot([4000, 8000], [mafe, mafe], 'black')

        if hasattr(PlotPars, 'directory'):
            if not os.path.exists(PlotPars.directory):
                os.mkdir(PlotPars.directory)
            plot = PlotPars.directory + '/' + plot
        if hasattr(PlotPars, 'figure_format'):
            plot = plot + '.' + PlotPars.figure_format
        plt.savefig(plot, bbox_inches='tight')
        # plt.close()

    if hasattr(Ref, 'name'):
        ref_star = Ref.name
    else:
        ref_star = None

    dfe = mfe1 - mfe2
    edfe = np.sqrt(efe1 ** 2 / nfe1 + efe2 ** 2 / nfe2)

    x = {'afe': round(mafe, 3), 'err_afe': round(eafe, 3),
         'afe1': round(mfe1, 3), 'err_afe1': round(efe1, 3), 'nfe1': nfe1,
         'afe2': round(mfe2, 3), 'err_afe2': round(efe2, 3), 'nfe2': nfe2,
         'slope_ep': slope_ep,
         'zero_ep': zero_ep,  # Added by LM
         'err_slope_ep': err_slope_ep,
         'slope_rew': slope_rew,
         'zero_rew': zero_rew,  # Added by LM
         'err_slope_rew': err_slope_rew,
         'reference': ref_star}
    Star.iron_stats = x

    if not silent:
        print("FeI  : {0:6.3f} +/- {1:5.3f} (n={2:3.0f})". \
              format(mfe1, efe1, nfe1))
        print("FeII : {0:6.3f} +/- {1:5.3f} (n={2:3.0f})". \
              format(mfe2, efe2, nfe2))

    return True


def solve_one(Star, SolveParsInit, Ref=object, PlotPars=object):
    sp = SolvePars()
    sp.__dict__ = SolveParsInit.__dict__.copy()
    if not hasattr(Star, 'model_atmosphere_grid'):
        logger.info('Star has no model yet. Calculating.')
        Star.get_model_atmosphere(sp.grid)
    if not hasattr(Star, 'model_atmosphere'):
        print('Unable to find a starting model atmosphere for this star')
        return None
    if Star.model_atmosphere_grid != sp.grid:
        logger.info('Inconsistent model atmosphere grids ' +
                    '(Star and SolvePars). ' +
                    'Fixing problem now.')
        Star.get_model_atmosphere(sp.grid)

    if hasattr(Ref, 'name'):
        if not hasattr(Ref, 'model_atmosphere_grid'):
            logger.info('Ref star has no model yet. Calculating.')
            Ref.get_model_atmosphere(sp.grid)
        if Ref.model_atmosphere_grid != sp.grid:
            logger.info('Inconsistent model atmosphere grids ' +
                        '(Ref star and SolvePars). ' +
                        'Fixing problem now.')
            Ref.get_model_atmosphere(sp.grid)

    dtv, dgv, dvv, stop_iter = [], [], [], False
    if hasattr(Star, 'converged'):
        if not Star.converged:
            Star.converged = False
    else:
        Star.converged = False
    Star.stop_iter = sp.niter
    if sp.niter == 0:
        Star.converged = True

    print('it Teff logg [Fe/H]  vt           [Fe/H]')
    print('-- ---- ---- ------ ----      --------------')

    for i in range(sp.niter + 1):
        if sp.step_teff <= 1 and sp.step_logg <= 0.01 \
                and sp.step_vt <= 0.01:
            if not stop_iter:
                Star.converged = False
                if SolveParsInit.niter > 0:
                    print('-- Begin final loop')
            stop_iter = True

        if i > 0:
            if Star.iron_stats['slope_ep'] > 0:
                Star.teff += sp.step_teff
            else:
                Star.teff -= sp.step_teff
            if Star.teff > 7000:
                Star.teff = 7000
            if Star.iron_stats['slope_rew'] > 0:
                Star.vt += sp.step_vt
            else:
                Star.vt -= sp.step_vt
            if Star.vt < 0:
                Star.vt = 0
            dfe = Star.iron_stats['afe1'] - Star.iron_stats['afe2']
            if dfe > 0:
                Star.logg += sp.step_logg
            else:
                Star.logg -= sp.step_logg
            if Star.logg > 5.0:
                Star.logg = 5.0

            if hasattr(Ref, 'name'):
                Star.feh = Ref.feh + Star.iron_stats['afe']
            else:
                Star.feh = Star.iron_stats['afe'] - sp.solar_afe
            if Star.feh > 1.0:
                Star.feh = 1.0
            if Star.feh > 0.5 and sp.grid != 'over':
                Star.feh = 0.5

            Star.get_model_atmosphere(sp.grid)

        if i + 1 == sp.niter or sp.niter == 0:
            plot = Star.name
            if hasattr(Ref, 'name'):
                plot = Star.name + '-' + Ref.name
                if Star.name == Ref.name:
                    plot = None
                    Star.converged = ''
        else:
            plot = None

        is_done = iron_stats(Star, Ref=Ref, plot=plot, PlotPars=PlotPars)

        print("{0:2.0f} {1:4.0f} {2:4.2f} {3:6.3f} {4:4.2f}" \
              " ---> {5:6.3f}+/-{6:5.3f}". \
              format(i, Star.teff, Star.logg, Star.feh, Star.vt,
                     Star.iron_stats['afe'], Star.iron_stats['err_afe']))

        dtv.append(Star.teff)
        dgv.append(Star.logg)
        dvv.append(Star.vt)

        if i >= 4:
            if np.std(dtv[-5:]) <= 0.8 * sp.step_teff and \
                            np.std(dgv[-5:]) <= 0.8 * sp.step_logg and \
                            np.std(dvv[-5:]) <= 0.8 * sp.step_vt:
                print('-- Converged at iteration ' + str(i) + \
                      ' of ' + str(sp.niter))
                if stop_iter:
                    plot = Star.name
                    if hasattr(Ref, 'name'):
                        plot = Star.name + '-' + Ref.name
                    iron_stats(Star, Ref=Ref, plot=plot, PlotPars=PlotPars)
                    Star.converged = True
                    Star.stop_iter = i
                    break
                sp.step_teff = sp.step_teff / 2
                sp.step_logg = sp.step_logg / 2
                sp.step_vt = sp.step_vt / 2
                if sp.step_teff < 1 and sp.step_teff > 0:
                    sp.step_teff = 1
                if sp.step_logg < 0.01 and sp.step_logg > 0:
                    sp.step_logg = 0.01
                if sp.step_vt < 0.01 and sp.step_vt > 0:
                    sp.step_vt = 0.01

    if not Star.converged:
        if hasattr(Ref, 'name'):
            if Star.name == Ref.name or SolveParsInit.niter == 0:
                print('--')
            else:
                print('-- Did not achieve final convergence.')
        else:
            print('-- Did not achieve final convergence.')

    print('------------------------------------------------------')

    if hasattr(Ref, 'name'):
        print('   D[Fe/H]    ||    D[Fe/H] Fe I   |   D[Fe/H] Fe II')
    else:
        print('    A(Fe)     ||      A(Fe I)      |     A(Fe II)   ')

    print("{0:6.3f} {1:6.3f} || {2:6.3f} {3:6.3f} {4:3d} " \
          "| {5:6.3f} {6:6.3f} {7:3d}". \
          format(Star.iron_stats['afe'], Star.iron_stats['err_afe'],
                 Star.iron_stats['afe1'], Star.iron_stats['err_afe1'],
                 Star.iron_stats['nfe1'],
                 Star.iron_stats['afe2'], Star.iron_stats['err_afe2'],
                 Star.iron_stats['nfe2']))
    print('------------------------------------------------------')

    Star.sp_err = {'teff': 0, 'logg': 0, 'afe': 0, 'vt': 0}
    if ((Star.converged and sp.errors) or \
                (sp.niter == 0 and sp.errors and Star.converged != '')):
        errors.error_one(Star, sp, Ref)
        Star.err_teff = int(Star.sp_err['teff'])
        Star.err_logg = Star.sp_err['logg']
        Star.err_feh = Star.sp_err['afe']
        Star.err_vt = Star.sp_err['vt']
        print("Solution with formal errors:")
        print("Teff    = {0:6d} +/- {1:5d}". \
              format(int(Star.teff), int(Star.sp_err['teff'])))
        print("log g   = {0:6.3f} +/- {1:5.3f}". \
              format(Star.logg, Star.sp_err['logg']))
        if hasattr(Ref, 'name'):
            print("D[Fe/H] = {0:6.3f} +/- {1:5.3f}". \
                  format(Star.iron_stats['afe'], Star.sp_err['afe']))
        else:
            print("A(Fe)   = {0:6.3f} +/- {1:5.3f}". \
                  format(Star.iron_stats['afe'], Star.sp_err['afe']))
        print("vt      = {0:6.2f} +/- {1:5.2f}". \
              format(Star.vt, Star.sp_err['vt']))
        print('------------------------------------------------------')


def solve_one_variant(Star, SolveParsInit, Ref=object, PlotPars=object):
    """ different optimization algorithm introduced by LMalavolta
    """
    sp = SolvePars()
    sp.__dict__ = SolveParsInit.__dict__.copy()
    if not hasattr(Star, 'model_atmosphere_grid'):
        logger.info('Star has no model yet. Calculating.')
        Star.get_model_atmosphere(sp.grid)
    if not hasattr(Star, 'model_atmosphere'):
        print('Unable to find a starting model atmosphere for this star')
        return None
    if Star.model_atmosphere_grid != sp.grid:
        logger.info('Inconsistent model atmosphere grids ' +
                    '(Star and SolvePars). ' +
                    'Fixing problem now.')
        Star.get_model_atmosphere(sp.grid)

    if hasattr(Ref, 'name'):
        if not hasattr(Ref, 'model_atmosphere_grid'):
            logger.info('Ref star has no model yet. Calculating.')
            Ref.get_model_atmosphere(sp.grid)
        if Ref.model_atmosphere_grid != sp.grid:
            logger.info('Inconsistent model atmosphere grids ' +
                        '(Ref star and SolvePars). ' +
                        'Fixing problem now.')
            Ref.get_model_atmosphere(sp.grid)

    dtv, dgv, dvv, stop_iter = [], [], [], False
    if hasattr(Star, 'converged'):
        if not Star.converged:
            Star.converged = False
    else:
        Star.converged = False
    Star.stop_iter = sp.niter
    if sp.niter == 0:
        Star.converged = True

    print('it Teff logg [Fe/H]  vt           [Fe/H]')
    print('-- ---- ---- ------ ----      --------------')

    # LM: no difference with the original algorithm until now
    plot = None
    is_done = iron_stats(Star, Ref=Ref, plot=plot, PlotPars=PlotPars)
    Star.sp_err = {'teff': 0, 'logg': 0, 'afe': 0, 'vt': 0}
    errors.error_one(Star, sp, Ref)

    ## Here: pseudo code in fortran
    for iQP in xrange(0, len(sp.QP)):

        plot = None

        # rejection of outliers
        afe1_model = Star.iron_stats['zero_ep'] + Star.iron_stats['slope_ep'] * Star.fe1['ep']
        afe2_model = Star.iron_stats['zero_ep'] + Star.iron_stats['slope_ep'] * Star.fe2['ep']
        Star.fe1['flag'] = (
        np.abs(afe1_model - Star.fe1['ab']) < sp.QP_sigma_rejection[iQP] * Star.iron_stats['err_afe1'])
        Star.fe2['flag'] = (
        np.abs(afe2_model - Star.fe2['ab']) < sp.QP_sigma_rejection[iQP] * Star.iron_stats['err_afe2'])

        # teff = teff_out ;  logg = logg_out ;  vmic = vmic_out ;  gfeh = gfeh_out
        teff_cycle_status = True
        vmic_cycle_status = True
        logg_cycle_status = True
        gfeh_cycle_status = True

        # ??
        # teff_var = teff_out
        # logg_var = logg_out
        # vmic_var = vmic_out
        # gfeh_var = gfeh_out
        Star.converged = True

        count_iters = 0
        count_recurs = 0

        while (teff_cycle_status or vmic_cycle_status or logg_cycle_status or gfeh_cycle_status):

            Star_copy = copy.deepcopy(Star)

            # if teff_fix then
            # teff_var = teff_inp
            # teff_err1 = 30.0
            # teff_cycle_status = .false.
            # else

            teff_range = np.arange(
                    max(sp.teff_min, Star.teff - sp.step_teff_alt * sp.inter_halfsteps),
                    min(sp.teff_max, Star.teff + sp.step_teff_alt * (sp.inter_halfsteps + 1)),
                        sp.step_teff_alt)

            interp_coeff = np.zeros(len(teff_range), dtype=np.double)
            sel_done = np.zeros(len(teff_range), dtype=bool)

            for i_teff, v_teff in enumerate(teff_range):
                Star_copy.teff = v_teff
                Star_copy.get_model_atmosphere(sp.grid)
                is_done = iron_stats(Star_copy, Ref=Ref, plot=plot, PlotPars=PlotPars)
                if is_done:
                    interp_coeff[i_teff] = Star_copy.iron_stats['slope_ep']
                    sel_done[i_teff] = True

            ind_zero = np.argmin(np.abs(interp_coeff))

            # LM To improve the parameter determination, only the points close to the zero
            # LM of the slope curve are selected
            sel_teff = sel_done & (np.abs(interp_coeff - interp_coeff[ind_zero]) < sp.step_teff_alt * (sp.iter_poly_order + 1))

            if np.sum(sel_teff) <= sp.iter_poly_order + 1:
                teff_iter = Star.teff
            else:
                poly_coeff = np.polyfit(interp_coeff[sel_teff], teff_range[sel_teff], sp.iter_poly_order)
                teff_iter = min(max(np.polyval(poly_coeff, 0.), sp.teff_min), sp.teff_max)

            # Temperature of the Star copy object is set to the value at the start of the iteration
            Star_copy.teff = Star.teff

            # Same algorithm applied to microturbulent velocity
            vmic_range = np.arange (
                    max(sp.vmic_min, Star.vt - sp.step_vmic_alt * sp.inter_halfsteps),
                    min(sp.vmic_max, Star.vt + sp.step_vmic_alt * (sp.inter_halfsteps + 1)),
                        sp.step_vmic_alt)

            interp_coeff = np.zeros(len(vmic_range), dtype=np.double)
            sel_done = np.zeros(len(vmic_range), dtype=bool)

            for i_vmic, v_vmic in enumerate(vmic_range):
                Star_copy.vt = v_vmic
                Star_copy.get_model_atmosphere(sp.grid)
                is_done = iron_stats(Star_copy, Ref=Ref, plot=plot, PlotPars=PlotPars)
                if is_done:
                    interp_coeff[i_vmic] = Star_copy.iron_stats['slope_ep']
                    sel_done[i_vmic] = True

            ind_zero = np.argmin(np.abs(interp_coeff))

            # LM To improve the parameter determination, only the points close to the zero
            # LM of the slope curve are selected
            sel_vmic = sel_done & (np.abs(interp_coeff - interp_coeff[ind_zero]) < sp.step_vmic_alt * (sp.iter_poly_order + 1))

            if np.sum(sel_vmic) <= sp.iter_poly_order + 1:
                vmic_iter = Star.vt
            else:
                poly_coeff = np.polyfit(interp_coeff[sel_vmic], vmic_range[sel_vmic], sp.iter_poly_order)
                vmic_iter = min(max(np.polyval(poly_coeff, 0.), sp.vmic_min), sp.vmic_max)

            # Microturbulent velocity of the Star copy object is set to the value at the start of the iteration
            Star_copy.vt = Star.vt

            # Same algorithm applied to gravity
            logg_range = np.arange (
                    max(sp.logg_min, Star.logg - sp.step_logg_alt * sp.inter_halfsteps),
                    min(sp.logg_max, Star.logg + sp.step_logg_alt * (sp.inter_halfsteps + 1)),
                        sp.step_logg_alt)

            interp_coeff = np.zeros(len(logg_range), dtype=np.double)
            sel_done = np.zeros(len(logg_range), dtype=bool)

            for i_logg, v_logg in enumerate(logg_range):
                Star_copy.logg = v_logg
                Star_copy.get_model_atmosphere(sp.grid)
                is_done = iron_stats(Star_copy, Ref=Ref, plot=plot, PlotPars=PlotPars)
                if is_done:
                    interp_coeff[i_logg] = Star_copy.iron_stats['afe2'] - Star.iron_stats['afe1']
                    sel_done[i_logg] = True

            ind_zero = np.argmin(np.abs(interp_coeff))

            # LM To improve the parameter determination, only the points close to the zero
            # LM of the slope curve are selected
            sel_logg = sel_done & (np.abs(interp_coeff - interp_coeff[ind_zero]) < sp.step_logg_alt * (sp.iter_poly_order + 1))

            if np.sum(sel_logg) <= sp.iter_poly_order + 1:
                logg_iter = Star.logg
            else:
                poly_coeff = np.polyfit(interp_coeff[sel_logg], logg_range[sel_logg], sp.iter_poly_order)
                logg_iter = min(max(np.polyval(poly_coeff, 0.), sp.logg_min), sp.logg_max)

            # Get the new abundance with all the other photospheric parameters varied to the new values
            Star_copy.teff = teff_iter
            Star_copy.logg = logg_iter
            Star_copy.vt   = vmic_iter

            Star_copy.get_model_atmosphere(sp.grid)
            is_done = iron_stats(Star_copy, Ref=Ref, plot=plot, PlotPars=PlotPars)
            gfeh_iter = Star_copy.feh

            # Check how far the new values are with respect to the initial values while keeping into account the associated errors
            teff_crit = min(max(Star.sp_err['teff'] / sp.QP[iQP], sp.teff_crit[0]), sp.teff_crit[1])
            logg_crit = min(max(Star.sp_err['logg'] / sp.QP[iQP], sp.logg_crit[0]), sp.logg_crit[1])
            vmic_crit = min(max(Star.sp_err['vt'] / sp.QP[iQP], sp.vmic_crit[0]), sp.vmic_crit[1])
            gfeh_crit = min(max(Star.sp_err['afe'] / sp.QP[iQP], sp.gfeh_crit[0]), sp.gfeh_crit[1])

            if np.abs(teff_iter - Star.teff) < teff_crit and sp.teff_min < teff_iter < sp.teff_max: teff_cycle_status = False
            if np.abs(logg_iter - Star.logg) < logg_crit and sp.logg_min < logg_iter < sp.logg_max: logg_cycle_status = False
            if np.abs(vmic_iter - Star.vt  ) < vmic_crit and sp.vmic_min < vmic_iter < sp.vmic_max: vmic_cycle_status = False
            if np.abs(gfeh_iter - Star.feh ) < gfeh_crit and sp.gfeh_min < gfeh_iter < sp.gfeh_max: gfeh_cycle_status = False

            # LM iteration output values are composed by 2/3 of the new value and 1/3 of the previous one, in a zig-zag
            # LM scheme to avoid degeneracies
            Star.teff = (Star.teff + 2. * teff_iter) / 3.
            Star.logg = (Star.logg + 2. * logg_iter) / 3.
            Star.vt   = (Star.vt   + 2. * vmic_iter) / 3.
            Star.feh  = (Star.feh  + 2. * gfeh_iter) / 3.

            Star.get_model_atmosphere(sp.grid)

            is_done = iron_stats(Star, Ref=Ref, plot=plot, PlotPars=PlotPars)
            errors.error_one(Star, sp, Ref)
            # Get the error on the new temperature

            if not (teff_cycle_status or vmic_cycle_status or logg_cycle_status or gfeh_cycle_status) \
                    and count_iters < sp.count_iters_limit:
                count_iters += 1
                teff_cycle_status = True
                vmic_cycle_status = True

            count_recurs += 1
            Star.stop_iter = count_recurs

            print("{0:2.0f} {1:4.0f} {2:4.2f} {3:6.3f} {4:4.2f}"
                " ---> {5:6.3f}+/-{6:5.3f}".
                format(count_recurs, Star.teff, Star.logg, Star.feh, Star.vt,
                        Star.iron_stats['afe'], Star.iron_stats['err_afe']))

            if count_recurs > sp.count_recurs_limit:
                print ' Excessive number of iterations, exiting from this iQP cycle '
                Star.converged = False
                break

        print
        print("{0:2.0f} {1:4.0f} {2:4.2f} {3:6.3f} {4:4.2f}" \
              " ---> {5:6.3f}+/-{6:5.3f}". \
              format(iQP, Star.teff, Star.logg, Star.feh, Star.vt,
                     Star.iron_stats['afe'], Star.iron_stats['err_afe']))
        print

    plot = Star.name
    if hasattr(Ref, 'name'):
        plot = Star.name + '-' + Ref.name
        if Star.name == Ref.name:
            plot = None
            Star.converged = ''
    iron_stats(Star, Ref=Ref, plot=plot, PlotPars=PlotPars)

    # LM back to original program
    if not Star.converged:
        if hasattr(Ref, 'name'):
            if Star.name == Ref.name or SolveParsInit.niter == 0:
                print('--')
            else:
                print('-- Did not achieve final convergence.')
        else:
            print('-- Did not achieve final convergence.')

    print('------------------------------------------------------')

    if hasattr(Ref, 'name'):
        print('   D[Fe/H]    ||    D[Fe/H] Fe I   |   D[Fe/H] Fe II')
    else:
        print('    A(Fe)     ||      A(Fe I)      |     A(Fe II)   ')

    print("{0:6.3f} {1:6.3f} || {2:6.3f} {3:6.3f} {4:3d} " \
          "| {5:6.3f} {6:6.3f} {7:3d}". \
          format(Star.iron_stats['afe'], Star.iron_stats['err_afe'],
                 Star.iron_stats['afe1'], Star.iron_stats['err_afe1'],
                 Star.iron_stats['nfe1'],
                 Star.iron_stats['afe2'], Star.iron_stats['err_afe2'],
                 Star.iron_stats['nfe2']))
    print('------------------------------------------------------')

    Star.sp_err = {'teff': 0, 'logg': 0, 'afe': 0, 'vt': 0}
    if ((Star.converged and sp.errors) or \
                (sp.niter == 0 and sp.errors and Star.converged != '')):
        errors.error_one(Star, sp, Ref)
        Star.err_teff = int(Star.sp_err['teff'])
        Star.err_logg = Star.sp_err['logg']
        Star.err_feh = Star.sp_err['afe']
        Star.err_vt = Star.sp_err['vt']
        print("Solution with formal errors:")
        print("Teff    = {0:6d} +/- {1:5d}". \
              format(int(Star.teff), int(Star.sp_err['teff'])))
        print("log g   = {0:6.3f} +/- {1:5.3f}". \
              format(Star.logg, Star.sp_err['logg']))
        if hasattr(Ref, 'name'):
            print("D[Fe/H] = {0:6.3f} +/- {1:5.3f}". \
                  format(Star.iron_stats['afe'], Star.sp_err['afe']))
        else:
            print("A(Fe)   = {0:6.3f} +/- {1:5.3f}". \
                  format(Star.iron_stats['afe'], Star.sp_err['afe']))
        print("vt      = {0:6.2f} +/- {1:5.2f}". \
              format(Star.vt, Star.sp_err['vt']))
        print('------------------------------------------------------')


def solve_all(Data, SolveParsInit, output_file, reference_star=None,
              PlotPars=object, ultra_verbose=False):  # LM Added ultra_verbose
    print('------------------------------------------------------')
    print('Initializing ...')
    start_time = datetime.datetime.now()
    print('- Date and time: ' + start_time.strftime('%d-%b-%Y, %H:%M:%S'))
    print('- Model atmospheres: ' + SolveParsInit.grid)
    print('- Star data: ' + Data.star_data_fname)
    print('- Line list: ' + Data.lines_fname)
    print('------------------------------------------------------')
    if reference_star:
        Ref = Star(reference_star)
        Ref.get_data_from(Data)
    else:
        Ref = None
    fout = open(output_file, 'w')
    if SolveParsInit.errors:
        fout.write('id,teff,logg,feh_model,vt,feh,err_feh_,' +
                   'feh1,err_feh1,nfe1,feh2,err_feh2,nfe2,'
                   'slope_ep,err_slope_ep,slope_rew,err_slope_rew,'
                   'stop_iter,converged,'
                   'err_teff,err_logg,err_feh,err_vt\n')
    else:
        fout.write('id,teff,logg,feh_model,vt,feh,err_feh,' +
                   'feh1,err_feh1,nfe1,feh2,err_feh2,nfe2,'
                   'slope_ep,err_slope_ep,slope_rew,err_slope_rew,'
                   'stop_iter,converged,'
                   'err_teff,err_logg,err_feh_,err_vt\n')
    for star_id in Data.star_data['id']:
        print('')
        print('*' * len(star_id))
        print(star_id)
        print('*' * len(star_id))
        s = Star(star_id)
        try:
            s.get_data_from(Data)
        except:
            logger.warning('No data found for ' + s.name + \
                           '. Excluded from output file.')
            print('Data not found.')
            # fout.write("{0},,,,,,,,,,"\
            #           ",,,,,,,,,,,,\n".\
            #           format(s.name))
            continue
        if ma.count(Data.lines[star_id]) == 0:
            print('Line data not found.')
            continue
        sp = SolvePars()
        sp.__dict__ = SolveParsInit.__dict__.copy()
        if reference_star:
            if s.name == Ref.name:
                sp.niter = 0
                print('Reference star. No calculations needed.')
                # continue
        if hasattr(s, 'converged') and sp.check_converged:
            if s.converged == 'True':
                print('Already converged.')
                continue
                # sp.niter = 0
                # s.converged = True
        if s.name in sp.ignore:
            print('Asked to ignore.')
            continue

        solve_one(s, sp, Ref, PlotPars=PlotPars)

        if sp.niter == 0:
            s.converged = ''

        fout.write("{0},{1:4.0f},{2:5.3f},{3},{4:4.2f},{5},{6:5.3f}," \
                   "{7},{8:5.3f},{9}," \
                   "{10},{11:5.3f},{12},{13:.6f},{14:.6f}," \
                   "{15:.6f},{16:.6f},{17},{18}," \
                   "{19:3.0f},{20:5.3f},{21:5.3f},{22:4.2f}\n". \
                   format(s.name, s.teff, s.logg, str(round(s.feh, 3)), s.vt,
                          str(round(s.iron_stats['afe'], 3)),
                          s.iron_stats['err_afe'],
                          str(round(s.iron_stats['afe1'], 3)),
                          s.iron_stats['err_afe1'],
                          s.iron_stats['nfe1'],
                          str(round(s.iron_stats['afe2'], 3)),
                          s.iron_stats['err_afe2'],
                          s.iron_stats['nfe2'],
                          s.iron_stats['slope_ep'],
                          s.iron_stats['err_slope_ep'],
                          s.iron_stats['slope_rew'],
                          s.iron_stats['err_slope_rew'],
                          s.stop_iter,
                          s.converged,
                          s.sp_err['teff'], s.sp_err['logg'],
                          s.sp_err['afe'], s.sp_err['vt']
                          ))

        if ultra_verbose:
            print "[Fe/H](Fe I)  = {0:5.3f} +/- {1:5.3f}". \
                format(s.iron_stats['afe1'], s.iron_stats['err_afe1'])
            print "[Fe/H](Fe II) = {0:5.3f} +/- {1:5.3f}". \
                format(s.iron_stats['afe2'], s.iron_stats['err_afe2'])
            print "A(FeI) vs. EP slope  = {0:.6f}".format(s.iron_stats['slope_ep'])
            print "A(FeI) vs. REW slope = {0:.6f}".format(s.iron_stats['slope_rew'])

            print "Final stellar parameters:"
            print "Teff = {0:4.2f} +/- {1:4.2f} K".format(s.teff, s.sp_err['teff'])
            print "logg = {0:4.2f} +/- {1:4.2f} dex".format(s.logg, s.sp_err['logg'])
            print " vt = {0:4.2f} +/- {1:4.2f} km/s".format(s.vt, s.sp_err['vt'])
            print "[Fe/H]= {0:5.2f} +/- {1:4.2f} dex".format(s.feh, s.sp_err['afe'])

    fout.close()

    print('')
    print('------------------------------------------------------')
    end_time = datetime.datetime.now()
    print('- Date and time: ' + end_time.strftime('%d-%b-%Y, %H:%M:%S'))
    delta_t = (end_time - start_time).seconds
    hours, remainder = divmod(delta_t, 3600)
    minutes, seconds = divmod(remainder, 60)
    print('- Time elapsed: %sH %sM %sS' % (hours, minutes, seconds))
    print('Done!')
    print('------------------------------------------------------')
    print('')


def make_single_solution_table(solution_files, single_solution_file):
    """Takes q2.specpars.solve_all outputs and creates a single final one

    Files must be in the order in which they were computed!
    """
    # solution_files = ['starsDec_solution1.csv', 'starsDec_solution2.csv']
    # single_solution_file = 'starsDec_solution.csv'
    fout = open(single_solution_file, 'w')
    with open(solution_files[0], 'r') as f:
        lines = f.readlines()
    for line in lines:
        sid = line[0:line.index(',')]
        if 'True' in line or 'id,teff' in line:
            # nline = line[0:line.rfind('\n')]
            # linew = nline[0:nline.rfind('\n')]
            # fout.write(linew+'\n')
            # print line
            fout.write(line)
        else:
            for i in range(1, len(solution_files)):
                with open(solution_files[i], 'r') as f2:
                    lines2 = f2.readlines()
                for line2 in lines2:
                    sid2 = line2[0:line2.index(',')]
                    # nline2 = line2[0:line2.rfind(',')]
                    # line2w = nline2[0:nline2.rfind(',')]
                    if 'True' in line2 and sid == sid2:
                        # fout.write(line2w+'\n')
                        fout.write(line2)
    fout.close()


def fancy_ironstats_plot(Star):
    """Makes bokeh hover-ing plots

    Function written to look for outliers and investigate line-to-line scatter
    """
    if not hasattr(Star, 'iron_stats'):
        logger.error('Star object (' + Star.name + ') has no ironstats attribute.')
        return None
    ww = np.concatenate((Star.fe1['ww'], Star.fe2['ww']))
    ep = np.concatenate((Star.fe1['ep'], Star.fe2['ep']))
    ew = np.concatenate((Star.fe1['ew'], Star.fe2['ew']))
    rew = np.concatenate((Star.fe1['rew'], Star.fe2['rew']))
    if Star.iron_stats['reference']:
        ab = np.concatenate((Star.fe1['difab'], Star.fe2['difab']))
        y_axis_label = '[Fe/H]'
    else:
        ab = np.concatenate((Star.fe1['ab'], Star.fe2['ab']))
        y_axis_label = 'A(Fe)'
    ws = [str(round(w, 1)) for w in ww]
    colors = np.concatenate((["blue"] * len(Star.fe1['ww']),
                             ["green"] * len(Star.fe2['ww'])))

    TOOLS = "pan,wheel_zoom,box_zoom,reset,hover"
    output_notebook()

    title = Star.name
    if getattr(Star, 'iron_stats')['reference']:
        title += ' - ' + getattr(Star, 'iron_stats')['reference']

    p1 = figure(title=title, plot_width=650, plot_height=300,
                x_axis_label='EP (eV)',
                y_axis_label=y_axis_label,
                tools=TOOLS)

    abst = [str(round(xab, 3)) for xab in ab]
    source = ColumnDataSource(
        data=dict(
            ws=ws,
            ep=ep,
            rew=rew,
            ab=ab,
            abst=abst,
            ew=ew,
            colors=colors,
        )
    )

    p1.scatter('ep', 'ab', size=10, color='colors',
               source=source, marker='square')

    hover = p1.select(dict(type=HoverTool))
    hover.tooltips = OrderedDict([
        ("Wavelength, EP", "@ws A, @ep eV"),
        ("EW, REW", "@ew mA, @rew"),
        ("Abundance", "@abst"),
    ])

    show(p1)

    p2 = figure(title='', plot_width=650, plot_height=300,
                x_axis_label='REW',
                y_axis_label=y_axis_label,
                tools=TOOLS)

    p2.scatter('rew', 'ab', size=10, color='colors',
               source=source, marker='square')

    hover = p2.select(dict(type=HoverTool))
    hover.tooltips = OrderedDict([
        ("Wavelength, EP", "@ws A, @ep eV"),
        ("EW, REW", "@ew mA, @rew"),
        ("Abundance", "@abst"),
    ])

    show(p2)
