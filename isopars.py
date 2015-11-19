import sqlite3
import numpy as np
import logging
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.interpolate import griddata
from config import *
import os
import datetime
from star import Star


logger = logging.getLogger(__name__)

class SolvePars:
    def __init__(self, key_parameter_known='plx',
                 db='yy02.sql3', feh_offset = 0,
                 nsigma=5, window_len_age=13):
        self.key_parameter_known = key_parameter_known
        self.get_isochrone_points_db = db
        self.feh_offset = feh_offset
        self.get_isochrone_points_nsigma = nsigma
        self.smooth_window_len_age = window_len_age
        self.smooth_window_len_mass = 0
        self.smooth_window_len_logl = 0
        self.smooth_window_len_mv = 0
        self.smooth_window_len_r = 0
        self.smooth_window_len_logg = 0

class PlotPars:
    def __init__(self, figure_format='png', directory="", make_figures=True):
        self.age_xlim = [0, 14]
        self.mass_xlim = None
        self.logl_xlim = None
        self.mv_xlim = None
        self.r_xlim = None
        self.logg_xlim = None
        self.directory = directory
        self.figure_format = figure_format
        self.title_inside = None
        self.make_figures = make_figures

def pdf(pdf_x, ips, prob, par, smooth_window_len):
    '''Calculates a probability distribution function (PDF) for parameter par
    given the x-values for the PDF, the isochrone points ips, and their 
    probability. Return PDF and smoothed PDF (using smooth_window_len) if
    possible (otherwise returns two non-smoothed PDFs), as well as a stats
    dictionary with mean, std, most probable value, etc.
    '''
    dx = 0.5*(pdf_x[1] - pdf_x[0])
    pdf_y = []
    for x in pdf_x:
        pdf_y.append(sum(prob[np.logical_and(ips[par] >= x-dx,
                                             ips[par] <  x+dx)]))
    pdf_y = np.array(pdf_y)
    pdf_y = pdf_y/simps(pdf_y, pdf_x)

    try:
        pdf_y_smooth = smooth(pdf_y, smooth_window_len)
        pdf_y_smooth = pdf_y_smooth/simps(pdf_y_smooth, pdf_x)
    except:
        pdf_y_smooth = pdf_y
        logger.warning('Unable to smooth '+par+' PDF.')

    stats = get_stats(pdf_x, pdf_y_smooth)

    if stats['most_probable']:
        print "{0:10s} = {1:6.3f} | {2:6.3f} - {3:6.3f} | "\
              "{4:6.3f} - {5:6.3f} | {6:6.3f} +/- {7:6.3f}"\
              .format(par,
                      stats['most_probable'],
                      stats['lower_limit_1sigma'],
                      stats['upper_limit_1sigma'],
                      stats['lower_limit_2sigma'],
                      stats['upper_limit_2sigma'],
                      stats['mean'], stats['std'])
    else:
        print "{0:10s} =        |        -        |  "\
              "      -        | {1:6.3f} +/- {2:6.3f}"\
              .format(par, stats['mean'], stats['std'])
        logger.warning("--- Unable to calculate PDF stats for "+par)

    return pdf_y, pdf_y_smooth, stats

def get_stats(pdf_x, pdf_y_smooth):
    stats = {}
    stats['most_probable'] = \
      np.mean(np.array(pdf_x)[pdf_y_smooth == max(pdf_y_smooth)])
    stats['mean'] = simps(pdf_y_smooth*pdf_x, pdf_x)
    stats['std'] = np.sqrt(simps(pdf_y_smooth*(pdf_x-stats['mean'])**2,\
                                 pdf_x))

    k = pdf_x <= stats['most_probable']
    pdf_y_left = 0.5*pdf_y_smooth[k]/simps(pdf_y_smooth[k], pdf_x[k])
    pdf_x_left = pdf_x[k]
    areas_left = []
    for x in pdf_x_left:
        areas_left.append(simps(pdf_y_left[pdf_x_left <= x],
                                pdf_x_left[pdf_x_left <= x]))
    areas_left = np.array(areas_left)
    if np.mean(areas_left) == 0:
        logger.warning("Left side of distribution is empty")
        stats['most_probable'] = None
        stats['lower_limit_1sigma'] = None
        stats['lower_limit_2sigma'] = None
        stats['upper_limit_1sigma'] = None
        stats['upper_limit_2sigma'] = None
        return stats

    k = pdf_x >= stats['most_probable']
    pdf_y_right = 0.5*pdf_y_smooth[k]/simps(pdf_y_smooth[k], pdf_x[k])
    pdf_x_right = pdf_x[k]
    areas_right = []
    for x in pdf_x_right:
        areas_right.append(simps(pdf_y_right[pdf_x_right <= x],
                                 pdf_x_right[pdf_x_right <= x]))
    areas_right = np.array(areas_right)

    try:
        stats['lower_limit_1sigma'] = \
          np.mean(griddata(areas_left, pdf_x_left, 0.158))
        stats['lower_limit_2sigma'] = \
          np.mean(griddata(areas_left, pdf_x_left, 0.022))
        stats['upper_limit_1sigma'] = \
          np.mean(griddata(areas_right, pdf_x_right, 0.341))
        stats['upper_limit_2sigma'] = \
          np.mean(griddata(areas_right, pdf_x_right, 0.477))
    except:
        stats['lower_limit_1sigma'] = -9.999
        stats['lower_limit_2sigma'] = -9.999
        stats['upper_limit_1sigma'] = -9.999
        stats['upper_limit_2sigma'] = -9.999

    return stats

def solve_one(Star, SolvePars, PlotPars):
    '''Calculates most likely parameters of Star using isochrone points
    '''
    ips = get_isochrone_points(Star, SolvePars.feh_offset,
                               SolvePars.get_isochrone_points_db,
                               SolvePars.get_isochrone_points_nsigma,
                               SolvePars.key_parameter_known)
    if ips == None:
        logger.warning('Could not get any isochrone points.')
        return None
    print 'Using {0} isochrone points'.format(len(ips['age']))
    logger.info('Using {0} Y2 isochrone points'.format(len(ips['age'])))
    Star.isokeyparameterknown = SolvePars.key_parameter_known
    Star.isonpoints = len(ips['age'])
    ips['t'] = 10**ips['logt']
    ips['r'] = 10**(0.5*(np.log10(ips['mass'])-ips['logg']+4.437))

    if SolvePars.key_parameter_known == 'logg':
        prob = np.exp(-1*((ips['t']-Star.teff)/          \
                          (1.414214*Star.err_teff))**2)* \
               np.exp(-1*((ips['logg']-Star.logg)/       \
                          (1.414214*Star.err_logg))**2)* \
               np.exp(-1*((ips['feh']-Star.feh)/         \
                          (1.414214*Star.err_feh))**2)
    if SolvePars.key_parameter_known == 'plx':
        prob = np.exp(-1*((ips['t']-Star.teff)/          \
                          (1.414214*Star.err_teff))**2)* \
               np.exp(-1*((ips['mv']-Star.M_V)/       \
                          (1.414214*Star.err_M_V))**2)* \
               np.exp(-1*((ips['feh']-Star.feh)/         \
                          (1.414214*Star.err_feh))**2)

    #age
    ages = 0.1+np.arange(150)*0.1
    pdf_age_x = ages[np.logical_and(ages >= min(ips['age'])-0.2,
                                   ages <= max(ips['age'])+0.2)]
    pdf_age_y, pdf_age_y_smooth, Star.isoage = \
      pdf(pdf_age_x, ips, prob, 'age', SolvePars.smooth_window_len_age)
    Star.pdf_age = {'x': pdf_age_x, 'y': pdf_age_y, 'ys': pdf_age_y_smooth}

    #mass
    masses = 0.4+np.arange(211)*0.01
    pdf_mass_x = masses[np.logical_and(masses >= min(ips['mass'])-0.02,
                                       masses <= max(ips['mass'])+0.02)]
    pdf_mass_y, pdf_mass_y_smooth, Star.isomass = \
      pdf(pdf_mass_x, ips, prob, 'mass', SolvePars.smooth_window_len_mass)

    #luminosity
    logls = -1.0+np.arange(401)*0.01
    pdf_logl_x = logls[np.logical_and(logls >= min(ips['logl'])-0.02,
                                      logls <= max(ips['logl'])+0.02)]
    pdf_logl_y, pdf_logl_y_smooth, Star.isologl = \
      pdf(pdf_logl_x, ips, prob, 'logl', SolvePars.smooth_window_len_logl)

    #absolute magnitude
    mvs = -3.0+np.arange(1601)*0.01
    pdf_mv_x = mvs[np.logical_and(mvs >= min(ips['mv'])-0.02,
                                  mvs <= max(ips['mv'])+0.02)]
    pdf_mv_y, pdf_mv_y_smooth, Star.isomv = \
      pdf(pdf_mv_x, ips, prob, 'mv', SolvePars.smooth_window_len_mv)

    #radius
    rs = 0.4+np.arange(211)*0.01
    pdf_r_x = rs[np.logical_and(rs >= min(ips['r'])-0.02,
                                rs <= max(ips['r'])+0.02)]
    pdf_r_y, pdf_r_y_smooth, Star.isor = \
      pdf(pdf_r_x, ips, prob, 'r', SolvePars.smooth_window_len_r)

    #logg
    if SolvePars.key_parameter_known == 'plx':
        loggs = np.arange(501)*0.01
        pdf_logg_x = loggs[np.logical_and(loggs >= min(ips['logg'])-0.05,
                                          loggs <= max(ips['logg'])+0.05)]
        pdf_logg_y, pdf_logg_y_smooth, Star.isologg = \
          pdf(pdf_logg_x, ips, prob, 'logg', SolvePars.smooth_window_len_logg)

    if not PlotPars.make_figures:
        return

    if not os.path.exists(PlotPars.directory) and PlotPars.directory != "":
        os.mkdir(PlotPars.directory)

    if Star.isoage:
        plt.figure(figsize=(7, 4))
        plt.rc("axes", labelsize=15, titlesize=12)
        plt.rc("xtick", labelsize=14)
        plt.rc("ytick", labelsize=14)
        plt.rc("lines", markersize=10, markeredgewidth=2)
        plt.rc("lines", linewidth=2)
        plt.rc("xtick.major", size=6, width=1)
        plt.rc("ytick.major", size=6, width=1)
        plt.xlim([0,15])
        plt.xlabel('Age (Gyr)')
        if PlotPars.age_xlim:
            plt.xlim(PlotPars.age_xlim)
        plt.ylabel('Probability density')
        k2 = np.logical_and(pdf_age_x >= Star.isoage['lower_limit_2sigma'],
                            pdf_age_x <= Star.isoage['upper_limit_2sigma'])
        k1 = np.logical_and(pdf_age_x >= Star.isoage['lower_limit_1sigma'],
                            pdf_age_x <= Star.isoage['upper_limit_1sigma'])
        plt.fill_between(pdf_age_x[k2], 0 , pdf_age_y_smooth[k2],
                         color='0.8', hatch="/")
        plt.fill_between(pdf_age_x[k1], 0 , pdf_age_y_smooth[k1],
                         color='0.6', hatch="X")
        plt.plot([Star.isoage['most_probable'], Star.isoage['most_probable']],
                 [0, max(pdf_age_y_smooth)], 'g--')
        plt.plot(pdf_age_x, pdf_age_y_smooth, 'g')

        #exclude eveything after __ in Star.name in legend:
        starname = Star.name.split("__")[0]

        if PlotPars.title_inside != None:
            starname = PlotPars.title_inside

        plt.text(0.92*plt.xlim()[1], 0.86*plt.ylim()[1], starname,
                 horizontalalignment='right', size=16)
        fig_name = os.path.join(PlotPars.directory,
                                Star.name.replace(' ', '_')+\
                                '_isoage_'+SolvePars.key_parameter_known)
        plt.savefig(fig_name+'.'+PlotPars.figure_format, bbox_inches='tight')
        plt.close()

    plt.figure(figsize=(6, 14))
    plt.rc("axes", labelsize=15, titlesize=12)
    plt.rc("xtick", labelsize=14)
    plt.rc("ytick", labelsize=14)
    plt.rc("lines", markersize=10, markeredgewidth=2)
    plt.rc("lines", linewidth=2)
    plt.rc("xtick.major", size=6, width=1)
    plt.rc("ytick.major", size=6, width=1)
    plt.subplots_adjust(hspace=0.4)

    npanels = 5
    if SolvePars.key_parameter_known == 'plx':
        npanels += 1

    for panel in np.arange(npanels)+1:
        ax = plt.subplot(npanels, 1, panel)
        ax.get_yaxis().set_visible(False)
        if panel == 1:
            pdf_x, pdf_y, pdf_y_smooth = \
              pdf_age_x, pdf_age_y, pdf_age_y_smooth
            par = Star.isoage
            ax.set_xlabel('Age (Gyr)')
            if PlotPars.age_xlim:
                ax.set_xlim(PlotPars.age_xlim)
        if panel == 2:
            pdf_x, pdf_y, pdf_y_smooth = \
              pdf_mass_x, pdf_mass_y, pdf_mass_y_smooth
            par = Star.isomass
            ax.set_xlabel('Mass ($M_\odot$)')
            if PlotPars.mass_xlim:
                ax.set_xlim(PlotPars.mass_xlim)
        if panel == 3:
            pdf_x, pdf_y, pdf_y_smooth = \
              pdf_logl_x, pdf_logl_y, pdf_logl_y_smooth
            par = Star.isologl
            ax.set_xlabel('$\log\,(L/L_\odot)$')
            if PlotPars.logl_xlim:
                ax.set_xlim(PlotPars.logl_xlim)
        if panel == 4:
            pdf_x, pdf_y, pdf_y_smooth = \
              pdf_mv_x, pdf_mv_y, pdf_mv_y_smooth
            par = Star.isomv
            ax.set_xlabel('$M_V$')
            if PlotPars.mv_xlim:
                ax.set_xlim(PlotPars.mv_xlim)
        if panel == 5:
            pdf_x, pdf_y, pdf_y_smooth = \
              pdf_r_x, pdf_r_y, pdf_r_y_smooth
            par = Star.isor
            ax.set_xlabel('Radius ($R_\odot$)')
            if PlotPars.r_xlim:
                ax.set_xlim(PlotPars.r_xlim)
        if panel == 6:
            pdf_x, pdf_y, pdf_y_smooth = \
              pdf_logg_x, pdf_logg_y, pdf_logg_y_smooth
            par = Star.isologg
            ax.set_xlabel('$\log g$ [cgs]')
            if PlotPars.logg_xlim:
                ax.set_xlim(PlotPars.logg_xlim)
        ax.plot(pdf_x, pdf_y, color='0.8')
        if par:
            ax.plot([par['most_probable'], par['most_probable']],
                    [0, max(pdf_y_smooth)], 'g--')
        ax.plot(pdf_x, pdf_y_smooth, 'g')

    fig_name = os.path.join(PlotPars.directory,
                            Star.name.replace(' ', '_')+\
                            '_isopar_'+SolvePars.key_parameter_known)
    plt.savefig(fig_name+'.'+PlotPars.figure_format, bbox_inches='tight')
    plt.close()

def solve_all(Data, SolvePars, PlotPars, output_file):
    print '------------------------------------------------------'
    print 'Initializing ...'
    start_time = datetime.datetime.now()
    print '- Date and time: '+start_time.strftime('%d-%b-%Y, %H:%M:%S')
    print '- Star data: '+Data.star_data_fname
    print '------------------------------------------------------'
    fout = open(output_file, 'wb')
    pars = ['age', 'mass', 'logl', 'mv', 'r']
    if SolvePars.key_parameter_known == 'plx':
        pars.append('logg')
    values = ['mp', 'll1s', 'ul1s', 'll2s', 'ul2s', 'mean', 'std']
    hd = 'id'
    for par in pars:
        for value in values:
            hd += ','+par+'_'+value
    fout.write(hd+'\n')
    for star_id in Data.star_data['id']:
        print ''
        print '*'*len(star_id)
        print star_id
        print '*'*len(star_id)
        s = Star(star_id)
        s.get_data_from(Data)
        try:
            solve_one(s, SolvePars, PlotPars)
        except:
            print 'Unable to find isochrone parameters.'
            print 'Input data might be missing or are too far from valid '+\
                  'isochrone points.'
            string = "{0}".format(s.name)
            for par in pars:
                string += ",,,,,,,"
            fout.write(string+"\n")
            continue
        string = "{0}".format(s.name)
        for par in pars:
            keys = ['most_probable',
                    'lower_limit_1sigma', 'upper_limit_1sigma',
                    'lower_limit_2sigma', 'upper_limit_2sigma']
            try:
                for key in keys:
                      string += ",{0:.3f}".format(getattr(s, 'iso'+par)[key])
            except:
                string += ",,,,,"
            try:
                string += ",{0:.3f},{1:.3f}".\
                          format(getattr(s, 'iso'+par)['mean'],\
                                 getattr(s, 'iso'+par)['std'])
            except:
                string += ",,"
        fout.write(string+"\n")
    fout.close()

    print ''
    print '------------------------------------------------------'
    end_time = datetime.datetime.now()
    print '- Date and time: '+end_time.strftime('%d-%b-%Y, %H:%M:%S')
    delta_t = (end_time - start_time).seconds
    hours, remainder = divmod(delta_t, 3600)
    minutes, seconds = divmod(remainder, 60)
    print '- Time elapsed: %sH %sM %sS' % (hours, minutes, seconds)
    print 'Done!'
    print '------------------------------------------------------'
    print ''

def get_isochrone_points(Star, feh_offset, db, nsigma, key_parameter_known):
    '''Looks in the db database for isochrone points within nsigma from
    the mean parameters of the Star and returns those values in a dict.
    '''
    conn = sqlite3.connect(os.path.join(ISOCHRONES_PATH, db))
    conn.row_factory = sqlite3.Row
    logtm = np.log10(Star.teff-nsigma*Star.err_teff)
    logtp = np.log10(Star.teff+nsigma*Star.err_teff)
    c = conn.cursor()
    if key_parameter_known != 'logg' and key_parameter_known != 'plx':
        logger.warning(key_parameter_known+\
                       ' is no a valid key parameter (use logg or plx).')
        return None
    if key_parameter_known == 'logg':
        x = c.execute('SELECT feh, age, mass, logt, logl, logg, mv ' +\
                      'FROM  fa, mtlgv ON fa.fa = mtlgv.fa WHERE '   +\
                      'logt >= ? AND logt <= ? AND '   +\
                      'feh  >= ? AND feh  <= ? AND '   +\
                      'logg >= ? AND logg <= ? ',
                      (logtm, logtp,
                       (Star.feh+feh_offset)-nsigma*Star.err_feh,
                       (Star.feh+feh_offset)+nsigma*Star.err_feh,
                       Star.logg-nsigma*Star.err_logg,
                       Star.logg+nsigma*Star.err_logg)
                     )
    if key_parameter_known == 'plx':
        try:
            Star.M_V = Star.v - 5 * np.log10(1000/Star.plx) + 5.
            Star.err_M_V = np.sqrt(Star.err_v**2 +\
              (np.log10(np.exp(1))**2)*25*(Star.err_plx/Star.plx)**2)
            logger.info('Absolute magnitude and error attributes '+\
                        'added to star object')
        except:
            logger.warning('Could not calculate absolute magnitude. '+\
                           'Star must have v and err_v attributes (vmag).')
            return None
        x = c.execute('SELECT feh, age, mass, logt, logl, logg, mv ' +\
                      'FROM  fa, mtlgv ON fa.fa = mtlgv.fa WHERE '   +\
                      'logt >= ? AND logt <= ? AND '   +\
                      'feh  >= ? AND feh  <= ? AND '   +\
                      'mv >= ? AND mv <= ? ',
                      (logtm, logtp,
                       (Star.feh+feh_offset)-nsigma*Star.err_feh,
                       (Star.feh+feh_offset)+nsigma*Star.err_feh,
                       Star.M_V-nsigma*Star.err_M_V,
                       Star.M_V+nsigma*Star.err_M_V)
                     )
    feh, age = [], []
    mass, logt, logl, logg, mv = [], [], [], [], []
    for xx in x.fetchall():
        feh.append(xx['feh'])
        age.append(xx['age'])
        mass.append(xx['mass'])
        logt.append(xx['logt'])
        logl.append(xx['logl'])
        logg.append(xx['logg'])
        mv.append(xx['mv'])
    conn.close()
    return {'feh' : np.array(feh)-feh_offset,
            'age' : np.array(age),
            'mass': np.array(mass),
            'logt': np.array(logt),
            'logl': np.array(logl),
            'logg': np.array(logg),
            'mv'  : np.array(mv)
            }

def smooth(x, window_len=11, window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[(window_len/2):-(window_len/2)]

def get_isochrone(age, feh, db):
    conn = sqlite3.connect(os.path.join(ISOCHRONES_PATH, db))
    conn.row_factory = sqlite3.Row
    c = conn.cursor()
    x = c.execute('SELECT mass, logt, logl, logg, mv '         +\
                  'FROM  fa, mtlgv ON fa.fa = mtlgv.fa WHERE ' +\
                  'age == ? AND feh == ?', (age, feh))
    mass, logt, logl, logg, mv = [], [], [], [], []
    for xx in x.fetchall():
        mass.append(xx['mass'])
        logt.append(xx['logt'])
        logl.append(xx['logl'])
        logg.append(xx['logg'])
        mv.append(xx['mv'])
    conn.close()
    if not len(mass):
        logger.warning('no {2} isochrone found for age={0} Gyr and [Fe/H]={1}'.
                       format(age, feh, db))
        return None
    else:
        return {'mass': np.array(mass), 'logt': np.array(logt),
                'logl': np.array(logl), 'logg': np.array(logg),
                'mv': np.array(mv)}
