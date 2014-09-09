import moog
from star import Star
import numpy as np
import datetime
import logging
from astropy.io import ascii
from scipy import interpolate
import os
from config import *


logger = logging.getLogger(__name__)


def all(Data, species_ids, output_file, reference=None, grid='odfnew'):
    print('------------------------------------------------------')
    print('Initializing ...')
    start_time = datetime.datetime.now()
    print('- Date and time: '+str(start_time))
    print('- Model atmospheres: '+grid)
    print('- Star data: '+Data.star_data_fname)
    print('- Line list: '+Data.lines_fname)
    if reference:
        print('- Reference star: '+reference)
    print('------------------------------------------------------')
    if reference:
        ref = Star(reference)
        ref.get_data_from(Data)
        ref.get_model_atmosphere(grid)
    else:
        ref = None
    fout = open(output_file, 'wb')
    header = 'id'
    for species_id in species_ids:
        header += ','+species_id+',err_'+species_id
        if reference:
            header += ',['+species_id+'],err_['+species_id+']'
    fout.write(header+'\n')
    for star_id in Data.star_data['id']:
        line = star_id
        print('')
        print('*'*len(star_id))
        print(star_id)
        print('*'*len(star_id))
        s = Star(star_id)
        try:
            s.get_data_from(Data)
            s.get_model_atmosphere(grid)
        except:
            print('No data available')
            logger.warning('Could not get all the necessary data')
            line += ','*(len(species_ids)*2)
            if reference:
                line += ','*(len(species_ids)*2)
            fout.write(line+'\n')
            continue
        one(s, species_ids, ref)
        for species_id in species_ids:
            print('\n'+species_id+'\n'+'-'*len(species_id))
            if not hasattr(s, species_id):
                print('No data available')
                logger.warning('There are no '+species_id+' abundances '+\
                               'for this star')
                line += ',,'
                if reference:
                    line += ',,'
                continue
            mab = np.mean(getattr(s, species_id)['ab'])
            sab = np.std(getattr(s, species_id)['ab'])
            print("ABS = {0:6.3f} +/- {1:6.3f}".format(mab, sab))
            line += ',{0:.3f},{1:.3f}'.format(mab, sab)
            if reference:
                da = getattr(s, species_id)['difab']
                mda = np.ma.masked_array(da, np.isnan(da))
                mdifab = np.mean(mda)
                sdifab = np.std(mda)
                print("DIF = {0:6.3f} +/- {1:6.3f}".format(mdifab, sdifab))
                line += ',{0:.3f},{1:.3f}'.format(mdifab, sdifab)
            else:
                mdifab = 0
            print('')
            llhd1 = 'Wavelength   ABS    RES '
            llhd2 = '----------  ----- ------'
            if reference:
                llhd1 += '   DIF    RES '
                llhd2 += '  -----  -----'
            print(llhd1+'\n'+llhd2)
            for wi, ab, difab in \
                zip(getattr(s, species_id)['ww'],
                    getattr(s, species_id)['ab'],
                    getattr(s, species_id)['difab']):
                if reference:
                    print("{0:10.4f} {1:6.3f} {2:6.3f} {3:6.3f} {4:6.3f}".
                           format(wi, ab, ab-mab, difab, difab-mdifab))
                else:
                    print("{0:10.4f} {1:6.3f} {2:6.3f}".
                           format(wi, ab, ab-mab))

        fout.write(line+'\n')
    fout.close()
    print('')
    print('------------------------------------------------------')
    end_time = datetime.datetime.now()
    print('- Date and time: '+str(end_time))
    print('- Time elapsed: '+str(end_time - start_time))
    print('Done!')
    print('------------------------------------------------------')
    print('')


def one(Star, species_ids, Ref=object, silent=True):
    logger.info('Working on: '+Star.name)
    for species_id in species_ids:
        species = getsp(species_id)
        if species == None:
            logger.warning('Not doing calculations for: '+species_id)
            continue
        logger.info('Working on: '+species_id)
        moog.abfind(Star, species, species_id)
        if not hasattr(Star, species_id):
            logger.warning('Did not calculate '+species_id+' abundances')
            continue

        if species_id == 'OI':
            if not silent:
                print('')
                print('777 nm oxygen abundances will be NLTE corrected')
            ao = []
            for wx in [7771.94, 7774.16, 7775.39]:
                k = np.where(abs(Star.OI['ww']-wx) < 0.05)
                if len(k[0]) == 1:
                    ao.append(np.mean(Star.OI['ab'][k]))
                else:
                    ao.append(0)
            aon = nlte_triplet(Star.teff, Star.logg, Star.feh, ao,
                               silent=silent)
            k= np.where(np.array(ao) > 0)
            getattr(Star, species_id)['ab'] = aon[k]

        if hasattr(Ref, 'name'):
            logger.info('Diferential analysis: '+Ref.name)
            if Star.name == Ref.name:
                logger.warning('Reference star object redefined!')
                Ref = Star
            if not hasattr(Ref, species_id):
                logger.info('Calculating reference star abundances: '+Ref.name)
                moog.abfind(Ref, species, species_id)

                if species_id == 'OI':
                    if not silent:
                        print('')
                        print('777 nm oxygen abundances will be NLTE '\
                              +'corrected (Reference)')
                    ao = []
                    for wx in [7771.94, 7774.16, 7775.39]:
                        k = np.where(abs(Ref.OI['ww']-wx) < 0.05)
                        if len(k[0]) == 1:
                            ao.append(np.mean(Ref.OI['ab'][k]))
                        else:
                            ao.append(0)
                    aon = nlte_triplet(Ref.teff, Ref.logg, Ref.feh, ao)
                    k= np.where(np.array(ao) > 0)
                    getattr(Ref, species_id)['ab'] = aon[k]

            else:
                logger.info('Reference star has '+species_id+\
                            ' abundances computed already: '+Ref.name)

            ws = getattr(Star, species_id)['ww']
            wr = getattr(Ref, species_id)['ww']
            ww = np.intersect1d(ws, wr)
            k  = [i for i, w in zip(range(len(ws)), ws) if w in ww]
            kr = [i for i, w in zip(range(len(wr)), wr) if w in ww]
            a = getattr(Star, species_id)['ab'][k] - \
                getattr(Ref, species_id)['ab'][kr]
            ax, ix = [], 0
            for wx in ws:
                if wx in ww:
                    ax.append(a[ix])
                    ix += 1
                else:
                    ax.append(None)
            getattr(Star, species_id)['difab'] = ax
        if not silent and len(species_ids) > 1:
            print(species_id + ' done')
    if not silent and len(species_ids) > 1:
        print('All species completed')


def getsp(species_id):
    d = {
         'CI'  :  6.0,
         'OI'  :  8.0,
         'NaI' : 11.0,
         'MgI' : 12.0,
         'AlI' : 13.0,
         'SiI' : 14.0,
         'SI'  : 16.0,
         'CaI' : 20.0,
         'ScI' : 21.0,
         'ScII': 21.1,
         'TiI' : 22.0,
         'TiII': 22.1,
         'VI'  : 23.0,
         'CrI' : 24.0,
         'CrII': 24.1,
         'MnI' : 25.0,
         'FeI' : 26.0,
         'FeII': 26.1,
         'CoI' : 27.0,
         'NiI' : 28.0,
         'CuI' : 29.0,
         'ZnI' : 30.0,
         'SrI' : 38.0,
         'YII' : 39.1,
         'ZrII': 40.1,
         'BaII': 56.1
         }
    try:
        species = d[species_id]
    except:
        logger.warning('species id not recognized: '+species_id)
        return None
    return species


def nlte_triplet(teff, logg, feh, ao, silent=True):
    grid = ascii.read(os.path.join(OTHER_PATH ,'nlte_triplet.csv'))

    t,g,f,dao0,dao1,dao2=[],[],[],[],[],[]
    for i in range(640):
        rg = range(i*7, i*7+7)
        x0 = interpolate.griddata(grid['ao'][rg], grid['dao0'][rg],\
                                  ao[0], method='cubic')
        x1 = interpolate.griddata(grid['ao'][rg], grid['dao1'][rg],\
                                  ao[1], method='cubic')
        x2 = interpolate.griddata(grid['ao'][rg], grid['dao2'][rg],\
                                  ao[2], method='cubic')
        x0, x1, x2 = float(x0), float(x1), float(x2)
        t.append(grid['teff'][rg[0]])
        g.append(grid['logg'][rg[0]])
        f.append(grid['feh'][rg[0]])
        dao0.append(x0)
        dao1.append(x1)
        dao2.append(x2)
    t = np.array(t)
    g = np.array(g)
    f = np.array(f)
    dao0 = np.array(dao0)
    dao1 = np.array(dao1)
    dao2 = np.array(dao2)

    tt,ff,dao00,dao11,dao22=[],[],[],[],[]
    for i in range(160):
        rg =range(i*4, i*4+4)
        x0 = interpolate.griddata(g[rg], dao0[rg], logg, method='cubic')
        x1 = interpolate.griddata(g[rg], dao1[rg], logg, method='cubic')
        x2 = interpolate.griddata(g[rg], dao2[rg], logg, method='cubic')
        x0, x1, x2 = float(x0), float(x1), float(x2)
        tt.append(t[rg[0]])
        ff.append(f[rg[0]])
        dao00.append(x0)
        dao11.append(x1)
        dao22.append(x2)
    tt = np.array(tt)
    ff = np.array(ff)
    dao00 = np.array(dao00)
    dao11 = np.array(dao11)
    dao22 = np.array(dao22)

    t,dao0,dao1,dao2=[],[],[],[]
    for i in range(16):
        rg =range(i*10, i*10+10)
        x0 = interpolate.griddata(ff[rg], dao00[rg], feh, method='cubic')
        x1 = interpolate.griddata(ff[rg], dao11[rg], feh, method='cubic')
        x2 = interpolate.griddata(ff[rg], dao22[rg], feh, method='cubic')
        x0, x1, x2 = float(x0), float(x1), float(x2)
        t.append(tt[rg[0]])
        dao0.append(x0)
        dao1.append(x1)
        dao2.append(x2)
    t = np.array(t)
    dao0 = np.array(dao0)
    dao1 = np.array(dao1)
    dao2 = np.array(dao2)

    x0 = interpolate.griddata(t, dao0, teff, method='cubic')
    x1 = interpolate.griddata(t, dao1, teff, method='cubic')
    x2 = interpolate.griddata(t, dao2, teff, method='cubic')
    x0, x1, x2 = float(x0), float(x1), float(x2)

    x0 = x0 - 0.0355
    x1 = x1 - 0.0180
    x2 = x2 - 0.0000
    
    if not silent:
        print('Wavelength (A) | A(O) LTE | Correction | A(O) NLTE')
        print("   7771.9      |  {0:6.3f}  |    {1:5.3f}   | {2:6.3f}".
              format(ao[0], x0, ao[0]-x0))
        print("   7774.2      |  {0:6.3f}  |    {1:5.3f}   | {2:6.3f}".
              format(ao[1], x1, ao[1]-x1))
        print("   7775.4      |  {0:6.3f}  |    {1:5.3f}   | {2:6.3f}".
              format(ao[2], x2, ao[2]-x2))
    ax = [ao[0]-x0, ao[1]-x1, ao[2]-x2]

    aon = np.ma.masked_array(ax,np.isnan(ax))

    if not silent:
        print("A(O) LTE  = {0:6.3f} +/- {1:5.3f}".
              format(np.mean(ao), np.std(ao)))
        print("A(O) NLTE = {0:6.3f} +/- {1:5.3f}".
              format(np.mean(aon), np.std(aon)))

    return aon
