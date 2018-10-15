import numpy as np
import logging
import modatm
from config import *
from tools import read_csv

def get_model(teff, logg, vt, feh, grid,  file_name='model.in', temporary_dir=''):
    """Creates a model atmosphere file for MOOG from the model_atmosphere
    attribute of a Star object.
    """
    model_atmosphere = modatm.interpolate(teff, logg, feh, grid)
    model_atmosphere_grid = grid

    with open(temporary_dir+'head.tmp', 'w') as f:
        f.write('KURUCZ\n')
        f.write('TEFF='+str(teff)+',LOGG='+str(logg)+
                ',[FE/H]='+str(feh)+','+model_atmosphere_grid+'\n')
        nd = len(model_atmosphere['T'])
        f.write('ND=       '+str(nd)+'\n')

    with open(temporary_dir+'body.tmp', 'w') as f:
        for idx in range(nd):
            f.write("{0:.8E} {1:.1F} {2:.3E} {3:.3E} {4:.3E}\n".format(\
                    model_atmosphere['RHOX'][idx],\
                    model_atmosphere['T'][idx],\
                    model_atmosphere['P'][idx],\
                    model_atmosphere['XNE'][idx],\
                    model_atmosphere['ABROSS'][idx])
                   )

    with open(temporary_dir+'tail.tmp', 'w') as f:
        f.write('%5.2F\n' %vt)
        if model_atmosphere_grid != 'marcs':
            path = os.path.join(MODATM_PATH, 'kurucz')
            fabund = open(os.path.join(path, 'p00.'+model_atmosphere_grid),
                          'r')
        else:
            path = os.path.join(MODATM_PATH, 'marcs')
            fabund = open(os.path.join(path, 'z+0.00'), 'r')

        line = fabund.readline()
        f.write(line[0:12]+' '+str(feh)+'\n')
        line = fabund.readline()
        while line:
            species = line[0:2]
            if model_atmosphere_grid == 'marcs':
                abund = float(line[3:9])+feh
                #alpha-element enhancement
                if species==' 8' or species=='10' or species=='12' or \
                   species=='14' or species=='16' or species=='18' or \
                   species=='20' or species=='22':
                    afe = -0.4*feh
                    if feh >=  0: afe=0.0
                    if feh <= -1: afe=0.4
                    abund = abund+afe
            else:
                abund = 12.+np.log10(np.power(10, float(line[3:9]))/0.92040)+ \
                        feh
            abund = str('%5.2F' %abund)
            f.write(species+' '+abund+'\n')
            line = fabund.readline()
        fabund.close()
        f.write('NMOL      22\n')
        f.write('  101.0   106.0   107.0   108.0   112.0  126.0\n')
        f.write('  606.0   607.0   608.0\n')
        f.write('  707.0   708.0\n')
        f.write('  808.0   812.0   822.0\n')
        f.write('  10108.0 60808.0\n')
        f.write('  6.1     7.1     8.1   12.1  22.1  26.1\n')

    file_list = [temporary_dir+'head.tmp', temporary_dir+'body.tmp', temporary_dir+'tail.tmp']
    with open(file_name, 'w') as outfile:
        for one_file in file_list:
            with open(one_file) as infile:
                outfile.write(infile.read())
    for one_file in file_list:
        os.unlink(one_file)
    #logger.info('Moog infile model atmosphere created: '+file_name)
