import numpy as np
import os
import logging
<<<<<<< HEAD
import tempfile
from config import *
=======
from .config import *
>>>>>>> astroChasqui/master

logger = logging.getLogger(__name__)

temp_dir = tempfile.mkdtemp(dir='./', prefix='.tmp_moog_')
temp_dir += '/'
print 'Remember to clean up hidden directory '+temp_dir+' if exit abnormally'


class Driver:
    """Set the options for your MOOG driver."""
    def __init__(self, mode='abfind'):
        self.mode = mode
        self.standard_out = 'moog.std'
        self.summary_out = 'moog.sum'
        self.model_in = 'model.in'
        self.lines_in = 'lines.in'
        self.plot = 0
        self.hfs_species = None
        # batch.par, MOOG and temporary files are stored in a random temprary directory
        # so that multipple instance on q2 can run in the same directory

    def create_file(self, file_name="batch.par", temporary_dir=''):
        """Creates the MOOG driver file."""
        self.file_name = file_name
        f = open(temporary_dir+file_name, 'w')
        if self.mode == 'abfind':
            if self.hfs_species:
                f.write('blends\n')
            else:
                f.write('abfind\n')
        else:
            f.write(self.mode+'\n')
        f.write('standard_out "'+self.standard_out+'"\n')
        f.write('summary_out  "'+self.summary_out+'"\n')
        f.write('model_in     "'+self.model_in+'"\n')
        f.write('lines_in     "'+self.lines_in+'"\n')
        f.write('atmosphere   1\n')
        f.write('molecules    1\n')
        f.write('lines        1\n')
        f.write('flux/int     0\n')
        f.write('damping      1\n')
        f.write('freeform     1\n')
        f.write('plot         '+str(self.plot)+'\n')
        if self.hfs_species:
            f.write('blenlimits\n')
            f.write(' 2.0 0.01 '+self.hfs_species+'\n')
        if self.mode == 'cog':
            f.write('coglimits\n')
            f.write('  -6.5 -3.5 0.1 0 0\n')
        f.close()


def create_model_in(Star, file_name='model.in', temporary_dir=''):
    """Creates a model atmosphere file for MOOG from the model_atmosphere
    attribute of a Star object.
    """
    try:
        Star.vt
    except:
        logger.error('Moog model_in file requires a microturbulence (vt)')
        return None
    if hasattr(Star, 'model_atmosphere'):
        Star.moog_model_in_name = temporary_dir + file_name
    else:
        logger.error('No model data to write to moog model_in file.')
        return None

    if hasattr(Star, 'feh_model'):
        feh = Star.feh_model
    else:
        feh = Star.feh

    with open(temporary_dir+'head.tmp', 'w') as f:
        f.write('KURUCZ\n')
        f.write('TEFF='+str(Star.teff)+',LOGG='+str(Star.logg)+
                ',[FE/H]='+str(feh)+','+Star.model_atmosphere_grid+'\n')
        nd = len(Star.model_atmosphere['T'])
        f.write('ND=       '+str(nd)+'\n')

    with open(temporary_dir+'body.tmp', 'w') as f:
        for idx in range(nd):
            f.write("{0:.8E} {1:.1F} {2:.3E} {3:.3E} {4:.3E}\n".format(\
                    Star.model_atmosphere['RHOX'][idx],\
                    Star.model_atmosphere['T'][idx],\
                    Star.model_atmosphere['P'][idx],\
                    Star.model_atmosphere['XNE'][idx],\
                    Star.model_atmosphere['ABROSS'][idx])
                   )

    with open(temporary_dir+'tail.tmp', 'w') as f:
        f.write('%5.2F\n' %Star.vt)
        if Star.model_atmosphere_grid != 'marcs':
            path = os.path.join(MODATM_PATH, 'kurucz')
            fabund = open(os.path.join(path, 'p00.'+Star.model_atmosphere_grid),
                          'r')
        else:
            path = os.path.join(MODATM_PATH, 'marcs')
            fabund = open(os.path.join(path, 'z+0.00'), 'r')

        line = fabund.readline()
        f.write(line[0:12]+' '+str(feh)+'\n')
        line = fabund.readline()
        while line:
            species = line[0:2]
            if Star.model_atmosphere_grid == 'marcs':
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
    with open(temporary_dir+file_name, 'w') as outfile:
        for one_file in file_list:
            with open(one_file) as infile:
                outfile.write(infile.read())
    for one_file in file_list:
        os.unlink(one_file)
    logger.info('Moog infile model atmosphere created: '+file_name)


def create_lines_in(Star, species=0, file_name='lines.in', add_error=False, temporary_dir=''):
    """Creates a line list file for MOOG"""
    """LM: Added flag to create a linelist with EWs increased by their errors"""
    if species > 0:
        idx = np.where(np.logical_and(Star.linelist['species'] == species,\
                                       Star.linelist['ew'] >= 0.1))[0] # LM
    else:
        #species = 0 means all species
        idx = np.where(np.logical_and(Star.linelist['species'] > species,\
                                       Star.linelist['ew'] >= 0.1))[0] # LM

    nlines = len(idx)
    if nlines == 0:
        logger.warning('No lines found for '+Star.name)
        return False
    else:
        logger.info(str(nlines)+' lines found for '+Star.name)
    gf_values = Star.linelist['gf'][idx]
    gf10 = [10**gfx for gfx in Star.linelist['gf'][idx] if gfx >= 0]
    if len(gf10) == len(Star.linelist['gf'][idx]):
        logger.info('all gf values for this species are positive --> 10^gf')
        #gf_values = gf10
        Star.linelist['gf'][idx] = gf10
    #Star.linelist['gf'][idx] = gf_values

    with open(temporary_dir+file_name, 'w') as f:
        f.write("MOOG linelist created by q2\n")
        if add_error:
            print 'Adding errors to EWS - be careful! '
            """LM: EWs are increased by their corresponding errors. Error on abundance due
            to the error in EW will be calculated by subtracting the unperturbed
            abundance from this value """
            for lidx in idx:
                f.write("{0:10.4f} {1:4.1f} {2:6.3f} {3:5.3f} 0 0 {4:5.1f}\n".format(\
                    Star.linelist['wavelength'][lidx],\
                    Star.linelist['species'][lidx],\
                    Star.linelist['ep'][lidx],\
                    Star.linelist['gf'][lidx],\
                    Star.linelist['ew'][lidx]+Star.linelist['ew_r'][lidx])
                   )
                # LM switched from '3 0' to '0 0', ie use the Unsold approximation
                # when Barklem values are not avaialble
        else:
            for lidx in idx:
                f.write("{0:10.4f} {1:4.1f} {2:6.3f} {3:5.3f} 0 0 {4:5.1f}\n".format(\
                    Star.linelist['wavelength'][lidx],\
                    Star.linelist['species'][lidx],\
                    Star.linelist['ep'][lidx],\
                    Star.linelist['gf'][lidx],\
                    Star.linelist['ew'][lidx])
                   )

    Star.linelist['gf'][idx] = gf_values

    logger.info('Moog line list created: '+file_name)
    return True


def abfind(Star, species, species_id):
    """Runs MOOG with abfind driver for a given Star and species

    Star is a star object; must have all attributes in place
    species could be 26.0 for Fe I, for example
    species_id is a string that will become a new attribute for the Star object
    Example: abfind(s, 26.1, 'fe2')
    s.fe2 #shows result from abfind
    MD is the moog driver object
    """
    """LM: added the computation of perturbed abundances"""
    if Star.use_errors:
        ab_pert = abfind_perturbed(Star, species)

    k = Star.linelist['species'] == species
    negs = [wx for wx in Star.linelist['wavelength'][k] if wx < 0]
    if len(negs) == 0:
        MD = Driver() #normal
    else:
        MD = Driver() #hfs
        MD.hfs_species = str(round(species))
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)
    #temporary files are saved into a temporary directoru
    #MD.standard_out = os.path.join('.q2', 'moog.std')
    #MD.summary_out = os.path.join('.q2', 'moog.sum')
    #MD.model_in = os.path.join('.q2', 'model.in')
    #MD.lines_in = os.path.join('.q2', 'lines.in')
    MD.create_file('batch.par',temporary_dir=temp_dir)

    create_model_in(Star, file_name=MD.model_in, temporary_dir=temp_dir)
    found_lines = create_lines_in(Star, species=species, file_name=MD.lines_in, temporary_dir=temp_dir)
    if not found_lines:
        logger.warning('Did not run abfind (no lines found)')
        return False
    logfile ='moog.log'
    os.chdir(temp_dir)
    os.system('MOOGSILENT > '+logfile+' 2>&1 ')
    os.chdir('../')
    f = open(temp_dir + MD.summary_out, 'r')
    line, stop = '', False
    while line[0:10] != 'wavelength':
        line = f.readline()
    if 'ID' in line:
        moogjul2014 = True
    else:
        moogjul2014 = False
    while not stop: #looping required for multiple iterations (molecules)
        ww, ep, ew, rew, ab, difab = [], [], [], [], [], []
        while line:
            line = f.readline()
            if line[0:7] == 'average': break
            linesplit = line.split()
            if float(linesplit[6]) > 999.: #exclude dummies (hfs)
                continue
            ww.append(float(linesplit[0]))
            if moogjul2014: #MOOGJUL2014 adds a new column 'ID' to moog.sum
                ep.append(float(linesplit[2]))
                ew.append(float(linesplit[4]))
                rew.append(float(linesplit[5]))
                ab.append(float(linesplit[6]))
            else: #older versions of MOOG don't have 'ID' but 'EP' in 2nd col
                ep.append(float(linesplit[1]))
                ew.append(float(linesplit[3]))
                rew.append(float(linesplit[4]))
                ab.append(float(linesplit[5]))
            difab.append(None)
        while line: #to break out of multiple iterations loop if done
            line = f.readline()
            if line[0:10] == 'wavelength':
                stop = False
                break
            stop = True
    f.close()
    os.unlink(temp_dir+MD.file_name)
    os.unlink(temp_dir+MD.model_in)
    os.unlink(temp_dir+MD.lines_in)
    os.unlink(temp_dir+MD.summary_out)
    os.unlink(temp_dir+MD.standard_out)
    os.unlink(temp_dir+logfile)
    if os.path.isfile(temp_dir+'fort.99'):
        os.unlink(temp_dir+'fort.99')

    if not Star.use_errors:
        ab_pert = np.array(ab)+0.01

    x = {'ww': np.array(ww), 'ep': np.array(ep), 'ew': np.array(ew),\
        'rew': np.array(rew), 'ab': np.array(ab), 'difab': np.array(difab),
        'ab_e': np.abs(ab_pert-np.array(ab))}
    setattr(Star, species_id, x)
    logger.info('Successfully ran abfind')
    return True

def abfind_perturbed(Star, species):
    """LM
    Added to perform the computation of perturbed EWs

    :return:
    """
    k = Star.linelist['species'] == species
    negs = [wx for wx in Star.linelist['wavelength'][k] if wx < 0]
    if len(negs) == 0:
        MD = Driver() #normal
    else:
        MD = Driver() #hfs
        MD.hfs_species = str(round(species))
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)

    MD.create_file('batch.par',temporary_dir=temp_dir)

    create_model_in(Star, file_name=MD.model_in, temporary_dir=temp_dir)
    found_lines = create_lines_in(Star, species=species, file_name=MD.lines_in, add_error=True, temporary_dir=temp_dir)
    if not found_lines:
        logger.warning('Did not run abfind (no lines found)')
        return False
    logfile ='moog.log'
    os.chdir(temp_dir)
    os.system('MOOGSILENT > '+logfile+' 2>&1 ')
    os.chdir('../')
    f = open(temp_dir + MD.summary_out, 'r')
    line, stop = '', False
    while line[0:10] != 'wavelength':
        line = f.readline()
    if 'ID' in line:
        moogjul2014 = True
    else:
        moogjul2014 = False
    while not stop: #looping required for multiple iterations (molecules)
        ww, ep, ew, rew, ab, difab = [], [], [], [], [], []
        while line:
            line = f.readline()
            if line[0:7] == 'average': break
            linesplit = line.split()
            if float(linesplit[6]) > 999.: #exclude dummies (hfs)
                continue
            ww.append(float(linesplit[0]))
            if moogjul2014: #MOOGJUL2014 adds a new column 'ID' to moog.sum
                ep.append(float(linesplit[2]))
                ew.append(float(linesplit[4]))
                rew.append(float(linesplit[5]))
                ab.append(float(linesplit[6]))
            else: #older versions of MOOG don't have 'ID' but 'EP' in 2nd col
                ep.append(float(linesplit[1]))
                ew.append(float(linesplit[3]))
                rew.append(float(linesplit[4]))
                ab.append(float(linesplit[5]))
            difab.append(None)
        while line: #to break out of multiple iterations loop if done
            line = f.readline()
            if line[0:10] == 'wavelength':
                stop = False
                break
            stop = True
    f.close()
    os.unlink(temp_dir+MD.file_name)
    os.unlink(temp_dir+MD.model_in)
    os.unlink(temp_dir+MD.lines_in)
    os.unlink(temp_dir+MD.summary_out)
    os.unlink(temp_dir+MD.standard_out)
    os.unlink(temp_dir+logfile)
    if os.path.isfile(temp_dir+'fort.99'):
        os.unlink(temp_dir+'fort.99')

    return np.array(ab)



def cog(Star, species, cog_id):
    """Runs MOOG with cog driver for a given Star and species

    Star is a star object; must have all attributes need by MOOG set.
    species could be 26.0 for Fe I, for example. cog_id is a string that
    will become a new attribute for the Star object. For example:
    >>>cog(s, 26.1, 'cog_fe2')
    s.cog_fe2 #shows result from cog
    MD is the moog driver object
    """
    k = Star.linelist['species'] == species
    #negs = [wx for wx in Star.linelist['wavelength'][k] if wx < 0]
    MD = Driver(mode='cog')
    MD.create_file(temporary_dir=temp_dir)
    create_model_in(Star,temporary_dir=temp_dir)
    found_lines = create_lines_in(Star, species=species, temporary_dir=temp_dir)
    if not found_lines:
        logger.warning('Did not run cog (no lines found)')
        return False
    os.chdir(temp_dir)
    os.system('MOOGSILENT > '+logfile+' 2>&1 ')
    os.chdir('../')

    f = open(temp_dir+MD.summary_out, 'r')
    line = f.readline()
    cog_obj = {}
    while line:
        line = f.readline()
        if line.startswith('wavelength'):
            npt = int(line.split('=')[5]) #number of cog points
            #wavelength = round(float(line.split('=')[1].split()[0]), 1)
            wavelength = float(line.split('=')[1].split()[0])
            line = f.readline()
            x, y = [], []
            for i in range(int(np.ceil(npt/5.))):
                line = f.readline()
                for j in range(len(line.split())/2):
                    x.append(float(line.split()[2*j].replace(',', '')))
                    y.append(float(line.split()[2*j+1]))
            cog_obj[wavelength] = {'loggf': np.array(x), 'logrw': np.array(y)}
    f.close()

    os.unlink(temp_dir+MD.file_name)
    os.unlink(temp_dir+MD.model_in)
    os.unlink(temp_dir+MD.lines_in)
    os.unlink(temp_dir+MD.summary_out)
    os.unlink(temp_dir+MD.standard_out)
    os.unlink(temp_dir+'moog.log')

    setattr(Star, cog_id, cog_obj)


def delete_tempdir():
    os.system('rm -rf ' + temp_dir)
