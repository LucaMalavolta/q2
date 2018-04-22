import os
import tempfile
import sys
sys.path.insert(0, '/Users/malavolta/Astro/CODE')
import q2

"""
def create_file():
    os.system('touch AAAA')

temp_dir = tempfile.mkdtemp(dir='./')
print temp_dir
os.chdir(temp_dir)
create_file()
os.chdir('..')
os.system('rm -r '+temp_dir)
"""

data = q2.Data('standards_stars.csv', 'standards_lines_with_errors.csv')
sp = q2.specpars.SolvePars()
sp.errors = True
q2.specpars.solve_all(data, sp, 'solution.csv', 'Sun', use_variant=True)
q2.moog.delete_tempdir()
