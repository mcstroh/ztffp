#! python
#
# Obtain ZTF forced photometry on a list of coordinates
#
# M.C. Stroh (Northwestern University)
#

from multiprocessing import Pool, cpu_count
import numpy as np
import pandas as pd

import ztffp

def ztf_forced_photometry(row):

    # Send ZTF forced photometry request
    ztf_file_names = ztffp.run_ztf_fp(days=20,
                                       ra=row.ra,
                                       decl=row.decl,
                                       source_name=row.source_name,
                                       directory_path='/tmp',
                                       verbose=False,
                                       do_plot=False)

    #
    # Do something intelligent with the output
    #     For this example we'll simply list the files
    for ztf_file_name in ztf_file_names:
        print(ztf_file_name)

    return


def batch_ztf_forced_photometry():

    # Populate a pandas dataframe with the positions you are interested in
    # This is a simple example using a couple of SN Ia discovered by the Young Supernova Experiment
    data = [['2021kcc', 194.0331000, -4.9606139],
            ['2021jze', 154.9723750, -3.7354000]]
    
    # The column names are used in the ztf_forced_photometry function defined above
    df = pd.DataFrame(data, columns=['source_name', 'ra', 'decl'])

    # Requires casting as an integer for pool argument / don't hit the ZTF server too hard
    n_workers = int(np.floor(cpu_count()/2))

    # Turn to list for pool mapping function below
    pool_vals = [x for _, x in df.iterrows()]

    # Now send to the list to be completed
    print(f"Processing {len(pool_vals)} ZTF forced photometry requests utilizing {n_workers} workers.")
    ztf_pool = Pool(n_workers)
    ztf_pool.map(ztf_forced_photometry, pool_vals)
    ztf_pool.close()
    ztf_pool.join()


if __name__ == "__main__":
    batch_ztf_forced_photometry()
