# ztf_fp.py
Library to streamline requesting Zwicky Transient Facility (ZTF) forced photometry.

## A simple command line example

```
python ztf_fp.py 256.7975042 58.0974194 -source_name 2021kjb
```
### Output
```
Sending ZTF request for (R.A.,Decl)=(256.7975042,58.0974194)
wget --http-user=ztffps --http-passwd=dontgocrazy! -O ztffp_493YUZW74Q.txt "https://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi?ra=256.797504&dec=58.097419&jdstart=2459269.628723&jdend=2459329.628723&email=<REDACTED>&userpass=<REDACTED>"

Waiting for the email (rechecking every 20 seconds).
Downloading file...
	wget --http-user=ztffps --http-password=dontgocrazy! -O forcedphotometry_req00015274_lc.txt "https://ztfweb.ipac.caltech.edu/ztf/ops/forcedphot/lc/0/15/req15274/cksum16bf0dfcd7fc6781336f6d0792a8874d/forcedphotometry_req00015274_lc.txt"
Downloading file...
	wget --http-user=ztffps --http-password=dontgocrazy! -O forcedphotometry_req00015274_log.txt "https://ztfweb.ipac.caltech.edu/ztf/ops/forcedphot/lc/0/15/req15274/cksum16bf0dfcd7fc6781336f6d0792a8874d/forcedphotometry_req00015274_log.txt"
Downloaded: 2021kjb_lc.txt
Downloaded: 2021kjb_log.txt
ZTF_r: 15 detections and 3 upper limits.
ZTF_g: 13 detections and 4 upper limits.
Creating ./2021kjb
     ZTF wget log: ./2021kjb/ztffp_493YUZW74Q.txt
     ZTF downloaded file: ./2021kjb/2021kjb_lc.txt
     ZTF downloaded file: ./2021kjb/2021kjb_log.txt
     ZTF figure: ./2021kjb/2021kjb_lc.png
```

## Requirements
### Software installed on your system
- python 3 (tested on python 3.7-3.9)
- wget (MacOS often is missing this. If so, check out brew, macports or anaconda.)
- Email address reachable with IMAP

### Python libraries
- astropy (for possible time and coordinate conversions)
- matplotlib (for possible plotting, but easily removable)
- numpy
- pandas 


## Installation
```
git clone https://github.com/mcstroh/ztf_fp.git
```

## Setup
The library requires your ZTF forced photometry credentials (email and password) and your credentials to access your email address where the links to the ZTF light curve are sent.

For bash/zsh equivalents (e.g., .bashrc, .bash_profile, .zshrc, etc)

```
### ZTF forced photometry service
ztf_email_address="john_doe@mydomain.com"
ztf_email_password="1234567890password"
ztf_email_imapserver="imap.mydomain.com"
ztf_user_password="01234ztf"
export ztf_email_address
export ztf_email_password
export ztf_email_imapserver
export ztf_user_password
```

or for csh/tcsh equivalents (e.g.,.cshrc, tcshrc, etc.)

```
setenv ztf_email_address john_doe@mydomain.com
setenv ztf_email_password 1234567890password
setenv ztf_email_imapserver imap.mydomain.com
setenv ztf_user_password 01234ztf
```

Note that if the email you registered with the ZTF forced photometry is an alias for another address, you will need to also define

```
ztf_user_address="john_doe_1234@mydomain.com"
export ztf_user_address
```
or the equivalent for csh and tcsh shells.

### Check your email configuration - IMPORTANT TO DO THIS FIRST!
It is important to verify that the script is able to use your email credentials so that you don't spam the ZTF forced photometry service with duplicate requests while you debug your email setup.

Run the following to see if the script can access your email with the provided credentials.

#### Input
```
python ztf_fp.py -emailtest
```
#### Successful output
```
Your email inbox was found and contains 600 messages.
If this is not correct, please check your settings.
```
#### Output if it definitely is not working
```
Your inbox was not located. Please check your settings.
```


## More about running on the command line

Sexagesimal coordinates are also supported:
```
python ztf_fp.py 17:07:11.40 +58:05:50.71 -source_name 2021kjb
```

By default, the script requests forced photometry from the 60 days prior to the moment you submit the job (equivalent to ```-days 60```), but this is easily changed. 

You may also specify date ranges in JD (```-jdstart``` and ```jdend``` arguments) or MJD (```-mjdstart``` and ```mjdend``` arguments).


### List of available options
```
usage: ztf_fp.py <ra> <decl> [-h] [-logfile [logfile]] [-plotfile [plotfile]] [-emailtest] [-source_name [source_name]] [-mjdstart [mjdstart]]
                             [-mjdend [mjdend]] [-jdstart [jdstart]] [-jdend [jdend]] [-ztf_all_jd] [-days [days]] [-emailcheck [emailcheck]]
                             [-skip_clean] [-directory_path [directory_path]] [-fivemindelay [fivemindelay]] [-skip_plot]
                             [ra] [decl]

Grab ZTF forced photometry on the given location.

positional arguments:
  ra                    Right ascension of the target. Can be provided in DDD.ddd or HH:MM:SS.ss formats.
  decl                  Declination of the target. Can be provided in +/-DDD.ddd or DD:MM:SS.ss formats.

optional arguments:
  -h, --help            show this help message and exit
  -logfile [logfile]    Log file to process instead of submitting a new job.
  -plotfile [plotfile]  Light curve file to plot instead of submitting and downloading a new job.
  -emailtest            Test your email settings. This is performed without sending a request to the ZTF server.
  -source_name [source_name]
                        Source name that will be used to name output files.
  -mjdstart [mjdstart]  Start of date range for forced photometry query. Overrides -jdstart.
  -mjdend [mjdend]      End of date range for forced photometry query. Overrides -jdstop
  -jdstart [jdstart]    Start of date range for forced photometry query.
  -jdend [jdend]        End of date range for forced photometry query.
  -ztf_all_jd           Use the full range of ZTF public dates.
  -days [days]          Number of days prior to jdend to query (or number of days prior to today if jdend is not given).
  -emailcheck [emailcheck]
                        How often to recheck your email for the ZTF results.
  -skip_clean           After completion skip placing all output files in the same directory.
  -directory_path [directory_path]
                        Path to directory for clean-up. Requires -directory option.
  -fivemindelay [fivemindelay]
                        How often (in seconds) to query the email after 5 minutes have elapsed.
  -skip_plot            Skip making the plot. Useful for automated and batch use cases, or if user wants to use their personal plotting code.
```

## Automated / batch job requests

This file supports being imported as a library and called for batch processing. Below is an example of an *external script* using the Pool module to call the function over a number of positions.

```
from multiprocessing import Pool, cpu_count
import ztf_fp

def ztf_forced_photometry(row):

    # Send ZTF forced photometry request
    ztf_file_names = ztf_fp.run_ztf_fp(days=20, ra=row.ra, decl=row.decl, source_name=row.name, directory_path='/tmp', verbose=True)
    
    #
    # Do something intelligent with the output
    #


def batch_ztf_forced_photometry():

    #
    # Populate a pandas dataframe with the positions you are interested in
    # This is a simple example using a couple of SN Ia discovered by the Young Supernova Experiment
    data = [['2021kcc',194.0331000,-4.9606139],['2021jze',154.9723750,-3.7354000]]
    
    # The column names are used in the ztf_forced_photometry function defined above
    df = pd.DataFrame(data, columns = ['name','ra','decl'])


    n_workers = int(np.floor(cpu_count()/2)) # Requires casting as an integer for pool argument / don't hit the ZTF server too hard
    pool_vals = [x for _, x in df.iterrows()] # Turn to list for pool mapping function below
    print(pool_vals)
    print(f"Processing {len(pool_vals)} ZTF forced photometry requests utilizing {n_workers} workers.")
    with Pool(n_workers) as ztf_pool:
         res = ztf_pool.map(ztf_forced_photometry, pool_vals)  
```

If you are not interested in the output files, you can direct the files to /tmp using the ```-directory_path``` option.
