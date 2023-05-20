#! python 
#
# Grab ZTF forced photometry on a field
#
# M.C. Stroh (Northwestern University)
#
# With contributions from: 
#     Wynn Jacobson-Galan
#     David Jones
#     Candice Stauffer
#
#

import argparse
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from datetime import datetime
import email
import imaplib
import numpy as np
import os
import pandas as pd
import random
import re
import shutil
import string
import subprocess
import time
import warnings
# Disable warnings from log10 when there are non-detections
warnings.filterwarnings("ignore")


import matplotlib
matplotlib.use('AGG') # Run faster, comment out for interactive plotting
import matplotlib.pyplot as plt

#
# Generic ZTF webdav login
#
_ztfuser = "ztffps"
_ztfinfo = "dontgocrazy!"


#
# Import ZTF email user information
#
def import_credentials() -> None:
    '''Load ZTF credentials from environmental variables.'''

    try:

        global _ztffp_email_address
        _ztffp_email_address = os.getenv("ztf_email_address", None)
        global _ztffp_email_password
        _ztffp_email_password = os.getenv("ztf_email_password", None)
        global _ztffp_email_server
        _ztffp_email_server = os.getenv("ztf_email_imapserver", None)
        # The email address associated with the ZTF FP server may be an alias,
        # so allow for that possiblility
        global _ztffp_user_address
        if 'ztf_user_address' in os.environ:
            _ztffp_user_address = os.getenv("ztf_user_address", None)
        # Assume ZTF and login email address are the same
        else:
            _ztffp_user_address = _ztffp_email_address
        global _ztffp_user_password
        _ztffp_user_password = os.getenv("ztf_user_password", None)

        # Success!
        return True

    except:
        print("ZTF credentials are not found in the environmental variables.")
        print("Please check the README file on github for setup instructions.")
    
        # Unsuccessful
        return False


def wget_check() -> bool:
    '''Check if wget is installed on the system.'''

    wget_installed: bool = False
    if shutil.which("wget") is None:
        wget_text = (f"wget is not installed on your system "
                     "(not the Python library). "
                     f"Please install wget before continuing.\n")
        print(wget_text)

    else:
        wget_installed = True

    return wget_installed


def random_log_file_name() -> str:
    '''Generate a random log file name.'''

    log_file_name: str | None = None
    while log_file_name is None or os.path.exists(log_file_name):
        random_chars = ''.join(random.choices(string.ascii_uppercase + string.digits,
                                              k=10))
        log_file_name = f"ztffp_{random_chars}.txt"
    
    return log_file_name


def download_ztf_url(url: str, verbose: bool = True) -> str | None:
    '''Download a ZTF files using wget.'''

    # Wget is required to download the ZTF forced photometry request submission
    wget_installed = wget_check()
    if wget_installed==False:
        return None


    wget_command = (f"wget --http-user={_ztfuser} "
                    f"--http-password={_ztfinfo} "
                    f"-O {url.split('/')[-1]} {url}")
    
    if verbose:
        print("Downloading file...")
        print(f'\t{wget_command}')

    subprocess.run(wget_command.split(), capture_output=True)

    return url.split('/')[-1]


def match_ztf_message(job_info, message_body, message_time_epoch, time_delta=10, new_email_matching=False, angular_separation=2):
    '''
    Check if the given email matches the information passed from the log file via job_info.
    If a similar request was passed prior to the submission, you may need to use the 
    only_require_new_email parameter because the information won't exactly match.
    In this case, look only for a close position, the relevant body text, and ensure that the 
    email was sent after the request.
    '''
    
    match = False

    #
    # Only continue if the message was received AFTER the job was submitted
    #
    if message_time_epoch < job_info['cdatetime'].to_list()[0]:

        return match


    message_lines = message_body.splitlines()

    for line in message_lines:

        #
        # Incomplete data product
        #
        if re.search("A request similar to yours is waiting to be processed", line):
            match = False
            break # Stop early if this isn't a finished data product


        if re.search("reqid", line):

            inputs = line.split('(')[-1]

            # Two ways
            # Processing has completed for reqid=XXXX ()
            test_ra = inputs.split('ra=')[-1].split(',')[0]
            test_decl = inputs.split('dec=')[-1].split(')')[0]
            if re.search('minJD', line) and re.search('maxJD', line):
                test_minjd = inputs.split('minJD=')[-1].split(',')[0]
                test_maxjd = inputs.split('maxJD=')[-1].split(',')[0]
            else:
                test_minjd = inputs.split('startJD=')[-1].split(',')[0]
                test_maxjd = inputs.split('endJD=')[-1].split(',')[0]
                
            if new_email_matching:

                # Call this a match only if parameters match
                if np.format_float_positional(float(test_ra), precision=6, pad_right=6).replace(' ','0') == job_info['ra'].to_list()[0] and \
                   np.format_float_positional(float(test_decl), precision=6, pad_right=6).replace(' ','0') == job_info['dec'].to_list()[0] and \
                   (np.format_float_positional(float(test_minjd), precision=6, pad_right=6).replace(' ','0') == job_info['jdstart'].to_list()[0] and \
                   np.format_float_positional(float(test_maxjd), precision=6, pad_right=6).replace(' ','0') == job_info['jdend'].to_list()[0]) or ( \
                       float(test_minjd) - time_delta < float(job_info['jdstart'].to_list()[0]) and float(test_minjd) + time_delta > float(job_info['jdstart'].to_list()[0]) and \
                       float(test_maxjd) - time_delta < float(job_info['jdend'].to_list()[0]) and float(test_maxjd) + time_delta > float(job_info['jdend'].to_list()[0])):

                   match = True
                
            else:

                # Check if new and positions are similar
                submitted_skycoord = SkyCoord(job_info["ra"], job_info["dec"], frame='icrs', unit='deg')
                email_skycoord = SkyCoord(test_ra, test_decl, frame='icrs', unit='deg')
                if submitted_skycoord.separation(email_skycoord).arcsecond < angular_separation and \
                    message_time_epoch > job_info['cdatetime'].to_list()[0]:

                    match = True

    return match


def read_job_log(file_name: str) -> pd.DataFrame:

    job_info = pd.read_html(file_name)[0]
    job_info['ra'] = np.format_float_positional(float(job_info['ra'].to_list()[0]),
                                                precision=6, pad_right=6).replace(' ','0')
    job_info['dec'] = np.format_float_positional(float(job_info['dec'].to_list()[0]),
                                                 precision=6, pad_right=6).replace(' ','0')
    job_info['jdstart'] = np.format_float_positional(float(job_info['jdstart'].to_list()[0]),
                                                     precision=6, pad_right=6).replace(' ','0')
    job_info['jdend'] = np.format_float_positional(float(job_info['jdend'].to_list()[0]),
                                                   precision=6, pad_right=6).replace(' ','0')
    job_info['isostart'] = Time(float(job_info['jdstart'].to_list()[0]),
                                format='jd', scale='utc').iso
    job_info['isoend'] = Time(float(job_info['jdend'].to_list()[0]),
                              format='jd', scale='utc').iso
    job_info['ctime'] = os.path.getctime(file_name) - time.localtime().tm_gmtoff
    job_info['cdatetime'] = datetime.fromtimestamp(os.path.getctime(file_name))

    return job_info


def test_email_connection(n_attempts = 5):
    '''
    Checks the given email address and password
    to see if a connection can be made.
    '''

    # Try a few times to be certain.
    for attempt in range(n_attempts):
    
        try:

            imap = imaplib.IMAP4_SSL(_ztffp_email_server)
            imap.login(_ztffp_email_address, _ztffp_email_password)

            status, messages = imap.select("INBOX")
            if status=='OK':
                found_message = ("Your email inbox was found and contains "
                                 f"{int(messages[0])} messages.\n"
                                 "If this is not correct, please check your settings.")
                print(found_message)
            else:
                print(f"Your inbox was not located. Please check your settings.")
        
            imap.close()
            imap.logout()

            # A successful connection was made
            return True

        # Connection could be broken
        except Exception:
            print("Encountered an exception when connecting to your email address. Trying again.")
            # Give a small timeout in the case of an intermittent connection issue.
            time.sleep(10)

    # No successful connection was made
    return False


def query_ztf_email(log_file_name: str,
                    source_name: str ='temp',
                    new_email_matching: bool = False,
                    verbose: bool = True):
    '''
    Checks the given email address for a message from ZTF.

    Parameters
    ----------
    log_file_name : str
        The name of the log file to check for a match.
    source_name : str, optional
        The name of the source that will be used for output files.
    new_email_matching : bool, optional
        If True, the email must be new.
    verbose : bool, optional
        If True, print out more information for logging.
    '''

    downloaded_file_names = None

    if not os.path.exists(log_file_name):

        print(f"{log_file_name} does not exist.")
        return -1


    # Interpret the request sent to the ZTF forced photometry server
    job_info = read_job_log(log_file_name)


    try:

        imap = imaplib.IMAP4_SSL(_ztffp_email_server)
        imap.login(_ztffp_email_address, _ztffp_email_password)

        status, messages = imap.select("INBOX")

        processing_match = False
        for i in range(int(messages[0]), 0, -1):

            if processing_match:
                break

            # Fetch the email message by ID
            res, msg = imap.fetch(str(i), "(RFC822)")
            for response in msg:
                if isinstance(response, tuple):
                    # Parse a bytes email into a message object
                    msg = email.message_from_bytes(response[1])
                    # decode the email subject
                    sender, encoding = email.header.decode_header(msg.get("From"))[0]
    
                    if not isinstance(sender, bytes) and re.search("ztfpo@ipac\.caltech\.edu", sender):
                                      
    
                        #
                        # Get message body
                        #
                        content_type = msg.get_content_type()
                        body = msg.get_payload(decode=True).decode()
    
                        this_date = msg['Date']
                        this_date_tuple = email.utils.parsedate_tz(msg['Date'])
                        local_date = datetime.fromtimestamp(email.utils.mktime_tz(this_date_tuple))
    
                        
                        #
                        # Check if this is the correct one
                        #
                        if content_type=="text/plain":
                            processing_match = match_ztf_message(job_info, body, local_date, new_email_matching)
                            subject, encoding = email.header.decode_header(msg.get("Subject"))[0]

                            if processing_match:
    
                                # Grab the appropriate URLs
                                lc_url = 'https' + (body.split('_lc.txt')[0] + '_lc.txt').split('https')[-1]
                                log_url = 'https' + (body.split('_log.txt')[0] + '_log.txt').split('https')[-1]
    
    
                                # Download each file
                                lc_initial_file_name = download_ztf_url(lc_url, verbose=verbose)
                                log_initial_file_name = download_ztf_url(log_url, verbose=verbose)    
    
                                # Rename files
                                lc_final_name = f"{source_name.replace(' ','')}_{lc_initial_file_name.split('_')[-1]}"
                                log_final_name = f"{source_name.replace(' ','')}_{log_initial_file_name.split('_')[-1]}"
                                os.rename(lc_initial_file_name, lc_final_name)
                                os.rename(log_initial_file_name, log_final_name)
                                downloaded_file_names = [lc_final_name, log_final_name]

        imap.close()
        imap.logout()

    # Connection could be broken
    except Exception:
        pass

    if downloaded_file_names is not None:

        for file_name in downloaded_file_names:
            if verbose:
                print(f"Downloaded: {file_name}")
    
    return downloaded_file_names


def ztf_forced_photometry(ra: int | float | str | None,
                          decl: int | float | str | None,
                          jdstart: float | None = None,
                          jdend: float | None = None,
                          days: int | float = 60,
                          send: bool = True,
                          verbose: bool = True) -> str | None:
    '''
    Submits a request to the ZTF Forced Photometry service.

    Parameters
    ----------
    ra : int, float, str, or None
        The right ascension of the source in decimal degrees or sexagesimal.
    decl : int, float, str, or None
        The declination of the source in decimal degrees or sexagesimal.
    jdstart : float, optional
        The start Julian date for the query.
        If None, the current date minus 60 days will be used.
    jdend : float, optional
        The end Julian date for the query.
        If None, the current date will be used.
    days : int or float, optional
        The number of days to query.
        This is only used if jdstart and jdend are None.
    send : bool, optional
        If True, the request will be sent to the ZTF Forced Photometry service.
    '''

    # Wget is required for the ZTF forced photometry request submission
    wget_installed = wget_check()
    if wget_installed==False:
        return None

    #
    # Set dates
    #
    if jdend is None:
        jdend = Time(datetime.utcnow(), scale='utc').jd

    if jdstart is None:
        jdstart = jdend - days

    if ra is not None and decl is not None:

        # Check if ra is a decimal
        try:
            # These will trigger the exception if they aren't float
            float(ra)
            float(decl)
            skycoord = SkyCoord(ra, decl, frame='icrs', unit='deg')
        
        # Else assume sexagesimal
        except Exception:
            skycoord = SkyCoord(ra, decl, frame='icrs', unit=(u.hourangle, u.deg))

        # Convert to string to keep same precision.
        # This will make matching easier in the case of submitting multiple jobs.
        jdend_str = np.format_float_positional(float(jdend), precision=6)
        jdstart_str = np.format_float_positional(float(jdstart), precision=6)
        ra_str = np.format_float_positional(float(skycoord.ra.deg), precision=6)
        decl_str = np.format_float_positional(float(skycoord.dec.deg), precision=6)

        log_file_name = random_log_file_name() # Unique file name

        if verbose:
            print(f"Sending ZTF request for (R.A.,Decl)=({ra},{decl})")
        
        wget_command = (f"wget --http-user={_ztfuser} "
                        f"--http-passwd={_ztfinfo} "
                        f"-O {log_file_name} "
                        "https://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi?"
                        f"ra={ra_str}&"
                        f"dec={decl_str}&"
                        f"jdstart={jdstart_str}&"
                        f"jdend={jdend_str}&"
                        f"email={_ztffp_user_address}&userpass={_ztffp_user_password}")

        # Replace .& with .0& to avoid wget error
        wget_command = wget_command.replace('.&', '.0&')

        if verbose:
            print(wget_command)

        if send:
            subprocess.run(wget_command.split(), capture_output=True)

        return log_file_name

    else:
        
        if verbose:
            print("Missing necessary R.A. or declination.")
        return None

def plot_ztf_fp(lc_file_name: str,
                file_format: str = '.png',
                threshold: int | float = 3.0,
                upperlimit: int | float = 5.0,
                verbose: bool = False):
    '''
    Create a simple ZTF forced photometry light curve.
    '''

    # Color mapping for figures
    filter_colors: dict = {'ZTF_g': 'g',
                           'ZTF_r': 'r',
                           'ZTF_i': 'darkorange'}


    try:
        ztf_fp_df = pd.read_csv(lc_file_name, delimiter=' ', comment='#')
    except:
        if verbose:
            print(f"Empty ZTF light curve file ({lc_file_name}). Check the log file.")
        return

    # Rename columns due to mix of , and ' ' separations in the files
    new_cols = {}
    for col in ztf_fp_df.columns:
        new_cols[col] = col.replace(',','')
    
    # Make a cleaned-up version
    ztf_fp_df.rename(columns=new_cols, inplace=True)
    ztf_fp_df.drop(columns=['Unnamed: 0'], inplace=True)

    # Create additional columns with useful calculations
    ztf_fp_df['mjd_midpoint'] = ztf_fp_df['jd'] - 2400000.5 - ztf_fp_df['exptime']/2./86400.
    ztf_fp_df['fp_mag'] = ztf_fp_df['zpdiff'] - 2.5*np.log10(ztf_fp_df['forcediffimflux'])
    ztf_fp_df['fp_mag_unc'] = 1.0857 * ztf_fp_df['forcediffimfluxunc']/ztf_fp_df['forcediffimflux']
    ztf_fp_df['fp_ul'] = ztf_fp_df['zpdiff'] - 2.5*np.log10(upperlimit * ztf_fp_df['forcediffimfluxunc'])


    fig = plt.figure(figsize=(12,6))

    # Iterate over filters
    for ztf_filter in set(ztf_fp_df['filter']):

        filter_df = ztf_fp_df[ztf_fp_df['filter']==ztf_filter]

        # Upper limit df
        ul_filter_df = filter_df[filter_df.forcediffimflux/filter_df.forcediffimfluxunc < threshold]

        # Detections df
        detection_filter_df = filter_df[filter_df.forcediffimflux/filter_df.forcediffimfluxunc >= threshold]

        if verbose:
            print(f"{ztf_filter}: {detection_filter_df.shape[0]} detections and {ul_filter_df.shape[0]} upper limits.")

        # Plot detections
        plt.plot(detection_filter_df.mjd_midpoint, detection_filter_df.fp_mag, color=filter_colors[ztf_filter], marker='o', linestyle='', zorder=3)
        plt.errorbar(detection_filter_df.mjd_midpoint, detection_filter_df.fp_mag, yerr=detection_filter_df.fp_mag_unc, color=filter_colors[ztf_filter], linestyle='', zorder=1)

        # Plot non-detections
        plt.plot(ul_filter_df.mjd_midpoint, ul_filter_df.fp_mag, color=filter_colors[ztf_filter], marker='v', linestyle='', zorder=2)
    

    # Final touches
    plt.gca().invert_yaxis()
    plt.grid(True)
    plt.tick_params(bottom=True, top=True, left=True, right=True, direction='in', labelsize='18', grid_linestyle=':')
    plt.ylabel('ZTF FP (Mag.)', fontsize='20')
    plt.xlabel('Time (MJD)', fontsize='20')
    plt.tight_layout()

    output_file_name = lc_file_name.rsplit('.', 1)[0] + file_format
    fig.savefig(output_file_name)
    plt.close(fig)

    return output_file_name

#
# Wrapper function so that other python code can call this
#
def run_ztf_fp(all_jd: bool = False,
               days: int | float = 60,
               decl: int | float | None = None,
               directory_path: str = '.',
               do_plot: bool = True,
               emailcheck: int = 20,
               fivemindelay: int =60,
               jdend: float | int | None = None,
               jdstart: float | int | None = None,
               logfile: str | None = None,
               mjdend: float | int | None = None,
               mjdstart: float | int | None = None,
               plotfile: str | None = None,
               ra: int | float | None = None,
               skip_clean: bool = False,
               source_name: str = 'temp',
               test_email: bool = False,
               new_email_matching: bool = False,
               verbose: bool = False):
    '''
    Wrapper function to run the ZTF Forced Photometry code.

    Parameters
    ----------
    all_jd : bool, optional
        If True, will run the code for all JDs in the given range. If False, will only run for the first JD in the range. The default is False.
    days : int or float, optional
        Number of days to run the code for. The default is 60.
    decl : int or float, optional
        Declination of the source in degrees. The default is None.
    directory_path : str, optional
        Path to the directory where the code will be run. The default is '.'.
    do_plot : bool, optional
        If True, will create a plot of the light curve. The default is True.
    emailcheck : int, optional
        Number of minutes between email checks. The default is 20.
    fivemindelay : int, optional
        Number of minutes to wait before checking for new data. The default is 60.
    jdend : float or int, optional
        Last JD to run the code for. The default is None.
    jdstart : float or int, optional
        First JD to run the code for. The default is None.
    logfile : str, optional
        Name of the log file. The default is None.
    mjdend : float or int, optional
        Last MJD to run the code for. The default is None.
    mjdstart : float or int, optional
        First MJD to run the code for. The default is None.
    plotfile : str, optional
        Name of the plot file. The default is None.
    ra : int or float, optional
        Right ascension of the source in degrees. The default is None.
    skip_clean : bool, optional
        If True, will skip the cleaning step. The default is False.
    source_name : str, optional
        Name of the source. The default is 'temp'.
    test_email : bool, optional
        If True, will test the email connection. The default is False.
    new_email_matching : bool, optional
        If True, will require the email is new. The default is False.
    verbose : bool, optional
        If True, will print out more information. The default is False.
    '''

    # Stop early if credentials were not found
    credentials_imported = import_credentials()
    if credentials_imported == False:
        return -1

    # Exit early if no sufficient conditions given to run
    run = False
    if (ra is not None and decl is not None) or (logfile is not None) or \
        (plotfile is not None) or (test_email==True):
        run = True

    # Go home early
    if run==False:
        print("Insufficient parameters given to run.")
        return -1

    # Perform an email test if necessary
    if test_email==True:

        # Use comments from the function for now
        email_connection_status = test_email_connection()
        return

    #
    # Change necessary variables based on what was provided
    #

    # Override jd values if mjd arguments are supplied
    if mjdstart is not None:
        jdstart = mjdstart + 2400000.5
    if mjdend is not None:
        jdend = mjdend + 2400000.5

    # Set to full ZTF range
    if all_jd:
        jdstart = 2458194.5
        jdend = Time(datetime.utcnow(), scale='utc').jd

    log_file_name = None
    if logfile is None and plotfile is None:

        log_file_name = ztf_forced_photometry(ra=ra,
                                              decl=decl,
                                              jdstart=jdstart,
                                              jdend=jdend,
                                              days=days)
    
    else:

        log_file_name = logfile
        plot_file_name = plotfile

    if log_file_name is not None:

        # Download via email
        downloaded_file_names = None
        time_start_seconds = time.time()
        while downloaded_file_names is None:

            if time.time() - time_start_seconds < emailcheck:
                if verbose:
                    print(f"Waiting for the email (rechecking every {emailcheck} seconds).")

            downloaded_file_names = query_ztf_email(log_file_name,
                                                    source_name=source_name,
                                                    new_email_matching=new_email_matching,
                                                    verbose=verbose)

            if downloaded_file_names == -1:
                if verbose:
                    print(f"{log_file_name} was not found.")
            elif downloaded_file_names is None:
                if emailcheck < fivemindelay and time.time() - time_start_seconds > 600: # After 5 minutes, change to checking every 1 minute
                    emailcheck = fivemindelay
                    if verbose:
                        print(f"Changing to re-checking every {emailcheck} seconds.")
                time.sleep(emailcheck)


    else:
        downloaded_file_names = [plot_file_name] 
    

    if downloaded_file_names[0] is not None:

        # Open LC file and plot it
        if do_plot:
            figure_file_name = plot_ztf_fp(downloaded_file_names[0], verbose=verbose)
        else:
            figure_file_name = None

    
    #
    # Clean-up
    #
    if skip_clean==False:
        output_directory = f"{directory_path}/{source_name}".replace('//','/')
        # Trim potential extra '/'
        if output_directory[-1:]=='/':
            output_directory = output_directory[:-1]

        # Create directory
        if not os.path.exists(output_directory):
            if verbose:
                print(f"Creating {output_directory}")
            os.makedirs(output_directory)

        # Move all files to this location

        # Wget log file
        output_files = list()
        if log_file_name is not None and os.path.exists(log_file_name):
            shutil.move(log_file_name, f"{output_directory}/{log_file_name.split('/')[-1]}")
            if verbose:
                print(f"{' '*5}ZTF wget log: {output_directory}/{log_file_name.split('/')[-1]}")
            output_files.append(f"{output_directory}/{log_file_name.split('/')[-1]}")


        # Downloaded files
        if isinstance(downloaded_file_names, list):
            for downloaded_file_name in downloaded_file_names:
                if os.path.exists(downloaded_file_name):
                    shutil.move(downloaded_file_name, f"{output_directory}/{downloaded_file_name.split('/')[-1]}")
                    if verbose:
                        print(f"{' '*5}ZTF downloaded file: {output_directory}/{downloaded_file_name.split('/')[-1]}")
                    output_files.append(f"{output_directory}/{downloaded_file_name.split('/')[-1]}")


        # Figure
        if figure_file_name is not None and os.path.exists(figure_file_name):
            shutil.move(figure_file_name, f"{output_directory}/{figure_file_name.split('/')[-1]}")
            if verbose:
                print(f"{' '*5}ZTF figure: {output_directory}/{figure_file_name.split('/')[-1]}")
            output_files.append(f"{output_directory}/{figure_file_name.split('/')[-1]}")

    if len(output_files)==0 or isinstance(output_files, list)==False:
        output_files = None


    # Useful for automation
    return output_files
