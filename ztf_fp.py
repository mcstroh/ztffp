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
import sys
import time

import matplotlib
matplotlib.use('AGG') # Run faster, comment out for interactive plotting
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore") # We'll get warnings from log10 when there are non-detections


#
# Generic ZTF webdav login
#
_ztfuser = "ztffps"
_ztfinfo = "dontgocrazy!"


#
# Import ZTF email user information
#
def import_credentials():

    try:

        global _ztffp_email_address
        _ztffp_email_address = os.environ["ztf_email_address"]
        global _ztffp_email_password
        _ztffp_email_password = os.environ["ztf_email_password"]
        global _ztffp_email_server
        _ztffp_email_server = os.environ["ztf_email_imapserver"]
        # The email address associated with the ZTF FP server may be an alias, so allow for that possiblility
        global _ztffp_user_address
        if 'ztf_user_address' in os.environ:
            _ztffp_user_address = os.environ["ztf_user_address"]
        # Assume ZTF and login email address are the same
        else:
            _ztffp_user_address = _ztffp_email_address
        global _ztffp_user_password
        _ztffp_user_password = os.environ["ztf_user_password"]

        # Success!
        return True


    except:
        print("ZTF credentials are not found in the environmental variables.\nPlease check the README file on github for information on how to set this up.")
    
        # Unsuccessful
        return False



def wget_check():

    if shutil.which("wget") is None:
       print("ztf_fp.py requires wget installed on your system (not the wget python library). Please install before continuing.")
       return False
    else:
        return True


def random_log_file_name():

    log_file_name = None
    while log_file_name is None or os.path.exists(log_file_name):
        log_file_name = f"ztffp_{''.join(random.choices(string.ascii_uppercase + string.digits, k=10))}.txt"
    
    return log_file_name



def download_ztf_url(url, verbose=True):

    # Wget is required to download the ZTF forced photometry request submission
    wget_installed = wget_check()
    if wget_installed==False:
        return None


    wget_command = f"wget --http-user={_ztfuser} --http-password={_ztfinfo} -O {url.split('/')[-1]} \"{url}\""
    
    if verbose:
        print("Downloading file...")
        print('\t' + wget_command)
    p = subprocess.Popen(wget_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = p.communicate()

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



def read_job_log(file_name):

    job_info = pd.read_html(file_name)[0]
    job_info['ra'] = np.format_float_positional(float(job_info['ra'].to_list()[0]), precision=6, pad_right=6).replace(' ','0')
    job_info['dec'] = np.format_float_positional(float(job_info['dec'].to_list()[0]), precision=6, pad_right=6).replace(' ','0')
    job_info['jdstart'] = np.format_float_positional(float(job_info['jdstart'].to_list()[0]), precision=6, pad_right=6).replace(' ','0')
    job_info['jdend'] = np.format_float_positional(float(job_info['jdend'].to_list()[0]), precision=6, pad_right=6).replace(' ','0')
    job_info['isostart'] = Time(float(job_info['jdstart'].to_list()[0]), format='jd', scale='utc').iso
    job_info['isoend'] = Time(float(job_info['jdend'].to_list()[0]), format='jd', scale='utc').iso
    job_info['ctime'] = os.path.getctime(file_name) - time.localtime().tm_gmtoff
    job_info['cdatetime'] = datetime.fromtimestamp(os.path.getctime(file_name))

    return job_info



#
# Test the connection to the users email server
#
def test_email_connection(n_attempts = 5):

    # Try a few times to be certain.
    for attempt in range(n_attempts):
    
        try:

            imap = imaplib.IMAP4_SSL(_ztffp_email_server)
            imap.login(_ztffp_email_address, _ztffp_email_password)

            status, messages = imap.select("INBOX")
            if status=='OK':
                print(f"Your email inbox was found and contains {int(messages[0])} messages.\nIf this is not correct, please check your settings.")
            else:
                print(f"Your inbox was not located. Please check your settings.")
        
            imap.close()
            imap.logout()

            # A successful connection was made
            return True


        # Connection could be broken
        except Exception:
            print("Encountered an exception when connecting to your email address. Trying again.")
            time.sleep(10) # Give a small timeout in the case of an intermittent connection issue.
            

    # No successful connection was made
    return False



#
# Look for email to download data products
# 
def query_ztf_email(log_file_name, source_name='temp', new_email_matching=False, verbose=True):

    downloaded_file_names = None

    if not os.path.exists(log_file_name):

        f"{log_file_name} does not exist."
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

                    


def ztf_forced_photometry(ra, decl, jdstart=None, jdend=None, days=60, send=True, verbose=True):

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
            

        # Convert to string to keep same precision. This will make matching easier in the case of submitting multiple jobs.
        jdend_str = np.format_float_positional(float(jdend), precision=6)
        jdstart_str = np.format_float_positional(float(jdstart), precision=6)
        ra_str = np.format_float_positional(float(skycoord.ra.deg), precision=6)
        decl_str = np.format_float_positional(float(skycoord.dec.deg), precision=6)


        log_file_name = random_log_file_name() # Unique file name

        if verbose:
            print(f"Sending ZTF request for (R.A.,Decl)=({ra},{decl})")
        
        wget_command = f"wget --http-user={_ztfuser} --http-passwd={_ztfinfo} -O {log_file_name} \"https://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi?" + \
                       f"ra={ra_str}&" + \
                       f"dec={decl_str}&" + \
                       f"jdstart={jdstart_str}&" +\
                       f"jdend={jdend_str}&" + \
                       f"email={_ztffp_user_address}&userpass={_ztffp_user_password}\""
        
        if verbose:
            print(wget_command)

        if send:
    
            p = subprocess.Popen(wget_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            stdout, stderr = p.communicate()

            if verbose:
                print(stdout.decode('utf-8'))

        
        return log_file_name


    else:
        
        if verbose:
            print("Missing necessary R.A. or declination.")
        return None



def plot_ztf_fp(lc_file_name, file_format='.png', threshold=3.0, upperlimit=5.0, verbose=False):

    # Color mapping for figures
    filter_colors = {'ZTF_g': 'g', 'ZTF_r': 'r', 'ZTF_i': 'darkorange'}


    # Use default naming convention
    plot_file_name = lc_file_name.replace('.txt', '.png')

    try:
        ztf_fp_df = pd.read_csv(lc_file_name, delimiter=' ', comment='#')
    except:
        if verbose:
            print(f"Empty ZTF light curve file ({lc_file_name}). Check the log file.")
        return

    # Rename columns due to mix of , and ' ' separations in the files
    new_cols = dict()
    for col in ztf_fp_df.columns:
        new_cols[col] = col.replace(',','')
    
    # Make a clean version
    ztf_fp_df.rename(columns=new_cols, inplace=True)
    ztf_fp_df.drop(columns=['Unnamed: 0'], inplace=True)

    #
    # Create additional columns with useful calculations
    #
    ztf_fp_df['mjd_midpoint'] = ztf_fp_df['jd'] - 2400000.5 - ztf_fp_df['exptime']/2./86400. # Do it all at once here
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
def run_ztf_fp(all_jd=False, days=60, decl=None, directory_path='.',
               do_plot=True, emailcheck=20, fivemindelay=60, jdend=None, 
               jdstart=None, logfile=None, mjdend=None, mjdstart=None, 
               plotfile=None, ra=None, skip_clean=False, source_name='temp', 
               test_email=False, new_email_matching=False, verbose=False):


    #
    # Stop early if credentials were not found
    #
    credentials_imported = import_credentials()
    if credentials_imported == False:
        return -1


    #
    # Exit early if no sufficient conditions given to run
    #
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

        log_file_name = ztf_forced_photometry(ra=ra, decl=decl, jdstart=jdstart, jdend=jdend, days=days)
    
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

            downloaded_file_names = query_ztf_email(log_file_name, source_name=source_name, new_email_matching=new_email_matching, verbose=verbose)

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


        #
        # Move all files to this location
        #
        
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



def main():

    # First 'fix' possible negative declinations which argparse can't handle on its own
    for i, arg in enumerate(sys.argv):
        if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg

    # Initialize argument parser.
    parser = argparse.ArgumentParser(prog="ztf_fp.py <ra> <decl>", description='Grab ZTF forced photometry on the given location.', formatter_class=argparse.RawTextHelpFormatter)

    # Necessary arguments
    parser.add_argument('ra', type=str, nargs='?', default=None,
                   help="Right ascension of the target. Can be provided in DDD.ddd or HH:MM:SS.ss formats.")
    parser.add_argument('decl', type=str, nargs='?', default=None,
                   help="Declination of the target. Can be provided in +/-DDD.ddd or DD:MM:SS.ss formats.")
    # OR (query has already been sent)
    parser.add_argument('-logfile', metavar='logfile', type=str, nargs='?', default=None, 
                   help='Log file to process instead of submitting a new job.')
    # OR (option to only build light curve)
    parser.add_argument('-plotfile', metavar='plotfile', type=str, nargs='?', default=None, 
                   help='Light curve file to plot instead of submitting and downloading a new job.')
    # OR (option to test email connection without sending a request to ZTF)
    parser.add_argument('-emailtest', action='store_true', dest='test_email',
                   help='Test your email settings. This is performed without sending a request to the ZTF server.')
    parser.set_defaults(test_email=False)


    # Additional arguments we can use
    parser.add_argument('-source_name', metavar='source_name', type=str, nargs='?', default='temp_source', 
                   help='Source name that will be used to name output files.')    
    parser.add_argument('-mjdstart', metavar='mjdstart', type=float, nargs='?', default=None, 
                   help='Start of date range for forced photometry query. Overrides -jdstart.')    
    parser.add_argument('-mjdend', metavar='mjdend', type=float, nargs='?', default=None, 
                   help='End of date range for forced photometry query. Overrides -jdstop')
    parser.add_argument('-jdstart', metavar='jdstart', type=float, nargs='?', default=None, 
                   help='Start of date range for forced photometry query.')    
    parser.add_argument('-jdend', metavar='jdend', type=float, nargs='?', default=None, 
                   help='End of date range for forced photometry query.')
    parser.add_argument('-ztf_all_jd', action='store_true', dest='all_jd',
                   help='Use the full range of ZTF public dates.')
    parser.set_defaults(all_jd=False)
    parser.add_argument('-days', metavar='days', type=int, nargs='?', default=60, 
                   help='Number of days prior to jdend to query (or number of days prior to today if jdend is not given).')    
    parser.add_argument('-emailcheck', metavar='emailcheck', type=float, nargs='?', default=20, 
                   help='How often to recheck your email for the ZTF results.')
    parser.add_argument('-skip_clean', action='store_false', dest='skip_clean',
                   help='After completion skip placing all output files in the same directory.')
    parser.set_defaults(skip_clean=False)
    parser.add_argument('-directory_path', metavar='directory_path', type=str, nargs='?', default='.', 
                   help='Path to directory for clean-up. Requires -directory option.') 
    parser.add_argument('-fivemindelay', metavar='fivemindelay', type=float, nargs='?', default=60, 
                   help='How often (in seconds) to query the email after 5 minutes have elapsed.')    
    parser.add_argument('-email_matching', action='store_false', dest='new_email_matching',
                   help='Use the alternative constraints requiring a new email, and similar positions, but no time matching.')
    parser.set_defaults(new_email_matching=True)
    parser.add_argument('-skip_plot', action='store_false', dest='do_plot',
                   help='Skip making the plot. Useful for automated and batch use cases, or if user wants to use their personal plotting code.')
    parser.set_defaults(do_plot=True)



    try:
        # Validate inputs
        args = parser.parse_args()
        run = True

    except Exception:

        run = False
    

    # Don't go further if there were problems with arguments or inputs
    if run:

        run_ztf_fp(**vars(args), verbose=True)



if __name__ == "__main__":
    main()

