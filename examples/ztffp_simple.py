#! python 
#
# Grab ZTF forced photometry on a single field
#
# M.C. Stroh (Northwestern University)
#


import argparse
import sys

import ztffp


def main():

    # First 'fix' possible negative declinations which
    # argparse can't handle on its own
    for i, arg in enumerate(sys.argv):
        if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg

    # Initialize argument parser.
    parser = argparse.ArgumentParser(prog="ztffp_simple.py <ra> <decl>",
                                     description='Grab ZTF forced photometry on the given location.',
                                     formatter_class=argparse.RawTextHelpFormatter)

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
        ztffp.run_ztf_fp(**vars(args), verbose=True)


if __name__ == "__main__":
    main()

