# -*- coding: utf-8 -*-
""" 
Script that deploys the static website to www.physics.mcgill.ca/~decotret.

Steps are:
1. Build the website;
2. Then, simply run:

    >>> python deploy.py 

This script requires:
    - Python 3.6+
    - paramiko
    - tqdm
"""
import argparse
import sys
import webbrowser
from contextlib import suppress
from getpass import getpass
from os import listdir, walk
from os.path import getsize, isfile, join

try:
    from paramiko import SSHClient, AutoAddPolicy, AuthenticationException
    from tqdm import tqdm
except ImportError:
    print("paramiko and tqdm are required for this script to run.")
    sys.exit(1)

USERNAME = "tbritt"
LOCAL_DIR = "_site"
REMOTE_DIR = f"/WWW/{USERNAME}/"

DESCRIPTION = """Update the personal website rendered by Hakyll. """


parser = argparse.ArgumentParser(
    description=DESCRIPTION, formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument(
    "-s",
    "--show",
    action="store_true",
    help="Navigate to website with default web browser after deployment",
)


def put_dir(client, source, target, exclude_exts=None):
    """ 
    Upload the contents of the source directory to the target path, including subdirectories. 
    
    Parameters
    ----------
    client : paramiko.SFTPClient
        
    source, target : str or path-like
        Source directory and target directory, respectively.
    
    Yields
    ------
    size : int
        Bytes transferred.
    message : str
        Message specifying the filename, and whether it was transferred or skipped.
    """
    if exclude_exts is None:
        exclude_exts = tuple()

    for item in listdir(source):
        src_path = join(source, item)
        dst_path = item

        if isfile(src_path):
            yield (
                client.put(src_path, dst_path).st_size,
                "transferred: " + str(src_path),
            )

        else:
            with suppress(IOError):
                client.mkdir(dst_path)
            client.chdir(dst_path)
            yield from put_dir(client, src_path, dst_path)
            client.chdir("..")


if __name__ == "__main__":

    arguments = parser.parse_args()

    password = getpass("CPM server password: ")

    with SSHClient() as client:

        client.set_missing_host_key_policy(AutoAddPolicy)
        try:
            client.connect(
                "gollum.physics.mcgill.ca", username=USERNAME, password=password
            )
            print("Connected to CPM server.")
        except AuthenticationException as e:
            print(str(e))
            sys.exit(1)

        # Step 0: Calculate the transfer size
        with client.open_sftp() as sftp_client:
            total_bytes = sum(
                getsize(join(root, file))
                for root, _, files in walk(LOCAL_DIR)
                for file in files
            )

            sftp_client.chdir(REMOTE_DIR)

            # Step 1 : upload content
            upload_stream = put_dir(sftp_client, source=LOCAL_DIR, target=REMOTE_DIR)

            with tqdm(
                desc="Upload to server", unit_scale=True, unit="B", total=total_bytes
            ) as pbar:
                for (bytes_transferred, fname) in upload_stream:
                    pbar.update(bytes_transferred)
                    pbar.write("\r {}".format(fname))

    print("Upload done!")

    if arguments.show:
        print("Opening web page externally...")
        webbrowser.open(f"http://www.physics.mcgill.ca/~{USERNAME}")
