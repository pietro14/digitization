#! /usr/bin/env python3

import boto3
import requests
from boto3sts import credentials as creds
import urllib.parse
import os

def presigned(fun, key, url, bucket, session, verbose):

    session = creds.assumed_session(session, endpoint=url,verify=True)
    s3 = session.client('s3', endpoint_url=url, config=boto3.session.Config(signature_version='s3v4'), verify=True)
    if fun == "get":
        url_out = s3.generate_presigned_url('get_object', 
                                        Params={'Bucket': bucket,
                                                'Key': key}, 
                                        ExpiresIn=3600)
    elif fun == "put":
        url_out = s3.generate_presigned_post(bucket, key, ExpiresIn=3600) #3600
    else:
        url_out = ''
    
    return url_out


if __name__ == "__main__":
    
    #python generate_presigned_urls.py -s <shortname> -i <input_folder> -o <output_folder_in_cygno-sim>
    
    from optparse import OptionParser
    parser = OptionParser(usage='usage: %prog\t [-ubsv] get/put Key')
    parser.add_option('-u','--url', dest='url', type='string', default='https://minio.cloud.infn.it/', 
                      help='url [https://minio.cloud.infn.it/];');
    parser.add_option('-b','--bucket', dest='bucket', type='string', default='cygno-sim', 
                      help='bucket [cygno-sim];');
    parser.add_option('-s','--session', dest='session', type='string', default='infncloud-iam', 
                      help='short name [infncloud-iam];');
    parser.add_option('-v','--verbose', dest='verbose', action="store_true", default=False, help='verbose output;');
    parser.add_option('-i','--input', dest='inputfold', type='string', default='./', help='input folder for digitization');
    parser.add_option('-o','--output', dest='outfold', type='string', default='test', help='tag for output in bucket');
    (options, args) = parser.parse_args()
    
    '''
    url_file = open('presigned_urls.txt','w')
    url_file.write('{\n')
    for infile in os.listdir(options.inputfold):
        if infile.endswith('.root'):
            outname = '{}/digi_{}'.format(options.outfold, infile) #CHANGE BASED ON HOW YOUR FILES ARE NAMED; MINE ARE BASED ON THE NAME OF THE INPUT FILE
            url = presigned('put', outname, options.url, options.bucket, options.session, options.verbose)
            url_file.write("'{}' \t : {},\n".format(infile,str(url)))
            #print(url)
    url_file.write('}')
    url_file.close()
    '''
    url_file = open('presigned_url_bkg.txt','w')
    url_file.write('{\n')
    for infile in os.listdir(options.inputfold):
        if infile.endswith('.root'):
            outname = '{}/digi_{}'.format(options.outfold, infile) #CHANGE BASED ON HOW YOUR FILES ARE NAMED; MINE ARE BASED ON THE NAME OF THE INPUT FILE
            url = presigned('put', outname, options.url, options.bucket, options.session, options.verbose)
            url_file.write("'{}' \t : {},\n".format(infile,str(url)))
            #print(url)
    url_file.write('}')
    url_file.close()