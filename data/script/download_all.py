#!/usr/bin/env python
from __future__ import print_function
import urllib2

def main():
    '''
    A quick script that dumps all columns of all the WIDE fields in the CFHTLenS survey
    to the current directory
    '''

    field_list = ['W1m0m0', 'W1m0m1', 'W1m0m2', 'W1m0m3', 'W1m0m4', 'W1m0p1',
                  'W1m0p2', 'W1m0p3', 'W1m1m0', 'W1m1m1', 'W1m1m2', 'W1m1m3',
                  'W1m1m4', 'W1m1p1', 'W1m1p2', 'W1m1p3', 'W1m2m0', 'W1m2m1',
                  'W1m2m2', 'W1m2m3', 'W1m2m4', 'W1m2p1', 'W1m2p2', 'W1m2p3',
                  'W1m3m0', 'W1m3m1', 'W1m3m2', 'W1m3m3', 'W1m3m4', 'W1m3p1',
                  'W1m3p2', 'W1m3p3', 'W1m4m0', 'W1m4m1', 'W1m4m2', 'W1m4m3',
                  'W1m4m4', 'W1m4p1', 'W1m4p2', 'W1m4p3', 'W1p1m0', 'W1p1m1',
                  'W1p1m2', 'W1p1m3', 'W1p1m4', 'W1p1p1', 'W1p1p2', 'W1p1p3',
                  'W1p2m0', 'W1p2m1', 'W1p2m2', 'W1p2m3', 'W1p2m4', 'W1p2p1',
                  'W1p2p2', 'W1p2p3', 'W1p3m0', 'W1p3m1', 'W1p3m2', 'W1p3m3',
                  'W1p3m4', 'W1p3p1', 'W1p3p2', 'W1p3p3', 'W1p4m0', 'W1p4m1',
                  'W1p4m2', 'W1p4m3', 'W1p4m4', 'W1p4p1', 'W1p4p2', 'W1p4p3',
                  'W2m0m0', 'W2m0m1', 'W2m0p1', 'W2m0p2', 'W2m0p3', 'W2m1m0',
                  'W2m1m1', 'W2m1p1', 'W2m1p2', 'W2m1p3', 'W2p1m0', 'W2p1m1',
                  'W2p1p1', 'W2p1p2', 'W2p1p3', 'W2p2m0', 'W2p2m1', 'W2p2p1',
                  'W2p2p2', 'W2p2p3', 'W2p3m0', 'W2p3m1', 'W2p3p1', 'W2p3p2',
                  'W2p3p3', 'W3m0m0', 'W3m0m1', 'W3m0m2', 'W3m0m3', 'W3m0p1',
                  'W3m0p2', 'W3m0p3', 'W3m1m0', 'W3m1m1', 'W3m1m2', 'W3m1m3',
                  'W3m1p1', 'W3m1p2', 'W3m1p3', 'W3m2m0', 'W3m2m1', 'W3m2m2',
                  'W3m2m3', 'W3m2p1', 'W3m2p2', 'W3m2p3', 'W3m3m0', 'W3m3m1',
                  'W3m3m2', 'W3m3m3', 'W3m3p1', 'W3m3p2', 'W3m3p3', 'W3p1m0',
                  'W3p1m1', 'W3p1m2', 'W3p1m3', 'W3p1p1', 'W3p1p2', 'W3p1p3',
                  'W3p2m0', 'W3p2m1', 'W3p2m2', 'W3p2m3', 'W3p2p1', 'W3p2p2',
                  'W3p2p3', 'W3p3m0', 'W3p3m1', 'W3p3m2', 'W3p3m3', 'W3p3p1',
                  'W3p3p2', 'W3p3p3', 'W4m0m0', 'W4m0m1', 'W4m0m2', 'W4m0p1',
                  'W4m1m0', 'W4m1m1', 'W4m1m2', 'W4m1p1', 'W4m1p2', 'W4m1p3',
                  'W4m2m0', 'W4m2p1', 'W4m2p2', 'W4m2p3', 'W4m3m0', 'W4m3p1',
                  'W4m3p2', 'W4m3p3', 'W4p1m0', 'W4p1m1', 'W4p1m2', 'W4p1p1',
                  'W4p2m0', 'W4p2m1', 'W4p2m2']         
    
    for field in field_list:
        filename = field + '.csv'
        
        url = 'http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/tap/sync?REQUEST=doQuery&LANG=ADQL&format=csv&query=SELECT%0D%0A*%0D%0AFROM%0D%0Acfht.clens%0D%0AWHERE%0D%0Acfht.clens.field+%3D+%27' + field + '%27'
    
        response = urllib2.urlopen(url)
    
        print('Downloading', filename)
        
        with open(filename, 'w') as f:
            f.write(response.read())

if __name__ == '__main__':
    main()
