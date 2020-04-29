#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:51:17 2020
@author: Ajit Johnson Nirmal
"""

import pandas as pd

def add_imageid (image,custom_imageid,save_dir):
    print ('loading image' + str(image))
    d = pd.read_csv(image) 
    imid = custom_imageid
    d['ImageId'] = imid
    print ('Image ID is ' + str(custom_imageid))
    # Save
    d.to_csv(save_dir+'/'+str(custom_imageid)+'.csv',index=False)



