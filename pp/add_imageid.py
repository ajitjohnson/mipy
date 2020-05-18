#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:51:17 2020
@author: Ajit Johnson Nirmal
"""

import pandas as pd

def add_imageid (image,custom_imageid,save_dir):
    """
    

    Parameters
    ----------
    image : list
        List of Image locations.
    custom_imageid : string
        Image ID to add to the image.
    save_dir : string
        Location to the directory to save.

    Returns
    -------
    None.
    
    Example
    -------
    image_loc = glob.glob("/Users/aj/Desktop/ptcl_tma/m_quat/Angioimmunoblastic Cell Lymphoma/*.csv")
    save_dir = '/Users/aj/Desktop/ptcl_tma/m_quat/' 
    for i in image_loc:
        custom_imageid = i.split("/")[-1:][0].split(".")[0]
        add_imageid(image=i, custom_imageid=custom_imageid, save_dir=save_dir)

    """
    # USe this for running in loop
    
    print ('loading image' + str(image))
    d = pd.read_csv(image) 
    imid = custom_imageid
    d['ImageId'] = imid
    print ('Image ID is ' + str(custom_imageid))
    # Save
    d.to_csv(save_dir+'/'+str(custom_imageid)+'.csv',index=False)



