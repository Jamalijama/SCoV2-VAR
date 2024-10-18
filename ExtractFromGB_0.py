import os
from ReadWriteFile import *

path = './data/not_aligned/'

file_list = os.listdir(path)
for file_obj in file_list:
    if file_obj.endswith('.gb'):
        ReadGbAndWrite(path, file_obj)
