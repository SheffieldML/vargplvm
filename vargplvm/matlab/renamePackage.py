#!/usr/bin/python

import os
from tempfile import mkstemp
from shutil import move


def replace(file_path, pattern, subst):
    #Create temp file
    fh, abs_path = mkstemp()
    new_file = open(abs_path,'w')
    old_file = open(file_path)
    for line in old_file:
        new_file.write(line.replace(pattern, subst))
    #close temp file
    new_file.close()
    os.close(fh)
    #fh.close()
    old_file.close()
    #Remove original file
    os.remove(file_path)
    #Move new file
    move(abs_path, file_path)


a = []
for dirname, dirnames, filenames in os.walk('.'):
    # print path to all subdirectories first.
    # for subdirname in dirnames:
    #    print os.path.join(dirname, subdirname)

    # print path to all filenames.
    for filename in filenames:
        if filename.endswith(".m"):
            print os.path.join(dirname, filename)
            a.append(os.path.join(dirname, filename))
print "Found " + str(len(a)) + " source files!"

for ff in a:
    replace(ff, "% SVARGPLVM\n", "% VARGPLVM\n")
#    input_file = file(ff, 'r') # substitute your path to the file
#    lines = input_file.readlines() # list of lines
#    output_file = file('num.txt', 'w') # substitute your path to the file
#    for line in lines:
#        nums = line.split(", ") # assumes numbers separated by ", ", gives list of numbers
#        retain = nums[:5] # keep first 5 numbers
#        new_line = ", ".join(retain)+"\n" # reassemble into a line
#        output_file.write(new_line) # write to file
#        output_file.close()





