#!/usr/bin/python3

import glob
import os
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='attaching a symlink into each config folder')
    parser.add_argument('floc', help='folder containing all the config files. Need to symlink')
    parser.add_argument('fsrc', help='source directory that contains all root files to symlink as rootFiles')

    args = parser.parse_args()
    floc=args.floc
    fsrc=args.fsrc

    if not os.path.exists(fsrc):
        print("source folder does not exist!")
        exit()

    cfgs=glob.glob(floc+"/**/*cfg",recursive=True)
    cfgDirs=["/".join(cfg.split("/")[:-1]) for cfg in cfgs]

    for cfgDir in cfgDirs:
        #print(f'ln -s {fsrc} {cfgDir+"/rootFiles"}')
        os.system(f'ln -s {fsrc} {cfgDir+"/rootFiles"}')

