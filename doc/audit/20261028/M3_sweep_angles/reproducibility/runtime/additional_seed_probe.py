"""Compatibility entry point for the improved free-x frequency search."""
import argparse
from frequency_surface import main
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case');p.add_argument('--frequencies',type=float,nargs='+');a=p.parse_args();main(a.case)
